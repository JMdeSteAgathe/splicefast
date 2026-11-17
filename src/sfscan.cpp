// ============================================================================
//  sfscan : lecteur BAM/CRAM de SpliceFast, en C++ / htslib
//
//  Une seule lecture du fichier produit, au choix ou ensemble :
//    --count  : jonctions (operations N) + reads NON coupes a chaque site
//               d'epissage annote                     (etape 'count')
//    --cov    : profondeur moyenne du corps de chaque intron et des bordures
//               exoniques                              (module 'retention')
//
//  Les deux calculs reproduisent EXACTEMENT les versions Python
//  (bam.scan_chunk et coverage.scan_chunk_cov) : memes filtres, memes
//  conventions de coordonnees globales, memes arrondis.
//
//  D'ou vient la vitesse
//  ---------------------
//  * decompression BGZF multi-thread par htslib (le vrai goulot en Python) ;
//  * l'idee de PanDepth -- un compteur par base -- appliquee a la grille de
//    segments : chaque bloc aligne coute deux recherches dans un tableau trie,
//    faites par "galop" a partir de la position du read precedent. Un BAM
//    etant trie, la reponse est presque toujours a 0-2 cases : acces memoire
//    chauds, sans jamais allouer de tableau de profondeur par chromosome
//    (1 Go pour chr1 chez PanDepth).
//  * table de hachage a adressage ouvert pour les jonctions.
//
//  Echanges avec Python : fichiers .npy (lus et ecrits ici sans dependance).
// ============================================================================
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

static const char *VERSION = "sfscan 2.1.0";

enum ExitCode { OK = 0, E_USAGE = 2, E_HEADER = 3, E_READ = 4, E_WRITE = 5, E_INPUT = 6, E_CORRUPT = 7 };

struct Fail : std::runtime_error {
    int code;
    Fail(int c, const std::string &m) : std::runtime_error(m), code(c) {}
};

// ---------------------------------------------------------------------------
// .npy minimal (little-endian, C order)
// ---------------------------------------------------------------------------
struct Npy {
    std::string descr;
    std::vector<size_t> shape;
    std::vector<char> raw;
    size_t count() const {
        size_t n = 1;
        for (size_t s : shape) n *= s;
        return n;
    }
};

static std::string dict_value(const std::string &h, const std::string &key) {
    size_t p = h.find("'" + key + "'");
    if (p == std::string::npos) return "";
    p = h.find(':', p);
    if (p == std::string::npos) return "";
    ++p;
    while (p < h.size() && h[p] == ' ') ++p;
    if (p >= h.size()) return "";
    if (h[p] == '\'') {
        size_t e = h.find('\'', p + 1);
        return h.substr(p + 1, e - p - 1);
    }
    if (h[p] == '(') {
        size_t e = h.find(')', p);
        return h.substr(p, e - p + 1);
    }
    size_t e = h.find_first_of(",}", p);
    return h.substr(p, e - p);
}

static Npy npy_load(const std::string &path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw Fail(E_INPUT, "impossible d'ouvrir " + path);
    char magic[8];
    f.read(magic, 8);
    if (!f || std::memcmp(magic, "\x93NUMPY", 6) != 0)
        throw Fail(E_INPUT, path + " : pas un fichier .npy");
    uint32_t hlen = 0;
    if (magic[6] == 1) {
        unsigned char b[2];
        f.read((char *)b, 2);
        hlen = b[0] | (b[1] << 8);
    } else {
        unsigned char b[4];
        f.read((char *)b, 4);
        hlen = b[0] | (b[1] << 8) | (b[2] << 16) | ((uint32_t)b[3] << 24);
    }
    std::string h(hlen, ' ');
    f.read(&h[0], hlen);
    Npy a;
    a.descr = dict_value(h, "descr");
    if (dict_value(h, "fortran_order").find("True") != std::string::npos)
        throw Fail(E_INPUT, path + " : ordre Fortran non gere");
    std::string sh = dict_value(h, "shape");
    size_t i = 0;
    while (i < sh.size()) {
        if (sh[i] >= '0' && sh[i] <= '9') {
            size_t j = i;
            while (j < sh.size() && sh[j] >= '0' && sh[j] <= '9') ++j;
            a.shape.push_back(std::stoull(sh.substr(i, j - i)));
            i = j;
        } else {
            ++i;
        }
    }
    size_t item = (a.descr == "<i8" || a.descr == "<f8") ? 8 : 4;
    a.raw.resize(a.count() * item);
    f.read(a.raw.data(), (std::streamsize)a.raw.size());
    if ((size_t)f.gcount() != a.raw.size()) throw Fail(E_INPUT, path + " : tronque");
    return a;
}

static std::vector<int64_t> npy_i64(const std::string &path, std::vector<size_t> *shape = nullptr) {
    Npy a = npy_load(path);
    if (a.descr != "<i8") throw Fail(E_INPUT, path + " : int64 attendu, trouve " + a.descr);
    std::vector<int64_t> v(a.count());
    std::memcpy(v.data(), a.raw.data(), a.raw.size());
    if (shape) *shape = a.shape;
    return v;
}

static void npy_save(const std::string &path, const char *descr, const std::vector<size_t> &shape,
                     const void *data, size_t nbytes) {
    std::string sh = "(";
    for (size_t i = 0; i < shape.size(); ++i) sh += std::to_string(shape[i]) + (shape.size() == 1 ? "," : (i + 1 < shape.size() ? ", " : ""));
    sh += ")";
    std::string h = std::string("{'descr': '") + descr + "', 'fortran_order': False, 'shape': " + sh + ", }";
    size_t total = 10 + h.size() + 1;
    h.append((64 - total % 64) % 64, ' ');
    h.push_back('\n');
    std::string tmp = path + ".part";
    FILE *f = std::fopen(tmp.c_str(), "wb");
    if (!f) throw Fail(E_WRITE, "ecriture impossible : " + tmp);
    unsigned char head[10] = {0x93, 'N', 'U', 'M', 'P', 'Y', 1, 0,
                              (unsigned char)(h.size() & 0xFF), (unsigned char)(h.size() >> 8)};
    bool ok = std::fwrite(head, 1, 10, f) == 10 && std::fwrite(h.data(), 1, h.size(), f) == h.size() &&
              (nbytes == 0 || std::fwrite(data, 1, nbytes, f) == nbytes);
    ok = (std::fflush(f) == 0) && ok;
    ok = (std::fclose(f) == 0) && ok;
    // Windows : rename() refuse d'ecraser un fichier existant (reliquat d'un
    // run interrompu) -> on le supprime d'abord.
    if (ok) std::remove(path.c_str());
    if (!ok || std::rename(tmp.c_str(), path.c_str()) != 0) throw Fail(E_WRITE, "ecriture impossible : " + path);
}

// ---------------------------------------------------------------------------
// Recherche par galop a partir d'un indice de depart.
//   UPPER=false : premier i tel que v[i] >= x   (lower_bound)
//   UPPER=true  : premier i tel que v[i] >  x   (upper_bound)
// ---------------------------------------------------------------------------
template <bool UPPER>
static inline size_t gallop(const int64_t *v, size_t n, int64_t x, size_t h) {
#define BEFORE(i) (UPPER ? (v[i] <= x) : (v[i] < x))
    if (h > n) h = n;
    size_t lo, hi;
    if (h < n && BEFORE(h)) {  // la reponse est apres h
        lo = h + 1;
        size_t step = 1, cand = h + 1;
        while (cand < n && BEFORE(cand)) {
            lo = cand + 1;
            step <<= 1;
            cand = h + step;
        }
        hi = cand < n ? cand : n;
    } else {  // la reponse est <= h
        hi = h;
        size_t step = 1;
        for (;;) {
            if (hi == 0) return 0;
            size_t cand = hi >= step ? hi - step : 0;
            if (BEFORE(cand)) {
                lo = cand + 1;
                break;
            }
            hi = cand;
            if (cand == 0) return 0;
            step <<= 1;
        }
    }
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (BEFORE(mid))
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
#undef BEFORE
}

// ---------------------------------------------------------------------------
// Table de hachage int64 -> int64 (adressage ouvert, sondage lineaire)
// ---------------------------------------------------------------------------
class JunctionMap {
    std::vector<int64_t> keys_, vals_;
    size_t mask_ = 0, size_ = 0;
    static inline uint64_t mix(uint64_t z) {
        z += 0x9E3779B97F4A7C15ULL;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    void rehash(size_t cap) {
        std::vector<int64_t> ok(cap, -1), ov(cap, 0);
        size_t m = cap - 1;
        for (size_t i = 0; i < keys_.size(); ++i) {
            if (keys_[i] < 0) continue;
            size_t p = mix((uint64_t)keys_[i]) & m;
            while (ok[p] >= 0) p = (p + 1) & m;
            ok[p] = keys_[i];
            ov[p] = vals_[i];
        }
        keys_.swap(ok);
        vals_.swap(ov);
        mask_ = m;
    }

  public:
    JunctionMap() { rehash(1 << 16); }
    inline void inc(int64_t k) {
        size_t p = mix((uint64_t)k) & mask_;
        for (;;) {
            int64_t cur = keys_[p];
            if (cur == k) {
                ++vals_[p];
                return;
            }
            if (cur < 0) {
                keys_[p] = k;
                vals_[p] = 1;
                if (++size_ * 10 > keys_.size() * 7) rehash(keys_.size() * 2);
                return;
            }
            p = (p + 1) & mask_;
        }
    }
    void sorted(std::vector<int64_t> &k, std::vector<int64_t> &c) const {
        std::vector<std::pair<int64_t, int64_t>> kv;
        kv.reserve(size_);
        for (size_t i = 0; i < keys_.size(); ++i)
            if (keys_[i] >= 0) kv.emplace_back(keys_[i], vals_[i]);
        std::sort(kv.begin(), kv.end());
        k.resize(kv.size());
        c.resize(kv.size());
        for (size_t i = 0; i < kv.size(); ++i) {
            k[i] = kv[i].first;
            c[i] = kv[i].second;
        }
    }
    size_t size() const { return size_; }
};

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------
struct Opts {
    std::string bam, reference, ref_names, ref_off;
    int threads = 1;
    int64_t max_reads = 0;
    // comptage
    std::string count_prefix, sites;
    int64_t anchor = 5, min_mapq = 0, exclude = 0x4 | 0x100 | 0x200 | 0x400 | 0x800;
    int64_t min_intron = 20, max_intron = 1000000, min_anchor_split = 5;
    // couverture
    std::string cov_out, cov_grid, cov_idx, cov_free;
    int64_t cov_mapq = 0, cov_exclude = 0x4 | 0x100 | 0x200 | 0x400 | 0x800;
};

static void usage() {
    std::fprintf(stderr,
        "%s -- lecteur BAM/CRAM de SpliceFast (appele par le module Python)\n\n"
        "  sfscan --bam FICHIER --ref-off ref_off.npy [--ref-names noms.txt] [-@ threads]\n"
        "         [--count PREFIXE [--sites sites.npy] [--anchor 5] [--min-mapq 0]\n"
        "                  [--exclude 3844] [--min-intron 20] [--max-intron 1000000]\n"
        "                  [--min-anchor-split 5]]\n"
        "         [--cov SORTIE.npy --cov-grid q.npy --cov-free free.npy --cov-idx idx.npy\n"
        "                  [--cov-mapq 0] [--cov-exclude 3844]]\n"
        "         [--reference genome.fa]   (CRAM)\n"
        "         [--max-reads N]           (mesure de debit)\n\n"
        "Sorties --count : PREFIXE.jkey.npy PREFIXE.jcount.npy PREFIXE.nonsplit.npy\n"
        "Sortie  --cov   : float32 (introns x 3) = [profondeur intron, base bordante, fenetre bordante]\n"
        "Bilan JSON sur la sortie standard. Codes : 2 usage, 3 entete, 4 lecture, 5 ecriture, 6 entree, 7 fichier tronque.\n",
        VERSION);
}

static Opts parse(int argc, char **argv) {
    Opts o;
    auto need = [&](int &i) -> std::string {
        if (i + 1 >= argc) throw Fail(E_USAGE, std::string("valeur manquante apres ") + argv[i]);
        return argv[++i];
    };
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help") { usage(); std::exit(0); }
        else if (a == "--version") { std::printf("%s (htslib %s)\n", VERSION, hts_version()); std::exit(0); }
        else if (a == "--bam") o.bam = need(i);
        else if (a == "--reference") o.reference = need(i);
        else if (a == "--ref-names") o.ref_names = need(i);
        else if (a == "--ref-off") o.ref_off = need(i);
        else if (a == "-@" || a == "--threads") o.threads = std::max(1, std::atoi(need(i).c_str()));
        else if (a == "--max-reads") o.max_reads = std::atoll(need(i).c_str());
        else if (a == "--count") o.count_prefix = need(i);
        else if (a == "--sites") o.sites = need(i);
        else if (a == "--anchor") o.anchor = std::atoll(need(i).c_str());
        else if (a == "--min-mapq") o.min_mapq = std::atoll(need(i).c_str());
        else if (a == "--exclude") o.exclude = std::strtoll(need(i).c_str(), nullptr, 0);
        else if (a == "--min-intron") o.min_intron = std::atoll(need(i).c_str());
        else if (a == "--max-intron") o.max_intron = std::atoll(need(i).c_str());
        else if (a == "--min-anchor-split") o.min_anchor_split = std::atoll(need(i).c_str());
        else if (a == "--cov") o.cov_out = need(i);
        else if (a == "--cov-grid") o.cov_grid = need(i);
        else if (a == "--cov-free") o.cov_free = need(i);
        else if (a == "--cov-idx") o.cov_idx = need(i);
        else if (a == "--cov-mapq") o.cov_mapq = std::atoll(need(i).c_str());
        else if (a == "--cov-exclude") o.cov_exclude = std::strtoll(need(i).c_str(), nullptr, 0);
        else throw Fail(E_USAGE, "option inconnue : " + a);
    }
    if (o.bam.empty() || o.ref_off.empty()) throw Fail(E_USAGE, "--bam et --ref-off sont obligatoires");
    if (o.count_prefix.empty() && o.cov_out.empty() && !o.max_reads)
        throw Fail(E_USAGE, "rien a faire : indiquer --count et/ou --cov");
    if (!o.cov_out.empty() && (o.cov_grid.empty() || o.cov_idx.empty() || o.cov_free.empty()))
        throw Fail(E_USAGE, "--cov exige --cov-grid, --cov-free et --cov-idx");
    return o;
}

// ---------------------------------------------------------------------------
int main(int argc, char **argv) {
    auto t0 = std::chrono::steady_clock::now();
    samFile *fp = nullptr;
    sam_hdr_t *hdr = nullptr;
    bam1_t *b = nullptr;
    int rc = OK;
    try {
        Opts o = parse(argc, argv);
        const bool do_count = !o.count_prefix.empty();
        const bool do_cov = !o.cov_out.empty();

        std::vector<int64_t> refoff = npy_i64(o.ref_off);

        // --- comptage : sites
        std::vector<int64_t> sites;
        if (do_count && !o.sites.empty()) sites = npy_i64(o.sites);
        const size_t nsite = sites.size();
        const bool have_sites = do_count && nsite > 0;
        std::vector<int32_t> sdiff(have_sites ? nsite + 2 : 0, 0);
        JunctionMap junc;

        // --- couverture : grille de segments
        std::vector<int64_t> q, seg_free, idx;
        std::vector<size_t> idx_shape;
        size_t nq = 0, nseg = 0;
        std::vector<int64_t> part;
        std::vector<int32_t> full;
        if (do_cov) {
            q = npy_i64(o.cov_grid);
            seg_free = npy_i64(o.cov_free);
            idx = npy_i64(o.cov_idx, &idx_shape);
            nq = q.size();
            nseg = nq ? nq - 1 : 0;
            if (seg_free.size() != nseg) throw Fail(E_INPUT, "cov-free : taille incoherente avec la grille");
            if (idx_shape.size() != 2 || idx_shape[1] != 8) throw Fail(E_INPUT, "cov-idx : forme (M, 8) attendue");
            for (int64_t v : idx)
                if (v < 0 || (size_t)v > nseg) throw Fail(E_INPUT, "cov-idx : indice hors grille");
            part.assign(nseg + 1, 0);
            full.assign(nseg + 2, 0);
        }

        // --- ouverture
        fp = sam_open(o.bam.c_str(), "r");
        if (!fp) throw Fail(E_READ, "ouverture impossible : " + o.bam);
        if (o.threads > 1 && hts_set_threads(fp, o.threads) != 0)
            std::fprintf(stderr, "sfscan: avertissement, threads non actives\n");
        if (!o.reference.empty()) hts_set_fai_filename(fp, o.reference.c_str());
        hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS,
                    SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
        hts_set_opt(fp, CRAM_OPT_DECODE_MD, 0);
        // Fichier tronque : inutile d'attendre, ce n'est pas une deconnexion.
        if (fp->format.compression == bgzf && fp->fp.bgzf) {
            int eof = bgzf_check_EOF(fp->fp.bgzf);
            if (eof == 0)
                throw Fail(E_CORRUPT, "fichier tronque (marqueur de fin BGZF absent) : " + o.bam +
                                          " -- copie incomplete ou alignement interrompu");
        }
        hdr = sam_hdr_read(fp);
        if (!hdr) throw Fail(E_READ, "entete illisible : " + o.bam);
        // --- correspondance des references BAM -> projet, PAR NOM.
        // Deux runs alignes avec des index differents (25 chromosomes d'un
        // cote, 194 sequences avec contigs non places / alt de l'autre) sont
        // compatibles tant que les chromosomes communs ont la meme longueur :
        // les reads des sequences absentes du projet sont simplement ignores
        // (l'annotation n'y a aucun gene). Une longueur differente = autre
        // assemblage (hg19 / hg38) : refus.
        const int nref = sam_hdr_nref(hdr);
        const size_t nproj = refoff.size() - 1;
        std::vector<int64_t> tid_off(nref, -1);
        int n_mapped = 0;
        if (!o.ref_names.empty()) {
            std::vector<std::string> pname;
            std::vector<int64_t> plen;
            std::ifstream nf(o.ref_names);
            std::string line;
            while (std::getline(nf, line)) {
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty()) continue;
                size_t t = line.find('\t');
                pname.push_back(line.substr(0, t));
                plen.push_back(t == std::string::npos ? -1 : std::atoll(line.c_str() + t + 1));
            }
            if (pname.size() != nproj)
                throw Fail(E_INPUT, "ref-names : " + std::to_string(pname.size()) + " noms pour " +
                                        std::to_string(nproj) + " offsets");
            auto find = [&](const std::string &n) -> int {
                std::vector<std::string> cands = {n, "chr" + n};
                if (n.compare(0, 3, "chr") == 0) cands.push_back(n.substr(3));
                if (n == "MT" || n == "M") cands.push_back("chrM");
                if (n == "chrM") cands.push_back("MT");
                for (const auto &c : cands)
                    for (size_t j = 0; j < pname.size(); ++j)
                        if (pname[j] == c) return (int)j;
                return -1;
            };
            std::string ignored;
            int n_ignored = 0;
            for (int i = 0; i < nref; ++i) {
                std::string n = sam_hdr_tid2name(hdr, i);
                int j = find(n);
                if (j < 0) {
                    if (n_ignored++ < 3) ignored += (ignored.empty() ? "" : ", ") + n;
                    continue;
                }
                int64_t L = (int64_t)sam_hdr_tid2len(hdr, i);
                if (plen[j] >= 0 && plen[j] != L)
                    throw Fail(E_HEADER, o.bam + " : " + n + " mesure " + std::to_string(L) +
                                             " pb mais " + std::to_string(plen[j]) +
                                             " dans le projet (autre assemblage du genome ?)");
                tid_off[i] = refoff[j];
                ++n_mapped;
            }
            if (n_mapped == 0)
                throw Fail(E_HEADER, o.bam + " : aucun chromosome en commun avec le projet "
                                             "(BAM aligne sur un autre genome ?)");
            if (n_ignored)
                std::fprintf(stderr, "sfscan: info, %d sequence(s) absente(s) du projet ignoree(s) (%s%s)\n",
                             n_ignored, ignored.c_str(), n_ignored > 3 ? ", ..." : "");
        } else {  // ancien mode : entete identique exige
            if ((size_t)nref != nproj)
                throw Fail(E_HEADER, o.bam + " : " + std::to_string(nref) + " references dans l'entete, " +
                                         std::to_string(nproj) + " dans le projet");
            for (int i = 0; i < nref; ++i) tid_off[i] = refoff[i];
            n_mapped = nref;
        }
        int64_t n_other_contig = 0;

        b = bam_init1();
        int64_t n_reads = 0, n_spliced = 0, n_cov_reads = 0, n_seen = 0;
        size_t hint_sites = 0, hint_seg = 0;
        const int64_t anchor2 = 2 * o.anchor;
        const int64_t LEN_MASK = (1LL << 24) - 1;
        const int64_t *S = have_sites ? sites.data() : nullptr;
        const int64_t *Q = do_cov ? q.data() : nullptr;

        auto add_site = [&](int64_t a, int64_t bb, size_t &hint) {
            size_t i0 = gallop<false>(S, nsite, a, hint);
            size_t i1 = gallop<true>(S, nsite, bb, i0);
            hint = i0;
            if (i1 > i0) {
                ++sdiff[i0];
                --sdiff[i1];
            }
        };
        // bases couvertes par [a, e) reparties sur les segments de la grille
        auto add_cov = [&](int64_t a, int64_t e, size_t &hint) {
            // ia = dernier indice tel que q[ia] <= a (-1 si aucun)
            size_t ua = gallop<true>(Q, nq, a, hint);
            size_t ub = gallop<true>(Q, nq, e, ua);
            hint = ua;
            int64_t ia = (int64_t)ua - 1, ib = (int64_t)ub - 1;
            const int64_t ns = (int64_t)nseg;
            if (ia == ib) {
                if (ia >= 0 && ia < ns) part[ia] += e - a;
                return;
            }
            if (ia >= 0 && ia < ns) part[ia] += Q[ia + 1] - a;
            if (ib >= 0 && ib < ns) part[ib] += e - Q[ib];
            int64_t lo = std::max<int64_t>(ia + 1, 0), hi = std::min<int64_t>(ib, ns);
            if (hi > lo) {
                ++full[lo];
                --full[hi];
            }
        };

        int r;
        while ((r = sam_read1(fp, hdr, b)) >= 0) {
            const bam1_core_t &c = b->core;
            if (o.max_reads && ++n_seen > o.max_reads) break;
            if (c.tid < 0) continue;
            if (c.tid >= nref) throw Fail(E_HEADER, "tid hors entete");
            if (tid_off[c.tid] < 0) {  // sequence absente du projet
                ++n_other_contig;
                continue;
            }
            const uint32_t *cig = bam_get_cigar(b);
            const uint32_t n_cig = c.n_cigar;
            const int64_t gpos = tid_off[c.tid] + c.pos;

            // ---------------- comptage (identique a bam.scan_chunk)
            if (do_count && !(c.flag & o.exclude) && c.qual >= o.min_mapq) {
                ++n_reads;
                if (n_cig == 1) {
                    int op = bam_cigar_op(cig[0]);
                    int64_t ln = bam_cigar_oplen(cig[0]);
                    if ((op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) && have_sites && ln >= anchor2)
                        add_site(gpos + o.anchor, gpos + ln - o.anchor, hint_sites);
                } else if (n_cig > 1) {
                    int64_t cur = gpos, blk = gpos, pend = 0;
                    bool spliced = false;
                    size_t h = hint_sites;
                    for (uint32_t i = 0; i < n_cig; ++i) {
                        int op = bam_cigar_op(cig[i]);
                        int64_t ln = bam_cigar_oplen(cig[i]);
                        if (op == BAM_CREF_SKIP) {
                            int64_t blen = cur - blk;
                            if (pend && blen >= o.min_anchor_split) {
                                junc.inc(pend);
                                spliced = true;
                            }
                            if (blen >= o.min_anchor_split && ln >= o.min_intron && ln <= o.max_intron)
                                pend = (cur << 24) | (ln & LEN_MASK);
                            else
                                pend = 0;
                            if (have_sites && blen >= anchor2) add_site(blk + o.anchor, cur - o.anchor, h);
                            cur += ln;
                            blk = cur;
                        } else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL) {
                            cur += ln;
                        }
                    }
                    if (pend && cur - blk >= o.min_anchor_split) {
                        junc.inc(pend);
                        spliced = true;
                    }
                    if (have_sites && cur - blk >= anchor2) add_site(blk + o.anchor, cur - o.anchor, h);
                    if (spliced) ++n_spliced;
                    // le curseur suit le debut des reads (monotone dans un BAM trie)
                    if (have_sites) hint_sites = gallop<false>(S, nsite, gpos, hint_sites);
                }
            }

            // ---------------- couverture (identique a coverage.scan_chunk_cov)
            if (do_cov && n_cig > 0 && !(c.flag & o.cov_exclude) && c.qual >= o.cov_mapq) {
                ++n_cov_reads;
                int64_t cur = gpos;
                size_t h = hint_seg;
                bool first = true;
                for (uint32_t i = 0; i < n_cig; ++i) {
                    int op = bam_cigar_op(cig[i]);
                    int64_t ln = bam_cigar_oplen(cig[i]);
                    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                        add_cov(cur, cur + ln, h);
                        if (first) {
                            hint_seg = h;
                            first = false;
                        }
                        cur += ln;
                    } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
                        cur += ln;
                    }
                }
            }
        }
        if (r < -1) throw Fail(E_READ, "erreur de lecture (disque deconnecte ou fichier corrompu ?) : " + o.bam);

        // ---------------- ecriture des resultats
        if (do_count) {
            std::vector<int64_t> k, cnt;
            junc.sorted(k, cnt);
            npy_save(o.count_prefix + ".jkey.npy", "<i8", {k.size()}, k.data(), k.size() * 8);
            npy_save(o.count_prefix + ".jcount.npy", "<i8", {cnt.size()}, cnt.data(), cnt.size() * 8);
            std::vector<int64_t> ns(nsite);
            int64_t acc = 0;
            for (size_t i = 0; i < nsite; ++i) {
                acc += sdiff[i];
                ns[i] = acc;
            }
            npy_save(o.count_prefix + ".nonsplit.npy", "<i8", {ns.size()}, ns.data(), ns.size() * 8);
        }
        if (do_cov) {
            // base[i] = bases partielles + (segments entierement couverts) x largeur
            std::vector<int64_t> base(nseg);
            int64_t run = 0;
            for (size_t i = 0; i < nseg; ++i) {
                run += full[i];
                base[i] = part[i] + run * (q[i + 1] - q[i]);
            }
            // sommes prefixes : bases / largeurs, parties non exoniques et exoniques
            std::vector<int64_t> cfb(nseg + 1, 0), cfw(nseg + 1, 0), ceb(nseg + 1, 0), cew(nseg + 1, 0);
            for (size_t i = 0; i < nseg; ++i) {
                const bool fr = seg_free[i] > 0;
                const int64_t w = q[i + 1] - q[i];
                cfb[i + 1] = cfb[i] + (fr ? base[i] : 0);
                cfw[i + 1] = cfw[i] + seg_free[i];
                ceb[i + 1] = ceb[i] + (fr ? 0 : base[i]);
                cew[i + 1] = cew[i] + (w - seg_free[i]);
            }
            const size_t M = idx_shape[0];
            std::vector<float> out(M * 3);
            for (size_t m = 0; m < M; ++m) {
                const int64_t *x = &idx[m * 8];
                // x = qi_s, qi_e, wi_l0, wi_l1, wi_r0, wi_r1, pt_l, pt_r
                int64_t s_sum = cfb[x[1]] - cfb[x[0]];
                int64_t s_len = cfw[x[1]] - cfw[x[0]];
                int64_t wb = (ceb[x[3]] - ceb[x[2]]) + (ceb[x[5]] - ceb[x[4]]);
                int64_t ww = (cew[x[3]] - cew[x[2]]) + (cew[x[5]] - cew[x[4]]);
                int64_t dl = (size_t)x[6] < nseg ? base[x[6]] : 0;
                int64_t dr = (size_t)x[7] < nseg ? base[x[7]] : 0;
                out[m * 3] = (float)((double)s_sum / (double)std::max<int64_t>(s_len, 1));
                out[m * 3 + 1] = (float)((double)(dl + dr) / 2.0);
                out[m * 3 + 2] = (float)((double)wb / (double)std::max<int64_t>(ww, 1));
            }
            npy_save(o.cov_out, "<f4", {M, 3}, out.data(), out.size() * 4);
        }

        double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        int64_t coff = 0;
        if (fp->format.compression == bgzf && fp->fp.bgzf) coff = bgzf_tell(fp->fp.bgzf) >> 16;
        std::printf("{\"n_reads\": %" PRId64 ", \"n_spliced\": %" PRId64 ", \"n_junctions\": %zu, "
                    "\"n_cov_reads\": %" PRId64 ", \"seconds\": %.3f, \"compressed_bytes\": %" PRId64
                    ", \"records\": %" PRId64 ", \"other_contig_records\": %" PRId64 "}\n",
                    n_reads, n_spliced, junc.size(), n_cov_reads, sec, coff, n_seen, n_other_contig);
    } catch (const Fail &f) {
        std::fprintf(stderr, "sfscan: %s\n", f.what());
        if (f.code == E_USAGE) usage();
        rc = f.code;
    } catch (const std::bad_alloc &) {
        std::fprintf(stderr, "sfscan: memoire insuffisante\n");
        rc = E_INPUT;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "sfscan: %s\n", e.what());
        rc = E_INPUT;
    }
    if (b) bam_destroy1(b);
    if (hdr) sam_hdr_destroy(hdr);
    if (fp && sam_close(fp) < 0 && rc == OK) rc = E_READ;
    return rc;
}
