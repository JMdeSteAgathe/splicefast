"""
Ecriture de classeurs .xlsx sans aucune dependance (zipfile + XML).

Pourquoi : les TSV s'ouvrent mal dans Excel des qu'un champ contient un
guillemet, une virgule decimale ou un retour a la ligne -- ce qui arrive tout
le temps avec les libelles PanelApp et OMIM, et decale les colonnes. Un vrai
classeur regle le probleme une fois pour toutes, avec en prime l'en-tete fige
et les filtres automatiques.

openpyxl ferait la meme chose, mais il n'est pas garanti present dans
l'environnement conda de l'utilisateur : 150 lignes ici evitent une dependance
de plus.
"""
from __future__ import annotations

import math
import zipfile
from pathlib import Path
from typing import Any, Iterable, List, Sequence, Tuple

MAX_ROWS = 1_048_575          # limite d'Excel (en-tete comprise)
MAX_CELL = 32_767             # longueur maxi d'une cellule


def _col(i: int) -> str:
    """0 -> A, 25 -> Z, 26 -> AA"""
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def _esc(x: str) -> str:
    out = (x.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
    # Excel refuse les caracteres de controle
    return "".join(c for c in out if c >= " " or c in "\t\n")


def _cell(ref: str, v: Any) -> str:
    if v is None or v == "":
        return ""
    if isinstance(v, bool):
        v = "yes" if v else "no"
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            return ""
        return f'<c r="{ref}"><v>{v!r}</v></c>' if isinstance(v, float) \
            else f'<c r="{ref}"><v>{v}</v></c>'
    s = str(v)
    if len(s) > MAX_CELL:
        s = s[:MAX_CELL - 3] + "..."
    return f'<c r="{ref}" t="inlineStr"><is><t xml:space="preserve">{_esc(s)}</t></is></c>'


def _sheet_xml(header: Sequence[str], rows: Iterable[Sequence[Any]]) -> str:
    n_col = len(header)
    parts = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
             '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">']
    widths = []
    for j, h in enumerate(header):
        widths.append(f'<col min="{j+1}" max="{j+1}" width="{min(60, max(9, len(str(h)) + 3))}" customWidth="1"/>')
    parts.append("<cols>" + "".join(widths) + "</cols>")
    parts.append("<sheetData>")
    parts.append('<row r="1">' + "".join(
        f'<c r="{_col(j)}1" t="inlineStr" s="1"><is><t>{_esc(str(h))}</t></is></c>'
        for j, h in enumerate(header)) + "</row>")
    r = 1
    for row in rows:
        r += 1
        if r > MAX_ROWS:
            break
        cells = "".join(_cell(f"{_col(j)}{r}", v) for j, v in enumerate(row[:n_col]))
        parts.append(f'<row r="{r}">{cells}</row>')
    parts.append("</sheetData>")
    if n_col:
        parts.append(f'<autoFilter ref="A1:{_col(n_col - 1)}{max(r, 1)}"/>')
    parts.append("</worksheet>")
    # volet fige sous l'en-tete
    parts.insert(1, '<sheetViews><sheetView workbookViewId="0" tabSelected="1">'
                    '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
                    '</sheetView></sheetViews>')
    return "".join(parts)


STYLES = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
          '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
          '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
          '<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
          '<fills count="2"><fill><patternFill patternType="none"/></fill>'
          '<fill><patternFill patternType="gray125"/></fill></fills>'
          '<borders count="1"><border/></borders>'
          '<cellStyleXfs count="1"><xf/></cellStyleXfs>'
          '<cellXfs count="2"><xf xfId="0"/>'
          '<xf xfId="0" fontId="1" applyFont="1"/></cellXfs>'
          '<cellStyles count="1">'
          '<cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
          '</styleSheet>')


def write_xlsx(path: Path, sheets: Sequence[Tuple[str, Sequence[str], Iterable[Sequence[Any]]]]) -> None:
    """sheets = [(nom, en-tete, lignes), ...]"""
    sheets = [s for s in sheets if s is not None]
    names, xmls = [], []
    for name, header, rows in sheets:
        nm = "".join(c for c in str(name) if c not in "[]:*?/\\")[:31] or "feuille"
        names.append(nm)
        xmls.append(_sheet_xml(header, rows))
    tmp = Path(str(path) + ".part")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        ct = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
              '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
              '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
              '<Default Extension="xml" ContentType="application/xml"/>'
              '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
              '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>']
        for i in range(len(names)):
            ct.append(f'<Override PartName="/xl/worksheets/sheet{i+1}.xml" '
                      f'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
        ct.append("</Types>")
        z.writestr("[Content_Types].xml", "".join(ct))
        z.writestr("_rels/.rels",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   '<Relationship Id="rId1" Target="xl/workbook.xml" '
                   'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument"/>'
                   '</Relationships>')
        wb = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
              '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
              'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>']
        rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
        for i, nm in enumerate(names):
            wb.append(f'<sheet name="{_esc(nm)}" sheetId="{i+1}" r:id="rId{i+1}"/>')
            rels.append(f'<Relationship Id="rId{i+1}" Target="worksheets/sheet{i+1}.xml" '
                        f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"/>')
        wb.append("</sheets></workbook>")
        rels.append(f'<Relationship Id="rId{len(names)+1}" Target="styles.xml" '
                    f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles"/>')
        rels.append("</Relationships>")
        z.writestr("xl/workbook.xml", "".join(wb))
        z.writestr("xl/_rels/workbook.xml.rels", "".join(rels))
        z.writestr("xl/styles.xml", STYLES)
        for i, x in enumerate(xmls):
            z.writestr(f"xl/worksheets/sheet{i+1}.xml", x)
    tmp.replace(path)
