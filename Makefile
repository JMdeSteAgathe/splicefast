# SpliceFast : compilation du lecteur C++ (sfscan)
#
#   make                          htslib trouvee par pkg-config (apt/conda)
#   make HTSLIB=/chemin/htslib    htslib compilee a la main
#   ./install_htslib.sh           tout compiler sans droits root (third_party/)
#
CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra -Wno-unused-parameter
BIN      := splicefast/bin/sfscan

# htslib compilee par ./install_htslib.sh : utilisee automatiquement
ifeq ($(HTSLIB),)
  ifneq ($(wildcard third_party/htslib/libhts.a),)
    HTSLIB := third_party/htslib
  endif
endif

ifdef HTSLIB
  HTS_CFLAGS := -I$(HTSLIB)/include -I$(HTSLIB)
  ifneq ($(wildcard $(HTSLIB)/libhts.a),)
    HTS_LIBS := $(HTSLIB)/libhts.a
  else
    HTS_LIBS := $(HTSLIB)/lib/libhts.a
  endif
  # dependances du libhts.a (bz2/lzma/curl/deflate selon la compilation)
  HTS_LIBS += $(shell sed -n 's/^HTSLIB_static_LDFLAGS *= *//p' $(HTSLIB)/htslib_static.mk 2>/dev/null)
  HTS_LIBS += $(shell sed -n 's/^HTSLIB_static_LIBS *= *//p' $(HTSLIB)/htslib_static.mk 2>/dev/null) -lz -lpthread -lm
else
  HTS_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
  HTS_LIBS   := $(shell pkg-config --libs htslib 2>/dev/null || echo -lhts) -lpthread
endif

.PHONY: all clean distclean test

all: $(BIN)

$(BIN): src/sfscan.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(HTS_CFLAGS) $< -o $@ $(HTS_LIBS)
	@echo "OK : $@"

test: $(BIN)
	python3 -m splicefast selftest -w /tmp/splicefast_selftest

clean:
	rm -f $(BIN)

distclean: clean
	rm -rf third_party
