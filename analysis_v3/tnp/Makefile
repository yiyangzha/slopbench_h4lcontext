# Build the v3 tag-and-probe programs against the lxplus system ROOT.
ROOTCONFIG := /usr/bin/root-config
CXXFLAGS := -O2 -g -Wall -Wextra $(shell $(ROOTCONFIG) --cflags) -I../common/include
LDLIBS := $(shell $(ROOTCONFIG) --libs) -Wl,-rpath,$(shell $(ROOTCONFIG) --libdir)
HEADERS := $(wildcard ../common/include/h4l/*.h)
PROGRAMS := bin/tnp_extract bin/tnp_histograms bin/fit_tnp

all: $(PROGRAMS)

bin/fit_tnp: src/fit_tnp.cpp $(HEADERS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS) -lMinuit2

bin/%: src/%.cpp $(HEADERS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

.PHONY: all
