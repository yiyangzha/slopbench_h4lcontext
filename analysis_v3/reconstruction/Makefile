# Build the v3 4-lepton reconstruction programs against the lxplus system ROOT.
ROOTCONFIG := /usr/bin/root-config
# The system compiler (the ROOT build's), independent of PATH.
CXX := /usr/bin/g++
CXXFLAGS := -O2 -g -Wall -Wextra $(shell $(ROOTCONFIG) --cflags) -I../common/include
LDLIBS := $(shell $(ROOTCONFIG) --libs) -Wl,-rpath,$(shell $(ROOTCONFIG) --libdir)
HEADERS := $(wildcard ../common/include/h4l/*.h)
# JHUGenMELA wrapper (production_v3/external/mela_wrapper, see its README.md);
# the /eos/user form is compiled in as the rpath.
MELA_WRAPPER := /eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l/production_v3/external/mela_wrapper
MELAFLAGS := -I$(MELA_WRAPPER) -L$(MELA_WRAPPER)/lib -Wl,-rpath,$(MELA_WRAPPER)/lib -lMelaWrapper
PROGRAMS := bin/h4l_reco bin/h4l_select

all: $(PROGRAMS)

bin/h4l_select: src/h4l_select.cpp $(HEADERS) $(MELA_WRAPPER)/MelaWrapper.h $(MELA_WRAPPER)/lib/libMelaWrapper.so
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $< $(MELAFLAGS) $(LDLIBS) -lMinuit2

bin/%: src/%.cpp $(HEADERS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

.PHONY: all
