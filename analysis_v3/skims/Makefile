# Build the v3 skim programs against the lxplus system ROOT.
ROOTCONFIG := /usr/bin/root-config
CXXFLAGS := -O2 -g -Wall -Wextra $(shell $(ROOTCONFIG) --cflags) -I../common/include
LDLIBS := $(shell $(ROOTCONFIG) --libs) -Wl,-rpath,$(shell $(ROOTCONFIG) --libdir)
HEADERS := $(wildcard ../common/include/h4l/*.h)
PROGRAMS := bin/skim_v3 bin/skim_histograms

all: $(PROGRAMS)

bin/%: src/%.cpp $(HEADERS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

.PHONY: all
