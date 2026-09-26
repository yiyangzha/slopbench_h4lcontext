# Build the v3 calibration programs against the lxplus system ROOT.
ROOTCONFIG := /usr/bin/root-config
CXXFLAGS := -O2 -g -Wall -Wextra $(shell $(ROOTCONFIG) --cflags) -I../common/include
LDLIBS := $(shell $(ROOTCONFIG) --libs) -Wl,-rpath,$(shell $(ROOTCONFIG) --libdir)
FITLIBS := -lMinuit2
HEADERS := $(wildcard ../common/include/h4l/*.h)
PROGRAMS := bin/calib_extract bin/zpeak_histograms bin/fit_zpeaks bin/template_fit bin/lambda_histograms bin/fit_lambda

all: $(PROGRAMS)

bin/fit_zpeaks: src/fit_zpeaks.cpp $(HEADERS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS) $(FITLIBS)

bin/template_fit: src/template_fit.cpp $(HEADERS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS) $(FITLIBS)

bin/fit_lambda: src/fit_lambda.cpp $(HEADERS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS) $(FITLIBS)

bin/%: src/%.cpp $(HEADERS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

.PHONY: all
