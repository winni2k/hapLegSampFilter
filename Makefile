### Makefile --- 

## Author: wkretzsch@gmail.com
## Version: $Id: Makefile,v 0.0 2014/09/12 14:05:24 winni Exp $
## Keywords: 
## X-URL: 


### Makefile ends here

BIN=hapLegSampFilter
all: $(BIN)
clean: 
	rm -f $(BIN) *.o

CXX= g++
CXXFLAGS = -Wall -std=c++11 -O3
LIBS = -lz -lboost_iostreams -lboost_program_options

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

$(BIN): $(BIN).o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)
