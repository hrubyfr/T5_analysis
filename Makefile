CXX = g++
CXXFLAGS = -O3 -Wall -pipe -g $(shell root-config --cflags)

LIBS = $(shell root-config --glibs)

TARGET = analyze_T5
OBJS = analyze_T5.o utils.o return_TOF_position.o buffer.o FitParameters.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)

analyze_T5.o: analyze_T5.cpp utils.h
	$(CXX) $(CXXFLAGS) -c $<

utils.o: utils.cpp utils.h
	$(CXX) $(CXXFLAGS) -c $<
return_TOF_position.o: return_TOF_position.cpp return_TOF_position.h
	$(CXX) $(CXXFLAGS) -c $<
buffer.o: buffer.cpp buffer.h
	$(CXX) $(CXXFLAGS) -c $<
FitParameters.o: FitParameters.cpp FitParameters.h
	$(CXX) $(CXXFLAGS) -c $<

clean: 
	rm -f $(OBJS) $(TARGET)
