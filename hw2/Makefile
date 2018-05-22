CXX=g++

TARGET = maxcut
INCLUDE_DIRS = ./

COMMON_FLAGS = $(foreach includedir, $(INCLUDE_DIRS), -I$(includedir))
CXXFLAGS = $(COMMON_FLAGS) -g -std=c++11 -O3 -Wall -fPIC 

########################################################
all: $(TARGET)

run:	
	./maxcut maxcut.in maxcut.out

$(TARGET): maxcut.o
	$(CXX)  -o  $@ $^ 

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $<


clean:
	rm -rf *.o $(TARGET)
