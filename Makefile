CPP=g++ -std=c++14
CPPFLAGS=-O3 -Wall -DVERBOSE
INCLUDES=-I/home/daniela/include/
LIB=/home/daniela/lib/libsdsl.a /home/daniela/lib/libdivsufsort.a /home/daniela/lib/libdivsufsort64.a
OBJECTS=pfordelta.o
BINS=build_lib

%.o: %.cpp
	@echo " [C++] Compiling $<"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

all: clean stats $(OBJECTS) $(BINS)

stats:
	@echo
	@echo " COMPILING pfordelta.a"
	@echo " ###################"
	@echo "  * Compiler flags: $(CPPFLAGS)"
	@echo "  * Include dirs: $(INCLUDES)"
	@echo "  * Lib dirs: $(LIB)"
	@echo

clean:
	@echo " [CLN] Removing object files"
	@rm -f $(OBJECTS) $(BINS)

build_lib: 
	@echo " [BLD] Building dlhipil.a"
	ar -rvcs pfordelta.a $(OBJECTS) $(LIB) 

test:
	@echo " [BLD] Building binary testPforDelta"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o testPforDelta test.cpp $(OBJECTS) $(LIB)
