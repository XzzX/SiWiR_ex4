SHELL        = /bin/bash
CXX          = mpic++
 
FLAGS        = -march=native -mtune=native -std=c++0x -Wall -Wextra -Wshadow -Werror -O3 -DNDEBUG
DEBUGFLAGS   = 
RELEASEFLAGS =
INCPATH      = 
LIBPATH      = 
LIBS         = 

# likwid
FLAGS       += -DUSE_LIKWID -pthread
INCPATH     += -I/usr/local/likwid-3.1.2/include/
LIBPATH     += -L/usr/local/likwid-3.1.2/lib/
LIBS        += -llikwid
 
COMMON       = Timer.h
 
all: cg

%.o: %.cpp $(COMMON)
	$(CXX) -c $(FLAGS) $(INCPATH) $<
 
cg: cg.cpp
	$(CXX) $(FLAGS) $(INCPATH) -o cg cg.cpp $(LIBPATH) $(LIBS)

run: cg
	./cg ${ARGS}
	gnuplot gnuplotSolution.plt

documentation: report/setup.tex report/main.tex
	cd report && pdflatex main.tex
	cd report && pdflatex main.tex
	cd report && pdflatex main.tex

clean:
	rm -f *.o cg

.PHONY : all clean
