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
 
all: heat

%.o: %.cpp $(COMMON)
	$(CXX) -c $(FLAGS) $(INCPATH) $<
 
heat: heat.cpp
	$(CXX) $(FLAGS) $(INCPATH) -o heat heat.cpp $(LIBPATH) $(LIBS)

run: heat
	./heat ${ARGS}
	gnuplot gnuplotSolution.plt

documentation: report/setup.tex report/main.tex
	cd report && pdflatex main.tex
	cd report && pdflatex main.tex
	cd report && pdflatex main.tex

clean:
	rm -f *.o heat

.PHONY : all clean
