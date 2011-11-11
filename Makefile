# Makefile for SparseMatrix
# Murtaza Gulamali (11/11/2011)

CC      = g++
INCLUDE = -Iinclude
LIBS    = -lm
DEFS    =
FLAGS   = -w -O3 ${INCLUDE} ${LIBS} ${DEFS}

SOURCES = SparseMatrix.cpp

.SUFFIXES:
.SUFFIXES: .cpp .o

OBJECTS = ${SOURCES:.cpp=.o}

all: demo

.cpp.o:
	${CC} -c ${FLAGS} $<

demo: ${OBJECTS}
	${CC} ${FLAGS} -o $@ ${OBJECTS} demo.cpp

clean:
	rm *.o

distclean: clean
	rm demo
