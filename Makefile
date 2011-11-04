# Simple Makefile for SparseMatrix
# Murtaza Gulamali (02/12/2006)

CC      = g++
INCLUDE =
LIBS    = -lm
DEFS    =
FLAGS   = -w -O3 ${INCLUDE} ${LIBS} ${DEFS}

SOURCES = SparseMatrix.cpp

.SUFFIXES:
.SUFFIXES: .cpp .o

OBJECTS = ${SOURCES:.cpp=.o}

all: SparseMatrixDemo

.cpp.o:
	${CC} -c ${FLAGS} $<

SparseMatrixDemo: ${OBJECTS}
	${CC} ${FLAGS} -o $@ ${OBJECTS} SparseMatrixDemo.cpp

clean:
	rm *.o

distclean: clean
	rm ${BIN}
