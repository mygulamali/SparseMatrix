# Makefile for SparseMatrix
# Murtaza Gulamali (11/11/2011)

CC   = g++
LIBS = -lm -lpthread
DEFS =

CPPFLAGS += -I${GTEST_DIR}/include
CXXFLAGS += -O3 -g -Wall -Wextra ${LIBS} ${DEFS}

GTEST_DIR = ./gtest-1.6.0

GTEST_HEADERS = ${GTEST_DIR}/include/gtest/*.h \
                ${GTEST_DIR}/include/gtest/internal/*.h

GTEST_SRCS_ = ${GTEST_DIR}/src/*.cc ${GTEST_DIR}/src/*.h ${GTEST_HEADERS}

SOURCES = SparseMatrix.cpp

.SUFFIXES:
.SUFFIXES: .cpp .o

OBJECTS = ${SOURCES:.cpp=.o}

all: test demo

.cpp.o:
	${CC} -c ${CPPFLAGS} ${CXXFLAGS} $<

test: gtest_main.a ${OBJECTS}
	${CC} ${CPPFLAGS} ${CXXFLAGS} -o $@ ${OBJECTS} gtest_main.a test.cpp

demo: ${OBJECTS}
	${CC} ${CXXFLAGS} -o $@ ${OBJECTS} demo.cpp

gtest-all.o: ${GTEST_SRCS_}
	${CC} ${CPPFLAGS} -I${GTEST_DIR} ${CXXFLAGS} -c ${GTEST_DIR}/src/gtest-all.cc

gtest_main.o: ${GTEST_SRCS_}
	${CC} ${CPPFLAGS} -I${GTEST_DIR} ${CXXFLAGS} -c ${GTEST_DIR}/src/gtest_main.cc

gtest.a: gtest-all.o
	${AR} ${ARFLAGS} $@ $^

gtest_main.a: gtest-all.o gtest_main.o
	${AR} ${ARFLAGS} $@ $^

clean:
	rm *.o *.a

distclean: clean
	rm demo test
