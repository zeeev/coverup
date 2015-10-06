######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=g++
FLAGS=-std=c++0x -g 
LD= -lm -lz  -lpthread


HEADERS := $(wildcard lib/*hpp)
SOURCE  := $(wildcard lib/*cpp)
OBJECTS := $(patsubst %.cpp,%.o,$(wildcard lib/*cpp))
BIN     := $(wildcard src/*cpp)

.PHONY: all

all:  mkbin bin 

gitinit:
	-cd tabixpp && git submodule init && git submodule update

bin: libsplit.a libgencode.a libtabix.a libhts.a libyagr.a
	$(CC) $(FLAGS) $(LD) $(BIN) *.a  -Ilib/  -I. -Itabixpp/htslib -Ltabixpp/htslib -lhts -o bin/coverUp
mkbin:
	-mkdir bin

maketabix: gitinit
	-cd tabixpp && make
$(OBJECTS):
	$(CC) -c $(FLAGS) $(SOURCE)

libhts.a: maketabix
	cp tabixpp/htslib/libhts.a .
libtabix.a: maketabix
	ar rcs libtabix.a tabixpp/tabix.o
libsplit.a: $(OBJECTS)
	ar rcs libsplit.a split.o
libgencode.a: $(OBJECTS)
	ar rcs libgencode.a genCodeClass.o
libyagr.a:
	ar rsc libyagar.a yagbv.o
clean:
	-rm *.a && rm *.o && rm data/gencode.v19.annotation.gtf.gindx
gitclean:
	-rm *~ && rm src/*~ && rm lib/*~