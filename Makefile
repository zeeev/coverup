######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################


CC=g++
FLAGS=-std=c++0x -g 
LD=-lsplit -lgencode -lm

HEADERS := $(wildcard lib/*hpp)
SOURCE  := $(wildcard lib/*cpp)
OBJECTS := $(patsubst %.cpp,%.o,$(wildcard lib/*cpp))
BIN     := $(wildcard src/*cpp)

.PHONY: all

all: mkbin bin clean gitclean

bin: libsplit.a libgencode.a
	$(CC) $(BIN) *.a -Ilib/ -I. -o bin/coverUp
mkbin:
	-mkdir bin

$(OBJECTS):
	$(CC) -c $(FLAGS) $(SOURCE)
libsplit.a: $(OBJECTS)
	ar rcs libsplit.a split.o
libgencode.a: $(OBJECTS)
	ar rcs libgencode.a genCodeClass.o
clean:
	-rm *.a && rm *.o && rm data/gencode.v19.annotation.gtf.gindx
gitclean:
	-rm *~ && rm src/*~ && rm lib/*~ 