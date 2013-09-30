#ViennaRNA package location
ViennaRNA=~/software/include/ViennaRNA
ViennaLIB=~/software/lib/libRNA.a

NAME=BHGbuilder

OBJ = $(NAME)_cmdline.o\
			main.o\
			move_set.o\
			$(NAME).o\
			RNAutils.o\
			flooding.o\
			graph_tech.o\
			RNAstruc.o\
			clustering.o

# should stay unchanged: 
CPP = g++
CC = gcc
CPPFLAGS =  -O2 -g -Wall -std=c++0x -fexceptions -Wno-write-strings
CFLAGS =  -O2 -g -Wall -fexceptions -Wno-write-strings
LFLAGS = -O2 -fopenmp

DIRS = -I $(ViennaRNA)

LIBS = $(ViennaLIB)

all: $(OBJ)
	$(CPP) $(LFLAGS) $(DIRS) $(OBJ) $(LIBS) -o $(NAME)
	rm -f $(OBJ)

$(NAME)_cmdline.h $(NAME)_cmdline.c: $(NAME).ggo
	gengetopt -i $(NAME).ggo

%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(DIRS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) $(DIRS) -c $<

clean:
	rm -f $(OBJ)
	rm -f $(NAME)
	rm -f $(NAME)_cmdline.c $(NAME)_cmdline.h

