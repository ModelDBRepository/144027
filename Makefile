#
# Makefile for theta1d
#
DIR = .
LIBS = -lm

# directories to look in to find objects
DERIVS =    $(DIR)
NETWORK =   $(DIR)
PARAMS =    $(DIR)
OBJECTS =   $(DIR)
FUNCTIONS = $(DIR)

CC = g++
# CC = /opt/gcc-3.2.3/bin/g++
CFLAGS = -I$(DERIVS) -I$(OBJECTS) -I$(NETWORK) -I$(PARAMS) -I$(FUNCTIONS) -Wno-deprecated
C = $(CC) $(CFLAGS) -O -c

SRC = main.c
OBJ = main.o derivs.o network.o params.o rgauss.o runge_kutta.o lib.o

main:	$(OBJ)
	$(CC) -O $(CFLAGS) $(OBJ) $(LIBS) -o theta1d

main.o:	main.c $(OBJECTS)/float2.h
	$(C) main.c

derivs.o: $(DERIVS)/derivs.h $(DERIVS)/derivs.c
	$(C) $(DERIVS)/derivs.c

network.o: $(NETWORK)/network.h $(NETWORK)/network.c
	$(C) $(NETWORK)/network.c

params.o: $(PARAMS)/params.h $(PARAMS)/params.c
	$(C) $(PARAMS)/params.c

rgauss.o: $(OBJECTS)/rgauss.h $(OBJECTS)/rgauss.c
	$(C) $(OBJECTS)/rgauss.c

runge_kutta.o: $(OBJECTS)/runge_kutta.h $(OBJECTS)/runge_kutta.c
	$(C) $(OBJECTS)/runge_kutta.c

lib.o: $(FUNCTIONS)/lib.h
	$(C) $(FUNCTIONS)/lib.c
