EXECUTABLE=fmg
SOURCES=multigrid.c residual.c smoothing.c matfill.c prolongation.c restriction.c solve.c init.c output.c vcycle.c addint.c main.c

CC=gcc
INCLIB=
LDLIB=

OPENMPLIB=-lgomp
MATHLIB=-lm
WARNFLAGS=-Wall -Wextra
OPENMPFLAG=-fopenmp
OPTFLAG=-O3


LDFLAGS= $(MATHLIB) 
CFLAGS=-c $(OPTFLAG) $(WARNFLAGS) $(HDF5FLAG)  




BIN=bin/
SRC=src/




#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CHEADER=$(addprefix $(SRC),$(HEADER))




all: $(CSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS)
	$(CC)  $(COBJECTS) $(LDLIB) $(LDFLAGS) -o $@ 

$(BIN)%.o: $(SRC)%.c $(CHEADER)
	$(CC) $(INCLIB) $(CFLAGS) $< -o $@ 


clean:
	rm $(COBJECTS) $(EXECUTABLE)
