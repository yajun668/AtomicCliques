# For cowboy
GUROBI_HOME =/opt/gurobi/9.5.1/linux64
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

GPP = /usr/bin/g++

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

OPT = -m64 -std=c++11 -Wextra

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
CFLAGS=-Wall -O2 -fopenmp
INSTANCE  = ./
INC     = $(GUROBI_HOME)/include
OPT_INC = $(OPT) -I$(INC)
LIBDIR  = $(GUROBI_HOME)/lib/
LIB     = -L$(LIBDIR) -lgurobi_c++ -lgurobi95 -lm -lpthread
SOURCES= main.cpp model.cpp common.cpp graph.cpp
HEADERS= model.h common.h graph.h
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(GPP) $(CFLAGS) $(OPT_INC) $(OBJECTS) -o $@ $(LIB)

%.o:  %.cpp  $(HEADERS)
	$(GPP) $(CFLAGS) $(OPT_INC) -c $< -o $@ 


	
clean :
	/bin/rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp
	/bin/rm -rf *~ *.mps *.ord *.sos *.sav *.net *.msg *.clp
