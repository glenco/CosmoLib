INCL = -I./include/ \
	-I../NR/include

OUT = libCosmoLib.a

SRC = Cosmo/powerEH.cpp \
Cosmo/powerEHv2.cpp \
cosmo.cpp \
halo.cpp \
nfw.cpp \
powerCDM.cpp \
utilities.cpp

#
OBJ = $(SRC:.cpp=.o)
# 
OPT = -O2
# 
DEBUG = -g 
# compiler  
CC = g++
#
RM = rm -fr
#

CFLAGS = $(INCL) $(DEBUG) $(OPT)
#
CLEAR = clear

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

$(OUT): $(OBJ)
	ar rcs $(OUT) $(OBJ)

clean:
	$(RM) $(OBJ) $(OUT) *~

help:
	$(CLEAR)	
	@echo 
