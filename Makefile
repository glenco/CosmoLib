NR_DIR = ../NR

INCL = -I./include/ \
	-I$(NR_DIR)/include \
	-I/usr/local/include

OUT = libCosmoLib.a

SRC = Cosmo/powerEH.cpp \
Cosmo/powerEHv2.cpp \
cosmo.cpp \
halo.cpp \
nfw.cpp \
powerCDM.cpp \
utilities.cpp \
powerCDMHM.cpp

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
	$(CC) $(OPT) $(CFLAGS) -c $< -o $@

all: $(OBJ)
	ar rcs $(OUT) $(OBJ)

clean:
	$(RM) $(OBJ) $(OUT) *~

