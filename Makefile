NR_DIR = /Users/cgiocoli/Dropbox/Astro/Codes/EclipseWorkspace/NR
GSL_DIR = /Users/cgiocoli/Library/CppLib/gsl-1.15/

INCL = -I./include/ \
       -I$(NR_DIR)/include \
       -I$(GSL_DIR)/include \
       -I$(GSL_DIR)/include/ 

OUT = libCosmoLib.a

SRC = Cosmo/powerEH.cpp \
Cosmo/powerEHv2.cpp \
cosmo.cpp \
halo.cpp \
nfw.cpp \
powerCDM.cpp \
utilities.cpp \
powerCDMHM.cpp   # <--- if you include GSL you can use also this Carlo

#
OBJ = $(SRC:.cpp=.o)
# 
OPT = -O2 -DGSL
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

all: $(OBJ)
	ar rcs $(OUT) $(OBJ)

clean:
	$(RM) $(OBJ) $(OUT) *~

