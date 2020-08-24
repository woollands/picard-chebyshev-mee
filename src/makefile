#
#  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
#  DATE WRITTEN:     Oct 2019
#  LAST MODIFIED:    Oct 2019
#  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
#  DESCRIPTION:      File to compile Picard-Chebyshev MEE code
#

CC = gcc -std=c99 -g

IDIR = ../include
ODIR = obj
LDIR = ../lib
LIBS = -lm

CFLAGS = -I$(IDIR)

_DEPS1 = const.h c_functions.h classical2mee.h matrix_loader.h mee2rv.h ecef2eci.h eci2ecef.h EGM2008.h \
inertial2radial.h
DEPS1 = $(patsubst %,$(IDIR)/%,$(_DEPS1))

_DEPS2 = clenshaw_curtis_ivpI.h lsq_chebyshev_fit.h chebyshev.h c_functions.h
DEPS2 = $(patsubst %,$(IDIR)/%,$(_DEPS2))

_OBJ1 = test_mee.o c_functions.o classical2mee.o matrix_loader.o mee2rv.o ecef2eci.o eci2ecef.o GM.o \
inertial2radial.o
OBJ1 = $(patsubst %,$(ODIR)/%,$(_OBJ1))

_OBJ2 = matrix_builder.o clenshaw_curtis_ivpI.o lsq_chebyshev_fit.o chebyshev.o c_functions.o clenshaw_curtis_ivpI.o
OBJ2 = $(patsubst %,$(ODIR)/%,$(_OBJ2))

$(ODIR)/%.o: %.c $(DEPS1) $(DEPS2)
	$(CC) -c -o $@ $< $(CFLAGS)

test_mee: $(OBJ1)
	gcc -g -o $@ $^ $(CFLAGS) $(LIBS)

matrix_builder: $(OBJ2)
	gcc -g -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
