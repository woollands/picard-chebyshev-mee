#
#  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
#  DATE WRITTEN:     Oct 2019
#  LAST MODIFIED:    Oct 2019
#  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
#  DESCRIPTION:      File to compile Picard-Chebyshev MEE code
#

# CC=clang
# CFLAGS=-Xpreprocessor -fopenmp -lomp
#
# obj= test_mee.o
#
# test: $(obj)
# 	 $(CC) ${CFLAGS} -o test $(obj)
#
# %.o:%.c
# 	 $(CC) ${CFLAGS} -c $<
#
# .PHONY: clean
#
# clean:
# 		rm *.o

IDIR = ../include
# CC = clang
CC = gcc -std=c99 -g
# CC = gcc-9 -std=c99 -g -fopenmp -lomp
# CC = /usr/local/opt/llvm/bin/clang -fopenmp -L/usr/local/opt/llvm/lib gcc -std=c99 -g
# CC = /usr/local/opt/llvm/bin/clang -fopenmp gcc -std=c99 -g
# CC = clang -Xpreprocessor -fopenmp -lomp
# clang -Xpreprocessor -fopenmp -lomp test_mee.c -o test

CFLAGS = -I$(IDIR)
#CFLAGS = -I$(IDIR) -Xpreprocessor -fopenmp -lomp

ODIR = obj
LDIR = ../lib

LIBS = -lm

_DEPS1 = const.h c_functions.h classical2mee.h matrix_loader.h mee2rv.h ecef2eci.h eci2ecef.h EGM2008.h \
inertial2radial.h perturbed_gravity.h radial_gravity.h
DEPS1 = $(patsubst %,$(IDIR)/%,$(_DEPS1))

_DEPS2 = clenshaw_curtis_ivpI.h lsq_chebyshev_fit.h chebyshev.h c_functions.h
DEPS2 = $(patsubst %,$(IDIR)/%,$(_DEPS2))

_OBJ1 = test_mee.o c_functions.o classical2mee.o matrix_loader.o mee2rv.o ecef2eci.o eci2ecef.o EGM2008.o \
inertial2radial.o perturbed_gravity.o radial_gravity.o
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
