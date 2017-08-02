IDIR =-Ilib/Eigen -Ilib -IC:/Ruby20/include/ruby-2.0.0 -ID:/Programs/OpenCV2.4/include -I.
LIBS=-lstdc++
CC=gcc -std=c++11
CFLAGS=$(IDIR) -O3

ODIR =obj


_DEPS = cobimg.h photomeasure.h measurement.h
DEPS = $(patsubst %,core/%,$(_DEPS))

_OBJ = cobimg.o photomeasure.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

obj/cobimg.o: core/cobimg.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)
    
obj/photomeasure.o: core/photomeasure.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 