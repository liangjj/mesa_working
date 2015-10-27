#   @(#)Makefile	5.1 11/6/94
#
#    Makefile MesaLib
#    
RANLIB = ranlib
OBJECTS = *.o

all: $(MESA_LIB) 
$(MESA_LIB): $(OBJECTS) 
	$(AR) $(ARFLAGS) $(MESA_LIB) $(OBJECTS)
	$(RANLIB) $(MESA_LIB)
clean:
	rm -f *.o
