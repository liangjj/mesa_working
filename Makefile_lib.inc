#
#@(#)Makefile	5.1   11/6/94
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
NOLINK = -c
ALL_LIBS = \
     $(SOURCE_LIB) \
     $(MOD_LIB) \
     $(POTENTIAL_LIB) \
     $(UTILITY_LIB) \
     $(BLIU_LIB) \
     $(MD_LIB)

ALL_SRCS= $(SRCS.F90:.F90=.o) $(SRCS.f90:.f90=.o) $(SRCS:.f=.o) $(SRCS.f:.f=.o) $(SRCS.F:.F=.o) 
ALL_LIB_SRCS= $(LIB_SRCS.F90:.F90=.o) $(LIB_SRCS.f90:.f90=.o) $(LIB_SRCS:.f=.o) $(LIB_SRCS.f:.f=.o) $(LIB_SRCS.F:.F=.o) 

COMPILE=$(F90) $(ALL_INCLUDE) $(NOLINK) $(F90FLAGS) $(CPPFLAGS)

LOAD=$(F90) $(ALL_INCLUDE) $(F90FLAGS) $(ALL_SRCS) $(ALL_LIBS) -o

.f.o:
	$(COMPILE) $<
.F.o:
	$(COMPILE) $<
.f90.o:
	$(COMPILE) $<
.F90.o:
	$(COMPILE) $<

GET = sccs get
REL = 

all: $(EXECUTABLE)
$(EXECUTABLE): $(ALL_LIB_SRCS) $(ALL_SRCS)
	$(LOAD) $(EXECUTABLE)

lib: $(SOURCE_LIB)
$(SOURCE_LIB): $(ALL_LIB_SRCS)
	$(AR) $(ARFLAGS) $(SOURCE_LIB) $(?:.f=.o) $(?:.f90=.o)
	$(RANLIB) $(SOURCE_LIB)

main: $(EXECUTABLE)
$(EXECUTABLE): $(ALL_LIB_SRCS) $(ALL_SRCS)
	$(LOAD) $(EXECUTABLE)

clean: 
	rm -fR *.a *.o *.mod