# Makefile 

# Variables
PROGRAM = proggis
CC      = mpif90 
CFLAGS  = -O3 
LIBS    = -llapack -lblas 

OBJS    =globaali.f90   \
	GPfemdvr.f90   \
	quadrature.f90 \
	lagrange.f90   \
	kanta.f90      \
	ptswts.f90     \
	init_wf.f90    \
	ke_lohkot.f90  \
	ke_props.f90   \
	mom_lohkot.f90 \
	mom_props.f90  \
	propagate.f90  \
	prop_so2.f90   \
	prop_so4.f90   \
	p_kinetic_p.f90 \
	p_kinetic_q.f90 \
	pot_propagation.f90 \
	kin_propagation.f90 \
	mom_propagation.f90 \
	ang_propagation.f90 \
	normalize.f90  \
	virhenormi.f90 \
	seiv_slice.f90 \
	seiv_wfn.f90

# Compile
$(PROGRAM): ${OBJS}
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)
# Clean
clean:
	rm -f $(PROGRAM) ${OBJS}
