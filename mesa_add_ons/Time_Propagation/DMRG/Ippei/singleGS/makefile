OFLAGS = --o2
LIBFLAGS = -llapack -lblas 
#MC = mpif90
MC = mpif77 #--tp4
MC = $(F90)

MPDFortran : system_parameters.o input_output.o cebct.o ebct.o Matrix_Exponential_Module.o MPDtools_module.o local_operations_module.o io_module.o observables_module.o propagation_module.o MPDFortran.o
	$(MC) $(OFLAGS) system_parameters.o input_output.o cebct.o ebct.o Matrix_Exponential_Module.o MPDtools_module.o local_operations_module.o io_module.o observables_module.o propagation_module.o MPDFortran.o -o MPDFortran $(LIBFLAGS) 

MPDFortran.o : MPDFortran.f90 
	$(MC) $(OFLAGS) -c MPDFortran.f90 

system_parameters.o : system_parameters.f90
	$(MC) $(OFLAGS) -c system_parameters.f90

MPDtools_module.o : MPDtools_module.f90
	$(MC) $(OFLAGS) -c MPDtools_module.f90

observables_module.o : observables_module.f90
	$(MC) $(OFLAGS) -c observables_module.f90

local_operations_module.o : local_operations_module.f90
	$(MC) $(OFLAGS) -c local_operations_module.f90

propagation_module.o : propagation_module.f90
	$(MC) $(OFLAGS) -c propagation_module.f90

io_module.o : io_module.f90
	$(MC) $(OFLAGS) -c io_module.f90

Matrix_Exponential_Module.o : Matrix_Exponential_Module.f90
	$(MC) $(OFLAGS) -c Matrix_Exponential_Module.f90

ebct.o : ebct.f90
	$(MC) $(OFLAGS) -c ebct.f90

cebct.o : cebct.f90
	$(MC) $(OFLAGS) -c cebct.f90
	
#prntcmn.o : prntcmn.f90
#	$(MC) $(OFLAGS) -c prntcmn.f90

#prntfmn.o : prntfmn.f90
#	$(MC) $(OFLAGS) -c prntfmn.f90
	
input_output.o : input_output.f90
	$(MC) $(OFLAGS) -c input_output.f90
