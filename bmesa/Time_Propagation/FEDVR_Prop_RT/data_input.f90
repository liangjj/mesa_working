  PROGRAM data_input_propagator
  IMPLICIT NONE
  CHARACTER (LEN=8)                        :: itoc, pr_main
  CHARACTER (LEN=80)                       :: chrkey, units, typke, diag_mod
  LOGICAL                                  :: logkey
  LOGICAL                                  :: genpts, nlse, plot, imtime
  LOGICAL                                  :: space, diag
  LOGICAL                                  :: absorb, proj
  INTEGER                                  :: spdim, prop_order
  INTEGER                                  :: plot_step, intkey, i, iostat
  CHARACTER(LEN=80)                        :: input_file_name, line
!
  WRITE(5,*) '          Terminal Input of Variables'
  WRITE(5,*) '          Input Filename'
  READ(5,*) input_file_name
  OPEN (UNIT=6,FILE=input_file_name,ACCESS='sequential',FORM='formatted',       &
        IOSTAT=IOSTAT,STATUS='new')
  WRITE(6,*), '$route'
  WRITE(6,*)
  WRITE(6,*) '$end'
  WRITE(6,*) '$nonstd'
  WRITE(5,*) 'Enter Route'
  line=' '
  DO WHILE(line /= '$end')
     READ(5,*) line
     WRITE(6,*) line
  END DO
  WRITE(6,*) '$end'
  WRITE(5,*) 'Enter Title Card'
  READ(5,*) line
  WRITE(6,*) line
  WRITE(6,*) "$end"
  WRITE(6,*) '$dvrprop_basis'
  WRITE(5,*) '$dvrprop_basis data'
  WRITE(5,*) 'number-of-space-variables, propagation-order'
  READ(5,*) spdim, prop_order
  WRITE(6,*) 'number-of-space-variables = ',spdim,' propagation-order = ',prop_order
  WRITE(5,*) 'coordinate system'
  READ(5,*) line
  WRITE(6,*) 'coordinate system = ',line
  DO i=1, spdim
     WRITE(5,*) 'label for space-variable-'//itoc(i)
     READ(5,*) line
     WRITE(6,*) 'space-variable-'//itoc(i)//' = ',line
  END DO
  WRITE(5,*) 'units, kinetic-energy-type, automate-points,     &
              non-linear-equations'     
  WRITE(5,*) 'plot, imaginary-time'     
  READ(5,*) units, typke, genpts, nlse, plot, imtime
  WRITE(6,*) 'units = ',units,' kinetic-energy-type = ',typke
  IF(genpts) THEN
     line='automate-points'
  END IF
  IF(nlse) THEN
     line = line//' non-linear-equations'
  END IF
  IF(plot) THEN
     line = line//' plot'
     WRITE(5,*) 'plot step'
     READ(5,*)  plot_step
     line=line//'plot-step='//itoc(plot_step)
  END IF
  IF(imtime) THEN
     line=line//' imaginary-time'
  END IF
  WRITE(6,*) line
  WRITE(5,*) 'no-spatial-hamiltonian, modify-diagonals, get-eigenpairs'
  READ(5,*) space, diag_mod, diag
  IF(space) THEN
     line='no-spatial-hamiltonian'
  END IF
  IF(diag_mod /='none') THEN
     line=line//' diagonal-modification=//diag_mod
  END IF
  IF(diag) THEN
     line=line//' get-eigenpairs'
  END IF
  WRITE(6,*) line
  WRITE(5,*) 'print parameter'
  READ(5,*) pr_main
  line='print='//pr_main
  WRITE(5,*) 'add-absorbing-poteitial'
  READ(5,*) absorb
  IF(absorb) THEN
     line=line//' add-absorbing-poteitial'
  END IF
  WRITE(5,*) 'projections'
  READ(5,*) proj
  IF(proj) THEN
     line=line//' projections'
  END IF
  WRITE(6,*) line
  stop
END PROGRAM data_input_propagator
