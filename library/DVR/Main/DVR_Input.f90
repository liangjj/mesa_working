!************************************************************************************
!deck Dvr_Input
!***begin prologue     Dvr_Input
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for spheroidal dvr basis sets
!***description        sets up the unperturbed dv hamiltonian for a two center
!***                   problem where charge Z_a is at -R/2 and Z_b at +R/2
!***references
!***routines called    
!***end prologue       Dvr_Input
  SUBROUTINE DVR_Input
  USE FEDVR_Global
  IMPLICIT NONE
  CHARACTER(LEN=80)                :: chrkey
  INTEGER                          :: intkey
  LOGICAL                          :: dollar
  LOGICAL                          :: logkey
  REAL(idp)                        :: fpkey
  INTEGER                          :: input
  INTEGER                          :: output
  INTEGER                          :: i
  COMMON /io/ input, output
!                                                                                                                                                   
!  Get the input and output file numbers which appear in the                                                                                        
!  Module input_output.f90 and put them into input and output which appears                                                                         
!  in the ONLY common block in the code.  This is needed to pass into library                                                                       
!  routines and for no other purpose.                                                                                                               
!                                                                                                                                                   
  input = inp
  output = iout
!                                                                                                                                                  
! ppen the input and output files                                                                                                             
!                                                                                                                                                   
  OPEN(input,file='dvr.inp',status='old')
  OPEN(output,file='dvr.out',status='unknown')
!
  IF ( dollar('$dvr_input',card,cpass,inp) ) THEN
!
!     Set the dvr keyword to the coordinate system.
!
       keyword = chrkey(card,'dvr_keyword','cartesian',' ')
!
!      Set the number of dimensions and labels.
!
       spdim = intkey(card,'number_of_spatial_dimensions',1,' ')
       ALLOCATE (reg_grid(1:spdim))
       DO i=1,spdim
          reg_grid(i)%label = chrkey(card,'coordinate_label_'//itoc(i),'x',' ')
       END DO
!
!      set the print variables for printing.
!
      prloc = 'print='//prnkey
      prloc(12)=chrkey(card,'print',prnkey(12),' ')
      IF(prloc(12) == 'all') THEN
         prn=.true.
      ELSE
         CALL setlog(prn,prloc,card,11)
      END IF
      mass=fpkey(card,'mass',massau,' ')
      units=chrkey(card,'units','atomic_units',' ')
      parity=chrkey(card,'parity','none',' ')
      IF(units == 'atomic_units') THEN
         mass=mass/massau
      END IF
      nodiag=logkey(card,'do_not_diagonalize',.false.,' ')
      diag=logkey(card,'diagonalize',.true.,' ')
      typwt=chrkey(card,'weight_type','one',' ')
      unit_weight=logkey(card,'use_unit_weight',.false.,' ')
      atomic=logkey(card,'atomic_calculation',.false.,' ')
      typke='dvr'
      drctv=chrkey(card,'type_calculation','all_integrals',' ')
      if(drctv == 'poisson') then
         dentyp=chrkey(card,'type_density','exponential',' ')
      end if
      proj=.false.
  END IF
!
  WRITE(iout,1) keyword
  IF (keyword == 'cartesian') THEN
      DO i = 1, spdim
         Call Read_Cartesian(reg_grid(i))
         Call Coordinate_Factors(reg_grid(i))
         Call Lobatto_Functions (reg_grid(i))
         Call FE_DVR_Matrices(reg_grid(i))
         Call Pe_Fedvr(reg_grid(i))
         Call H_0_Fedvr(reg_grid(i))
      END DO
  ELSE IF(keywrd == 'spherical') then
      Call Read_Spherical_Data
      DO i = 1, spdim
         Call Read_Spherical(reg_grid(i))
         Call Coordinate_Factors(reg_grid(i))
         Call Lobatto_Functions (reg_grid(i))
         Call FE_DVR_Matrices(reg_grid(i))
         Call Pe_Fedvr(reg_grid(i))
         Call H_0_Fedvr(reg_grid(i))
      END DO
  ELSE IF(keywrd == 'fourier') then
         Call Read_Fourier(reg_grid(1))
         Call Lobatto_Functions (reg_grid(i))
         Call FE_DVR_Matrices(reg_grid(1))
  ELSE IF (typwt == 'legendre') then
           Call Read_Legendre(reg_grid(1))
  ELSE IF (typwt == 'theta') then
           Call Read_Theta(reg_grid(1))
  ELSE IF (typwt == 'hermite') then
          Call Read_Hermite(reg_grid(1))
           Call Lobatto_Functions (reg_grid(i))
  ELSE IF(typwt == 'laguerre') then
           Call Read_LaGuerre(reg_grid(1))
           Call Lobatto_Functions (reg_grid(i))
  ELSE IF (keywrd == 'spheroidal') THEN
           Call Read_Spheroidal_Data
           DO i = 1, spdim
              Call Read_Spheroidal(reg_grid(i))
              Call Coordinate_Factors(reg_grid(i))
              Call Lobatto_Functions (reg_grid(i))
              Call FE_DVR_Matrices(reg_grid(i))
              Call Pe_Fedvr(reg_grid(i)%label,reg_grid(i))
              Call H_0_Fedvr(reg_grid(i))
           END DO
  ELSE
       Call lnkerr('bad keyword')
       write(iout,2)
  END IF
1  FORMAT(/,40X,'dvr functions and operators in ',a16,               &
                ' coordinate system'                                 &
2 FORMAT(/,1x,'no input found corresponding to keyword = ',a80)
  END SUBROUTINE DVR_Input
!*************************************************************************************
