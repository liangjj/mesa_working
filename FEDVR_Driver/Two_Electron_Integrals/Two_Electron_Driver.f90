!deck Two_Electron_Driver
!**begin prologue     Two_Electron_Driver
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            Two_Electron_Driver
!**description        Computes two electron integrals using the Poisson equation
!**                   
!**references
!**routines called        
!                       ----                    --------
!
!**end prologue       
  PROGRAM Two_Electron_Driver
  USE input_output
  USE Data_Module
  USE FEDVR_Derived_Types
  USE FEDVR_Shared
  USE Two_Electron_Shared
  USE Two_Electron_FEDVR_Module
  IMPLICIT NONE
  INTEGER                          :: input
  INTEGER                          :: output
  INTEGER                          :: intkey
  INTEGER                          :: n_radial
  LOGICAL                          :: dollar
  CHARACTER (LEN = 64 )            :: chrkey 

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
!  Open the input and output files
!
  OPEN(input,file='Two_Electron.inp',status='old')
  OPEN(output,file='Two_Electron.out',status='unknown')
!
  write(iout,1)
  write(iout,2)
  write(iout,1)
  IF ( dollar('$Two_Electron_Integrals',card,cpass,inp) )THEN
       keyword = chrkey(card,'coordinate_system','spherical',' ')
       representation = chrkey(card,'representation','spherical_harmonics',' ')
       len=lenth(keyword)
       IF (keyword(1:len) /= 'spherical' .or. keyword(1:len) /= 'spheroidal') THEN
           Call lnkerr('Quit.  Coordinate System not Spherical or Spheroidal')           
       END IF
       write(iout,3) keyword
       maximum_orbital_l = intkey(card,'maximum_orbital_l',0,' ')
       maximum_orbital_m = intkey(card,'maximum_orbital_m',maximum_orbital_l,' ')
       minimum_orbital_m = - maximum_orbital_l 
       maximum_total_L = intkey(card,'maximum_total_L',maximum_orbital_l + maximum_orbital_l,' ') 
       maximum_total_M = intkey(card,'maximum_total_M',maximum_total_L,' ')  
       minimum_total_M = - maximum_total_M

       write(iout,4) 
  END IF
  Call Setup_2e 
  Call IOsys('rewind all on FEDVR_Two_Electron_Integral_File read-and-write',0,0,0,' ')
  Call IOsys('close Two_Electron_Integral_File',0,0,0,' ')
  CLOSE(input)
  CLOSE(output)
!
1 Format('           **************************************************************' &
         '****************')
2 Format(15x,'Calculation of FEDVR Two Electron Integrals')
3 FORMAT(/,25x,'Coordinate System = ',a32)
4 Format(/,15x,'Input Parameters'//                     &
               'Maximum Orbital l = ',i5,/15x,          &
               'Maximum Orbital m = ',i5,/15x,          &
               'Minimum Orbital m = ',i5,/15x,          &
               'Maximum Total L   = ',i5,/15x,          &
               'Maximum Total M   = ',i5,/15x,          &
               'Minimum Total M   = ',i5)
  stop
END PROGRAM Two_Electron_Driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
