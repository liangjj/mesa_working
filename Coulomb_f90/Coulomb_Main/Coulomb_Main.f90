!deck @(#)Coulomb_Main
!***begin prologue     Coulomb_Main
!***date written       920525   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           
!***author             schneider, b. (nsf)
!***source             
!***purpose            positive energy regular and irregular coulomb functions
!***
!***description        regular and irregular coulomb function are computed
!***                   using the Barnett code.
!***references

!***routines called    iosys, util and mdutil
!***end prologue       Coulomb_Main

  Program Coulomb_Main
  USE Coulomb_Module
  IMPLICIT NONE
  LOGICAL                                  :: dollar
  REAL*8                                   :: fpkey
  INTEGER                                  :: input
  INTEGER                                  :: output
  INTEGER                                  :: intkey
  INTEGER                                  :: i
  INTEGER                                  :: j
  LOGICAL                                  :: logkey
  CHARACTER (LEN=80)                       :: chrkey
  COMMON /io/ input, output
!
  input = inp
  output = iout

!
!  Open the input and output files
!
  OPEN(input,file='Coulomb.inp',status='old')
  OPEN(output,file='Coulomb.out',status='unknown')
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  IF ( dollar('$coulomb',card,cpass,inp) ) then
!
!**********************************************************************c
!                calculate k or kappa and the value of the             c
!                             coulomb eta                              c
!**********************************************************************c
     charge = fpkey(card,'charge',-1.d0,' ')
     energy = fpkey(card,'energy',1.d0,' ')
     l_min = intkey(card,'smallest_angular_momentum',0,' ')
     l_max = intkey(card,'largest_angular_momentum',0,' ')
     r_min = fpkey(card,'smallest_r',0.d0,' ')
     r_max = fpkey(card,'largest_r',0.d0,' ')
     r_step = fpkey(card,'delta_r',0.d0,' ')
     number_of_r_values = ( r_max - r_min ) / r_step + 1
     quantities_returned = chrkey(card,'quantities_returned',          &
                                  'functions_and_derivatives',' ')
     type = chrkey(card,'function_type','coulomb',' ')
     ALLOCATE(rho(1:number_of_r_values))
     ALLOCATE(fc(1:l_max+1), dfc(1:l_max+1), gc(1:l_max+1), dgc(1:l_max+1))
     call fparr(card,'r_values',rho,number_of_r_values,' ')
     rho(1) = r_min
     DO i= 2, number_of_r_values
        rho(i) = rho(i-1) + r_step
     END DO
     k = SQRT(energy*two)
     eta_in = charge/k
     rho(:) = k * rho(:)
     WRITE(iout,3) energy, k, l_min, l_max, eta_in
     DO i = 1, number_of_r_values
        Call Coulfg (rho(i),eta_in,fc,gc,dfc,dgc) 
        IF (ifail /= 0) THEN
            write(iout,*) 'Warning message !!!'
        END IF
        write(iout,4) rho(i)
        DO j = l_min+1 , l_max+1
           fl  = fc(j)
           dfl = dfc(j)
           gl  = gc(j)
           dgl = dgc(j)
           wronskian = fl*dgl - dfl*gl
           write(iout,5) j-1, fl, dfl , gl , dgl, wronskian
        END DO
     END DO
  END IF
  WRITE(iout,1)
1 FORMAT('************************************************************************')
2 FORMAT(25X,'Coulomb Functions')
3 FORMAT(/,5X,'Energy                    = ',e15.8,1X,'k(kappa)      = ',e15.8,  &
         /,5x,'Smallest Angular Momentum = ',i3,                                 &
           1x,'Largest Angular Momentum  = ',i3,/,5x,'Eta  = ',e15.8)
4 FORMAT(/,30x,'rho = ',d15.8)
5 FORMAT(/,5x,'L = ',i4,2x,'R = ',e15.8,2x,'DR = ',e15.8,2x,'G = ',e15.8,2x,'DG = ',e15.8,  &
         /,15x,'Wronskian = ',e15.8)
  END PROGRAM Coulomb_Main
