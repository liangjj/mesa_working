!deck Input_DVR
!***begin prologue     Input_DVR
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           Input, DVR
!***author             schneider, b. i.(nsf)
!***source
!***purpose            input for DVR
!***
!***references
!***routines called    iosys, util and mdutil
!***modules used
!**                    Name                 Location
!                      ----                 -------
!
!                   dvrprop_global          Modules
!                   dvr_shared              Modules
!                   dvr_global              Modules
!
!***end prologue       input
!
!
  Subroutine Input_DVR(case)
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  USE potential
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  LOGICAL                                  :: dollar, logkey
  REAL*8                                   :: fpkey
  INTEGER                                  :: intkey, i, j
  INTEGER                                  :: case
!

  WRITE(iout,*)
!
!  Look for the keyword $dvr_basis in the input file and read all variables
!  that are associated with this keyword.
!
  IF ( dollar('$dvr_basis_case_'//itoc(case),card,cpass,inp) ) then
!
!         Set spatial dimensionality of problem and program 
!         options
!
        keywrd = chrkey(card,'key_word','h0',' ')
!       Is this a 1D, 2D or 3D calculation.
!
        spdim=intkey(card,'number_of_space_variables',1,' ')
        WRITE(iout,1) spdim
        system=chrkey(card,'coordinate_system','cartesian',' ')
!
!       The Module grid_global.f90 defines a whole lot of constants that are
!       used in the calculation.  Are we in atomic units?  If not then the
!       definitions of the constants give hbar in joule-sec and mass of the
!       electron in kilograms.
!
        potential_type = 'short_range'
        units=chrkey(card,'units','atomic_units',' ')
        IF(units == 'atomic_units') then
           hbar=1.d0
           mass=1.d0
        END IF
!
!       To automate the way points are generated, set genpts to .true.
!       How this is done is defined elsewhere.
!
!
!       The parameters below need to be set or the default values
!       will be used.  In most cases the default can be used as there
!       is typically no non-linear potential, the basis is DVR (fd is the
!       other possibility and the convergence criterion is tight enough)
!       
        diag=.true.
        typke='dvr'
        proj=.false.
        write(iout,2) system, units
!
!       Print options are set here.  Taking the default(none) will give minimal
!       print which is ok most of the time but not when debugging or getting
!       a feel for how the code runs.  Using the all-details option gives lots
!       of information.
!
!       One can modify the diagonal part of the one-body operators to
!       include the diagonal part of the one-body potential.  This is not
!       always important but when you are splitting the operator could result
!       in more accuracy to second order.
!
        diag_mod=chrkey(card,'diagonal_modification','none',' ')
!
!       Now we begin the real work.
!       For each dimension compute the minimal information to proceed.
!
        DO i=1, spdim
!
!          This is just a label to get to a keyword for data entry.
!
           coord(i)=chrkey(card,'space_variable_'//itoc(i),'x',' ')
        END DO   
  ELSE
        write(iout,3)
        stop
  END IF  
!
1 FORMAT(/,20X,'Set Up FEDVR Basis and Matrices',/20x,           &
               'Number of Variables = ',i1)
2 FORMAT(/,25X,'data',                                           &
         /,5X, 'coordinate system          = ',a32,              &
         /,5X,'units                      = ',a32)
3 FORMAT(/,5x,'no basis card section')
END Subroutine Input_Dvr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
