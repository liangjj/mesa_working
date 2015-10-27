!deck Prolate_Driver
!**begin prologue     Prolate_Driver
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            Prolate_Driver test code
!**description        
!**                   
!**references
!**routines called        
!                       ----                    --------
!
!**end prologue       Prolate_Driver
  PROGRAM Prolate_Driver
  USE accuracy
  USE input_output
  USE Associated_Legendre_Functions
  USE Lentz_Thompson
  IMPLICIT NONE
!                                                                                                  
  TYPE (CF_Legendre)                       :: CFL
  CHARACTER(LEN=8)                         :: itoc
  INTEGER                                  :: i
  namelist / input_data / title, l_max, m_max, n_points, normalize, Print_Functions, &
                          eps, R, recur
  namelist / electron_coordinates / x1, x2, y1, y2, z1, z2
!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90
!
!  Open the input and output files
!
   OPEN(inp,file='Input_Pro',status='old')
   OPEN(iout,file='Output_Pro',status='unknown')
!
!  Read all data except arguments from namelist
!
   READ(inp,nml=input_data)
!
   write(iout,1)
   write(iout,2) title
   write(iout,1)
   write(iout,3) l_max, m_max, n_points, eps, recur
!
!     Compute factorials
!
   ALLOCATE( Factor(0:l_max+m_max) )
   Call Factorials 
!
   Write(iout,4) R
   ALLOCATE(x(1))
   DO i = 1, n_points
      READ(inp,nml=electron_coordinates)
      rsqr = x1*x1 + y1*y1 + z1*z1
      r1   = SQRT( quarter*R*R+rsqr-R*z1 )
      r2   = SQRT( quarter*R*R+rsqr+R*z1 )
      rho1 = SQRT(x1*x1+y1*y1)
      xi1  = (r1+r2)/R
      eta1 = (r1-r2)/R
      varphi1 = ACOS(x1/rho1)
      IF(y1<0.0_idp) THEN
         varphi1 = two_pi - varphi1
      END IF
      rsqr = x2*x2 + y2*y2 + z2*z2
      r1    = SQRT( quarter*R*R+rsqr-R*z2 )
      r2    = SQRT( quarter*R*R+rsqr+R*z2 )
      rho2  = SQRT(x2*x2+y2*y2)
      xi2  = (r1+r2)/R
      eta2 = (r1-r2)/R
      varphi2 = ACOS(x2/rho2)
      IF(y2<0.0_idp) THEN
         varphi2 = two_pi - varphi2
      END IF
      Write(iout,5) x1, y1, z1, x2, y2, z2, xi1, eta1, xi2, eta2
      dx = x1-x2
      dy = y1-y2
      dz = z1-z2
      r_12 = SQRT( dx*dx + dy*dy + dz*dz )
      r_12_invs = 1.0_idp / r_12
      xi_small = MIN(xi1,xi2)
      xi_large = MAX(xi1,xi2)
      IF ( recur == 'Miller' ) THEN
           Leg%D%A%Dir=recur
      ELSE IF ( recur == 'Wronskian' ) THEN   
           Leg%D%B%Dir=recur
      END IF
      x(1) = xi_small
      directive = 'regular'
      ALLOCATE(Leg%R_LM%F(0:l_max,0:m_max))
      Call Legendre( R_LM=Leg%R_LM )
      ALLOCATE(Leg%R_LM%Xi%F_Small(0:l_max,0:m_max))
      Leg%R_LM%Xi%F_Small(:,:) =  Leg%R_LM%F(:,:)
      x(1) = eta1
      Call Legendre( R_LM=Leg%R_LM )
      ALLOCATE(Leg%R_LM%Eta%F_1(0:l_max,0:m_max))
      Leg%R_LM%Eta%F_1(:,:) =  Leg%R_LM%F(:,:)
      x(1) = eta2
      Call Legendre( R_LM=Leg%R_LM )
      ALLOCATE(Leg%R_LM%Eta%F_2(0:l_max,0:m_max))
      Leg%R_LM%Eta%F_2(:,:) =  Leg%R_LM%F(:,:)
      DEALLOCATE(Leg%R_LM%F)
      x(1) = xi_large
      directive = 'irregular'
      ALLOCATE(Leg%I_LM%F(0:l_max,0:m_max))
      IF (Leg%D%B%Dir == 'Wronskian' ) THEN
          ALLOCATE(Leg%R_L%F(0:l_max))
      END IF
      Call Legendre( I_LM=Leg%I_LM )
      ALLOCATE(Leg%I_LM%Xi%F_Large(0:l_max,0:m_max))
      Leg%I_LM%Xi%F_Large(:,:) =  Leg%I_LM%F(:,:)
      DEALLOCATE(Leg%I_LM%F)
      IF (Leg%D%B%Dir == 'Wronskian' ) THEN
          DEALLOCATE(Leg%R_L%F)
      END IF
      varphi_diff = varphi1 - varphi2
      csum_real = 0.0_idp
      csum_imag = 0.0_idp
      DO lsum =  0, l_max
         dl21  = dble(lsum + lsum + 1)
         DO msum = -lsum, lsum
            vardm  = dble(msum) * varphi_diff
            mabs = ABS(msum)
            meo  = MOD(mabs,2)
            facm =  1.0_idp
            IF (meo==1) THEN
                facm = -1.0_idp
            END IF
            temp = Factor(lsum-mabs) / Factor(lsum+mabs)
            temp = temp*temp
            temp = facm*dl21*temp
            temp = temp * Leg%R_LM%Xi%F_Small(lsum,mabs) * Leg%I_LM%Xi%F_Large(lsum,mabs)   &
                                                         *                                  &
                          Leg%R_LM%Eta%F_1(lsum,mabs)    * Leg%R_LM%Eta%F_2(lsum,mabs)     
            ctemp_real = temp * COS(vardm) 
            ctemp_imag = temp * SIN(vardm)
            csum_real = csum_real + ctemp_real
            csum_imag = csum_imag + ctemp_imag
         END DO
      END DO
      csum_real = two * csum_real / R
      csum_imag = two * csum_imag / R
      Write(iout,'(1x,a,es22.15)') '1/r_12 (Neuman) real = ',csum_real
      Write(iout,'(1x,a,es22.15)') '1/r_12 (Neuman) imag = ',csum_imag
      Write(iout,'(1x,a,es22.15)') '1/r_12 (direct)      = ',r_12_invs
   END DO
   DEALLOCATE (Leg%R_LM%Xi%F_Small, Leg%I_LM%Xi%F_Large, Leg%R_LM%Eta%F_1, Leg%R_LM%Eta%F_2 )
   DEALLOCATE( Factor )
   CLOSE(inp)
   CLOSE(iout)
!
1 Format('           **************************************************************' &
         '****************')
2 Format(15x,'Begin Prolate Test Calculation = ',a80)
3 Format(/,10x,'Maximum L                             = ',i6,1x,                     &
               'Maximum M                      = ',i6,                               &
         /,10x 'Number of Points                      = ',i6,1x,                     &
         /,10x,'Continued Fraction Convergence = ',1pe15.8,1x,                       &
               'Backward Recurrence Method            = ',a24) 
4 Format(/,15x,'Testing Spheroidal Expansion at Internuclear Distance = ',1pe15.8)
5 Format(/,10x,'Coordinates: x_1  = ',1x,1pe15.8,1x,'y_1   = ',1x,1pe15.8,1x,        &
                            'z_1  = ',1x,1pe15.8,/,10x,                              &
               '             x_2  = ',1x,1pe15.8,1x,'y_2   = ',1x,1pe15.8,1x,        &
                            'z_2  = ',1x,1pe15.8,/,10x,                              &
               '             xi_1 = ',1x,1pe15.8,1x,'eta_1 = ',1x,1pe15.8,/10x,      &
               '             xi_2 = ',1x,1pe15.8,1x,'eta_2 = ',1x,1pe15.8)
  stop
END PROGRAM Prolate_Driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

