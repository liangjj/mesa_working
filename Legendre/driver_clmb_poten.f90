!============================================================!
! This is a driver to show how to call the Legendre ruitnes  !
! in a practical problem. The example here is the Coloumb    !
! potential between two electrons in spheroidal coordinates. !
! The potential, 1/r_12, is calculated by using Neuman's     !
! expansion in terms of regular and irregular Legendre       !
! functions in the regions (|x|<=1 or x > 1).                !
! See text for detail.                                       !
!------------------------------------------------------------!
! Input parameters are specified in a "namelist" file.       !
!------------------------------------------------------------!
! LAPACK/BLAS  : NOT used                                    !
! author       : Dr. Xiaoxu Guan at Drake University         !
! last updated : April 13, 2010                              !
!============================================================!
  PROGRAM driver_clmb_poten
  USE input_output
  USE Associated_Legendre_Functions
  USE Lentz_Thompson
  use Special_Functions
  use Data_Module,                      only : two_pi
  IMPLICIT NONE
  CHARACTER(LEN=16)                         :: Directive
  CHARACTER(LEN=24)                         :: recur

  integer                                   :: i,lorder,morder,mabs,lsum,msum,meo
  real(KIND=idp)                            :: a,radius_moeq,xi1,xi2,eta1,eta2,varphi1,varphi2
  real(KIND=idp)                            :: x1,y1,z1,rho1,x2,y2,z2,rho2,xi_small,xi_large,&
                                               facm,vardm,dl21,temp,csum_real,csum_imag,     &
                                               varphi_diff,ctemp_real,ctemp_imag,dx,dy,dz,   &
                                               rsqr,r1,r2,r_12,r_12_invs
  real(KIND=idp),dimension(:,:),allocatable :: p_lm_xismall,q_lm_xilarge,p_lm_eta1,p_lm_eta2


!xx  namelist / input_data / title, l_max, m_max, n_points, upper, lower, &
!xx                          directive, input_values, normalize, eps, recur, x
  namelist/CONSTANTS/normalize,recur,l_max,m_max,radius_moeq
  namelist/ELECTRON1/x1,y1,z1
  namelist/ELECTRON2/x2,y2,z2



!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90
!
!  Open the input and output files
!
!xx  write(6,*) '                  Lets begin the calculation'
!xx  write(6,*) '  If you wish (do not wish) some information, type yes(no)'
!xx  read(5,*) ans
!xx  IF ( ans == 'yes' ) THEN
!xx       Call Info
!xx  ELSE IF (ans == 'no' ) THEN
!xx       OPEN(inp,file='Input',status='old')
  OPEN(inp, file='coulomb_poten_r12.inp',status='old')
  OPEN(iout,file='coulomb_poten_r12.out',status='unknown')
!
!      Read all data except arguments from namelist
!
!xx      READ(inp,nml=input_data)
!
!      If you simply want to read in arbitrary values of the argument set input_values
!      to .true.  Otherwise read in an upper and lower value, a step and a number of
!      points and the code will generate the arguments.  Note that some compilers do not
!      like allocated variables to appear in namelist statements even though the 
!      allocation is done before the variable is read in.  You might have to fix that.
!      The intel compiler is fine with it.
!
!xx     IF (input_values) THEN
!xx         ALLOCATE(x(1:n_points)) 
!xx         READ(inp,nml=input_data)
!xx         upper = zero
!xx         lower = zero
!xx         DO i = 1, n_points
!xx            upper = max(upper,x(i))
!xx            lower = min(lower,x(i))
!xx         END DO
!xx     ELSE
!xx         step = (upper - lower ) / n_points
!xx         n_points = n_points + int_one
!xx         ALLOCATE(x(1:n_points)) 
!xx         x(1) = lower
!xx         DO i = 2, n_points
!xx            x(i) = x(i-int_one) + step
!xx         END DO
!xx     END IF
!xx     write(iout,1)
!xx     write(iout,2) title
!xx     write(iout,1)
!xx     write(iout,3) l_max, m_max, n_points, lower, upper, directive, &
!xx                   normalize, eps, recur
!
!    Print the arguments.
!
!xx     title='Grid'
!xx     Call Print_Matrix(x,iout)
!
!    Here we have the option of using either the Miller algorithm (A) or the continued
!    fraction approach (B).
!

     read(inp,CONSTANTS)
! check parameters:
     if(normalize/=.false. .and. normalize/=.true.) then
        write(iout,'(1x,a)') 'The parameter normalize is neither .F. nor .T. !!!'
        write(iout,'(1x,a)') 'I stopped !!!'
        STOP
     end if

     if(TRIM(recur)/='Continued_Fraction'.and.TRIM(recur)/='Miller') then
        write(iout,'(1x,a)') 'The parameter recur is neither Continued_Fraction &
     & nor Miller !!!'
        write(iout,'(1x,a)') 'I stopped !!!'
        STOP
     end if

     if(l_max < 0 .or. m_max < 0) then
        write(iout,'(1x,a)') 'The parameter l_max or m_max < 0 !!!'
        write(iout,'(1x,a,i4)') 'l_max = ',l_max
        write(iout,'(1x,a,i4)') 'm_max = ',m_max
        write(iout,'(1x,a)') 'I stopped !!!'
        STOP
     end if

     if(radius_moeq < 0.0_idp) then
        write(iout,'(1x,a)') 'The parameter radius_moeq < 0.0 !!!'
        write(iout,'(1x,a,f10.4)') 'radius_moeq = ',radius_moeq
        write(iout,'(1x,a)') 'I stopped !!!'
        STOP
     end if

     read(inp,ELECTRON1)
     read(inp,ELECTRON2)

     a = 0.5_idp*radius_moeq

!     write(*,*) normalize
!     write(*,*) recur
!     write(*,*) radius_moeq
!     write(*,*) l_max
!     write(*,*) m_max


     n_points = 1
     ALLOCATE( x(1:n_points) )

! electron 1:
     rsqr = x1*x1 + y1*y1 + z1*z1
     r1   = SQRT( a*a+rsqr-radius_moeq*z1 )
     r2   = SQRT( a*a+rsqr+radius_moeq*z1 )
     rho1 = SQRT(x1*x1+y1*y1)

     xi1  = (r1+r2)/radius_moeq
     eta1 = (r1-r2)/radius_moeq
     varphi1 = ACOS(x1/rho1)
     if(y1<0.0_idp) varphi1 = two_pi - varphi1

! electron 2:
     rsqr = x2*x2 + y2*y2 + z2*z2
     r1    = SQRT( a*a+rsqr-radius_moeq*z2 )
     r2    = SQRT( a*a+rsqr+radius_moeq*z2 )
     rho2  = SQRT(x2*x2+y2*y2)

     xi2  = (r1+r2)/radius_moeq
     eta2 = (r1-r2)/radius_moeq
     varphi2 = ACOS(x2/rho2)
     if(y2<0.0_idp) varphi2 = two_pi - varphi2


! 1/r_12:
     dx = x1-x2
     dy = y1-y2
     dz = z1-z2
     r_12      = SQRT( dx*dx + dy*dy + dz*dz )
     r_12_invs = 1.0_idp / r_12


     xi_small = MIN(xi1,xi2)
     xi_large = MAX(xi1,xi2)


!x     l_max = 20
!x     m_max = 20
!x     recur     = 'Continued_Fraction'
!     directive = 'irregular'
!     directive = 'regular'
!x     normalize = .false.

     IF ( recur == 'Miller' ) THEN
          Leg%D%A%Dir=recur
     ELSE IF ( recur == 'Continued_Fraction' ) THEN   
          Leg%D%B%Dir=recur
     END IF
!
!    Calculate either regular Legendre only, irregular Legendre only or both.
!
!xx     IF ( directive == 'regular') THEN
!xx          Call Legendre( R_LM=Leg%R_LM )
!xx     ELSE IF ( directive == 'irregular') THEN
!xx          Call Legendre( I_LM=Leg%I_LM )
!xx     ELSE IF( directive == 'both') THEN
!xx          Call Legendre( R_LM=Leg%R_LM, I_LM=Leg%I_LM  )
!xx     END IF
!xx     CLOSE(inp)
!xx     CLOSE(iout)

! P_LM(xi_small) 
     x(1) = xi_small
     ALLOCATE( p_lm_xismall(0:l_max,0:m_max) )
     directive = 'regular'
     Call Legendre( R_LM=Leg%R_LM )
     p_lm_xismall(:,:) = Leg%R_LM%F(:,:)
     DEALLOCATE( Leg%R_LM%F )

!     do morder=0,m_max
!     do lorder=morder,l_max
!     write(100,'(1x,i3,2x,i3,4x,es22.15)') lorder,morder,p_lm_xismall(lorder,morder)
!     end do
!     end do

! Q_LM(xi_large) 
     x(1) = xi_large
     ALLOCATE( q_lm_xilarge(0:l_max,0:m_max) )
     directive = 'irregular'
     Call Legendre( I_LM=Leg%I_LM )
     q_lm_xilarge(:,:) = Leg%I_LM%F(:,:)

!     do morder=0,m_max
!     do lorder=morder,l_max
!     write(100,'(1x,i3,2x,i3,4x,es22.15)') lorder,morder,q_lm_xilarge(lorder,morder)
!     end do
!     end do


! P_LM(eta1)
     x(1) = eta1
     ALLOCATE( p_lm_eta1(0:l_max,0:m_max) )
     directive = 'regular'
     Call Legendre( R_LM=Leg%R_LM )
     p_lm_eta1(:,:) = Leg%R_LM%F(:,:)
     DEALLOCATE( Leg%R_LM%F )

!     do morder=0,m_max
!     do lorder=morder,l_max
!     write(100,'(1x,i3,2x,i3,4x,es22.15)') lorder,morder,p_lm_eta1(lorder,morder)
!     end do
!     end do


! P_LM(eta2)
     x(1) = eta2
     ALLOCATE( p_lm_eta2(0:l_max,0:m_max) )
     directive = 'regular'
     Call Legendre( R_LM=Leg%R_LM )
     p_lm_eta2(:,:) = Leg%R_LM%F(:,:)
     DEALLOCATE( Leg%R_LM%F )

!     do morder=0,m_max
!     do lorder=morder,l_max
!     write(100,'(1x,i3,2x,i3,4x,es22.15)') lorder,morder,p_lm_eta2(lorder,morder)
!     end do
!     end do


! Factor(k) = k!.
      ALLOCATE( Factor(0:l_max+m_max) )
      Call Factorials 

! Summation over l and m:
         varphi_diff = varphi1 - varphi2
         csum_real = 0.0_idp
         csum_imag = 0.0_idp
      do lsum =  0, l_max
         dl21  = dble(lsum + lsum + 1)
      do msum = -lsum, lsum
         vardm  = dble(msum) * varphi_diff
         mabs = ABS(msum)
         meo  = MOD(mabs,2)
                    facm =  1.0_idp
         if(meo==1) facm = -1.0_idp

         temp = Factor(lsum-mabs) / Factor(lsum+mabs)
         temp = temp*temp
         temp = facm*dl21*temp
         temp = p_lm_xismall(lsum,mabs) * q_lm_xilarge(lsum,mabs) &
              * p_lm_eta1(lsum,mabs)    * p_lm_eta2(lsum,mabs)    &
              * temp
         ctemp_real = temp * COS(vardm)
         ctemp_imag = temp * SIN(vardm)

         csum_real = csum_real + ctemp_real
         csum_imag = csum_imag + ctemp_imag
      end do
      end do

        csum_real = csum_real / a
        csum_imag = csum_imag / a
        write(iout,'(1x,a,es22.15)') '1/r_12 (Neuman) real = ',csum_real
        write(iout,'(1x,a,es22.15)') '1/r_12 (Neuman) imag = ',csum_imag
        write(iout,'(1x,a,es22.15)') '1/r_12 (direct)      = ',r_12_invs

!xx  ELSE
!xx     stop
!xx  END IF
!

      stop
      END PROGRAM driver_clmb_poten
