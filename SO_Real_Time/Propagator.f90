!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE Propagator
!**begin prologue     Propagator
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       Propagator
!
                      USE Data
                      USE Derived_Types
                      USE regional_umat
  IMPLICIT NONE
!
!
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE CN_Propagation                       
                       MODULE PROCEDURE CN_Length_Propagation,                           &
                                        CN_Velocity_Propagation
                            END INTERFACE CN_Propagation
!____________________________________________________________________________________________!
!____________________________________________________________________________________________!
!
                                  Contains
!********************************************************************************
!********************************************************************************                             
!deck propagator_2_real_time
!***begin prologue     propagator_2_real_time
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            second order propagators
!***
!***references
!***routines called
!***end prologue       propagator_2
!
  SUBROUTINE propagator_2_real_time(q)
  USE dvr_global
  USE dvrprop_global
  USE regional_umat
  IMPLICIT NONE
  INTEGER                                  :: q, i, even_odd
  REAL*8                                   :: tau
!
!
! Construct the "outer" propagator at delta_t/2
!
  even_odd = num_reg(q) - 2 * ( num_reg(q) / 2 )
  IF (even_odd == 0 ) THEN
      tau = deltat*.5d0
      DO i=1,num_reg(q),2
         call umat_reg(mat_reg_d(i,q)%eigval_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%cosine_t_mat,      &
                       mat_reg_d(i,q)%sine_t_mat,        &
                       si_d,ci_d,tau,nfun_reg(i,q),1)
      END DO
!
! Construct the "inner" propagator at delta_t
!
      tau = deltat
      DO i=2,num_reg(q)-1,2
         call umat_reg(mat_reg_d(i,q)%eigval_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%cosine_t_mat,      &
                       mat_reg_d(i,q)%sine_t_mat,        &
                       si_d,ci_d,tau,nfun_reg(i,q),1)
      END DO
      i = num_reg(q)
      IF(.not.absorb) THEN
         call umat_reg(mat_reg_d(i,q)%eigval_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%cosine_t_mat,      &
                       mat_reg_d(i,q)%sine_t_mat,        &
                       si_d,ci_d,tau,nfun_reg(i,q),1)
      ELSE
         call umat_reg(mat_reg_z(q)%eigval_mat_z,        &
                       mat_reg_z(q)%eigvec_mat_z_r,      &
                       mat_reg_z(q)%eigvec_mat_z_l,      &
                       mat_reg_d(i,q)%cosine_t_mat,      &
                       mat_reg_d(i,q)%sine_t_mat,        &
                       si_z,ci_z,tau,nfun_reg(i,q),1)
      END IF
  ELSE
      tau = deltat*.5d0
      DO i=1,num_reg(q)-1,2
         call umat_reg(mat_reg_d(i,q)%eigval_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%cosine_t_mat,      &
                       mat_reg_d(i,q)%sine_t_mat,        &
                       si_d,ci_d,tau,nfun_reg(i,q),1)
      END DO
      i = num_reg(q)
      IF(.not.absorb) THEN
          call umat_reg(mat_reg_d(i,q)%eigval_mat_d,     &
                        mat_reg_d(i,q)%eigvec_mat_d,     &
                        mat_reg_d(i,q)%eigvec_mat_d,     &
                        mat_reg_d(i,q)%cosine_t_mat,     &
                        mat_reg_d(i,q)%sine_t_mat,       &
                        si_d,ci_d,tau,nfun_reg(i,q),1)
      ELSE
          call umat_reg(mat_reg_z(q)%eigval_mat_z,       &
                        mat_reg_z(q)%eigvec_mat_z_r,     &
                        mat_reg_z(q)%eigvec_mat_z_l,     &
                        mat_reg_d(i,q)%cosine_t_mat,     &
                        mat_reg_d(i,q)%sine_t_mat,       &
                        si_z,ci_z,tau,nfun_reg(i,q),1)
      END IF
!
! Construct the "inner" propagator at delta_t
!
      tau = deltat
      DO i=2,num_reg(q),2
         call umat_reg(mat_reg_d(i,q)%eigval_mat_d,      &
                         mat_reg_d(i,q)%eigvec_mat_d,    &
                         mat_reg_d(i,q)%eigvec_mat_d,    &
                         mat_reg_d(i,q)%cosine_t_mat,    &
                         mat_reg_d(i,q)%sine_t_mat,      &
                         si_d,ci_d,tau,nfun_reg(i,q),1)
      END DO
  END IF
END SUBROUTINE propagator_2_real_time
!**********************************************************************
!***********************************************************************
!deck propagator_4_real_time
!***begin prologue     propagator_4_real_time
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            fourth order propagators
!***
!***references
!***routines called
!***end prologue       propagator_4
!
  Subroutine propagator_4_real_time(q)
  USE dvr_global
  USE dvrprop_global
  USE regional_umat
  IMPLICIT NONE
  INTEGER                                  :: q, i
  REAL*8, DIMENSION(4)                     :: tau
  REAL*8                                   :: p
!
  p=1.d0
  IF (prop_order == 4) THEN
      p = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
  ENDIF
!
! Construct the propagators
!
  tau(1) = p * deltat *.5d0
  tau(2) = p * deltat
  tau(3) = ( 1.d0 - 4.d0 * p ) * deltat *.5d0
  tau(4) = ( 1.d0 - 4.d0 * p ) * deltat
  DO i=1,num_reg(q),2
     call umat_reg(mat_reg_d(i,q)%eigval_mat_d,            &
                   mat_reg_d(i,q)%eigvec_mat_d,            &
                   mat_reg_d(i,q)%eigvec_mat_d,            &
                   mat_reg_d(i,q)%cosine_t_mat,            &
                   mat_reg_d(i,q)%sine_t_mat,              &
                   si_d,ci_d,tau(1),nfun_reg(i,q),1)
     call umat_reg(mat_reg_d(i,q)%eigval_mat_d,            &
                   mat_reg_d(i,q)%eigvec_mat_d,            &
                   mat_reg_d(i,q)%eigvec_mat_d,            &
                   mat_reg_d(i,q)%cosine_t_mat,            &
                   mat_reg_d(i,q)%sine_t_mat,              &
                   si_d,ci_d,tau(2),nfun_reg(i,q),2)
     call umat_reg(mat_reg_d(i,q)%eigval_mat_d,            &
                   mat_reg_d(i,q)%eigvec_mat_d,            &
                   mat_reg_d(i,q)%eigvec_mat_d,            &
                   mat_reg_d(i,q)%cosine_t_mat,            &
                   mat_reg_d(i,q)%sine_t_mat,              &
                   si_d,ci_d,tau(4),nfun_reg(i,q),3)
  END DO
  DO i=2,num_reg(q),2
     call umat_reg(mat_reg_d(i,q)%eigval_mat_d,            &
                   mat_reg_d(i,q)%eigvec_mat_d,            &
                   mat_reg_d(i,q)%eigvec_mat_d,            &
                   mat_reg_d(i,q)%cosine_t_mat,            &
                   mat_reg_d(i,q)%sine_t_mat,              &
                   si_d,ci_d,tau(2),nfun_reg(i,q),1)
     call umat_reg(mat_reg_d(i,q)%eigval_mat_d,            &
                   mat_reg_d(i,q)%eigvec_mat_d,            &
                   mat_reg_d(i,q)%eigvec_mat_d,            &
                   mat_reg_d(i,q)%cosine_t_mat,            &
                   mat_reg_d(i,q)%sine_t_mat,              &
                   si_d,ci_d,tau(4),nfun_reg(i,q),2)
  END DO
END SUBROUTINE   Subroutine propagator_4_real_time
!**********************************************************************
!**********************************************************************
!deck propagator_real_time
!***begin prologue     propagator_real_time
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            propagators
!***
!***description        The second order decomposition of the time exponential
!***                   is performed for a Hamiltonian H = H_a + H_b + V as
!***                   U_2(tau) = exp(-V * tau/2)   *
!***                                   u_2(tau)     *
!***                              exp(-V * tau/2)
!***
!***                   u_2(tau) = exp(-H_a * tau/2) * 
!***                              exp(-H_b * tau)   * 
!***                              exp(_H_a * tau/2) *
!***
!***                   The expression for u_2 requires the calculation of two 
!***                   independent propagators.
!***
!***                   The fourth order decomposition of the time exponential
!***                   is performed as
!***                   U_4(tau) = U_2(p*tau)       * 
!***                              U_2(p*tau)       * 
!***                              U_2((1-4*p)*tau) * 
!***                              U_2(p*tau)       * 
!***                              U_2(p*tau) 
!***                              p = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
!***
!***                   The expression for U_4 requires the calculation of four 
!***                   independent propagators.
!***references
!***routines called
!***end prologue       propagator_real_time
!
  SUBROUTINE propagator_real_time(prop,mat_reg,n_fun,num_reg)
  IMPLICIT NONE
  TYPE(REGIONAL_MATRICES),   DIMENSION(:)  :: mat_reg
  TYPE(COMPLEX_PROP)                       :: prop
  INTEGER, DIMENSION(:)                    :: n_fun                                  
  INTEGER                                  :: q, i, even_odd
  INTEGER                                  :: trips, ntrip
  REAL(idp)                                :: tau
!
!
! Special case of only one region
!
  p=1.d0
  trips = 1
  IF (prop_order == 4) THEN
      p = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
      trips = 2
  ENDIF
!
  tau(1) = p * deltat * .5d0
  tau(2) = 2.d0 * tau(1)
  tau(3) = ( 1.d0 - 4.d0 * p ) * deltat * .5d0
  tau(4) = 2.d0 * tau(3)
!
  write(iout,1) prop_order, p
!
  IF ( num_reg == 1) THEN
       write(iout,2)
       call umat_reg(mat_reg(1)%eigval_mat,             &
                     mat_reg(1)%eigvec_mat,             &
                     mat_reg(1)%exp_t_mat(:,:,1)        &
                     exp_tmp,deltat,nfun_reg(1))
  ELSE
      Write(iout,3) trips
!
      DO ntrip=1,trips
         write(iout,4) ntrip, tau(ntrip:ntrip+1)
!
!         Construct the "outer" propagator at p_fac * delta_t/2
!
          DO i=1,num_reg(q),2
             call umat_reg(mat_reg(i)%eigval_mat,             &
                           mat_reg(i)%eigvec_mat,             &
                           mat_reg(i)%exp_t_mat(:,:,ntrip),   &
                           exp_tmp,tau(ntrip),nfun_reg(i))
          END DO
!
!         Construct the "inner" propagator at p_loc * delta_t
!
          DO i=2,num_reg(q),2
             call umat_reg(mat_reg(i)%eigval_mat,             &
                           mat_reg(i)%eigvec_mat,             &
                           mat_reg(i)%exp_t_mat(:,:,ntrip),   &
                           exp_tmp,tau(ntrip+1),nfun_reg(i))
          END DO
!
       END DO
  END IF
1 Format(/,10x,'Constructing the Propagators for order = ',i2, &
         /,10x,'Exponential factor                     = ',e15.8)
2 FORMAT(/,10x,'There is only one region')
  FORMAT(/,10x,'Number of Propagator Passes            = ',i1)
4 FORMAT(/,10x,'Pass                                   = ',i1, &
         /,10x,'p_fac                                  = ',e15.8)
END SUBROUTINE propagator_real_time
!**********************************************************************
!**********************************************************************
!deck propagator_imaginary_time
!***begin prologue     propagator_imaginary_time
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            propagators
!***
!***description        The second order decomposition of the time exponential
!***                   is performed for a Hamiltonian H = H_a + H_b + V as
!***                   U_2(tau) = exp(-V * tau/2)   *
!***                                   u_2(tau)     *
!***                              exp(-V * tau/2)
!***
!***                   u_2(tau) = exp(-H_a * tau/2) * 
!***                              exp(-H_b * tau)   * 
!***                              exp(_H_a * tau/2) *
!***
!***                   The expression for u_2 requires the calculation of two 
!***                   independent propagators.
!***
!***                   The fourth order decomposition of the time exponential
!***                   is performed as
!***                   U_4(tau) = U_2(p*tau)       * 
!***                              U_2(p*tau)       * 
!***                              U_2((1-4*p)*tau) * 
!***                              U_2(p*tau)       * 
!***                              U_2(p*tau) 
!***                              p = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
!***
!***                   The expression for U_4 requires the calculation of four 
!***                   independent propagators.
!***references
!***routines called
!***end prologue       propagator_imaginary_time
!
  SUBROUTINE propagator_imaginary_time(prop,mat_reg,n_fun,num_reg)
  IMPLICIT NONE
  TYPE(REGIONAL_MATRICES),   DIMENSION(:)  :: mat_reg
  TYPE(REAL_PROP)                          :: prop
  INTEGER, DIMENSION(:)                    :: n_fun                                  
  INTEGER                                  :: q, i, even_odd
  INTEGER                                  :: trips, ntrip
  REAL(idp)                                :: tau
!
!
! Special case of only one region
!
  p=1.d0
  trips = 1
  IF (prop_order == 4) THEN
      p = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
      trips = 2
  ENDIF
!
  tau(1) = p * deltat * .5d0
  tau(2) = 2.d0 * tau(1)
  tau(3) = ( 1.d0 - 4.d0 * p ) * deltat * .5d0
  tau(4) = 2.d0 * tau(3)
!
  write(iout,1) prop_order, p
!
  IF ( num_reg == 1) THEN
       write(iout,2)
       call umat_reg(mat_reg(1)%eigval_mat,             &
                     mat_reg(1)%eigvec_mat,             &
                     mat_reg(1)%exp_t_mat(:,:,1)        &
                     exp_tmp,deltat,nfun_reg(1))
  ELSE
      Write(iout,3) trips
!
      DO ntrip=1,trips
         write(iout,4) ntrip, tau(ntrip:ntrip+1)
!
!         Construct the "outer" propagator at p_fac * delta_t/2
!
          DO i=1,num_reg(q),2
             call umat_reg(mat_reg(i)%eigval_mat,             &
                           mat_reg(i)%eigvec_mat,             &
                           mat_reg(i)%exp_t_mat(:,:,ntrip),   &
                           exp_tmp,tau(ntrip),nfun_reg(i))
          END DO
!
!         Construct the "inner" propagator at p_loc * delta_t
!
          DO i=2,num_reg(q),2
             call umat_reg(mat_reg(i)%eigval_mat,             &
                           mat_reg(i)%eigvec_mat,             &
                           mat_reg(i)%exp_t_mat(:,:,ntrip),   &
                           exp_tmp,tau(ntrip+1),nfun_reg(i))
          END DO
!
       END DO
  END IF
1 Format(/,10x,'Constructing the Propagators for order = ',i2, &
         /,10x,'Exponential factor                     = ',e15.8)
2 FORMAT(/,10x,'There is only one region')
  FORMAT(/,10x,'Number of Propagator Passes            = ',i1)
4 FORMAT(/,10x,'Pass                                   = ',i1, &
         /,10x,'p_fac                                  = ',e15.8)
END SUBROUTINE propagator_imaginary_time
!**********************************************************************
!**********************************************************************
!deck propagator_real_time
!***begin prologue     propagator_real_time
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            propagators
!***
!***description        The second order decomposition of the time exponential
!***                   is performed for a Hamiltonian H = H_a + H_b + V as
!***                   U_2(tau) = exp(-V * tau/2)   *
!***                                   u_2(tau)     *
!***                              exp(-V * tau/2)
!***
!***                   u_2(tau) = exp(-H_a * tau/2) * 
!***                              exp(-H_b * tau)   * 
!***                              exp(_H_a * tau/2) *
!***
!***                   The expression for u_2 requires the calculation of two 
!***                   independent propagators.
!***
!***                   The fourth order decomposition of the time exponential
!***                   is performed as
!***                   U_4(tau) = U_2(p*tau)       * 
!***                              U_2(p*tau)       * 
!***                              U_2((1-4*p)*tau) * 
!***                              U_2(p*tau)       * 
!***                              U_2(p*tau) 
!***                              p = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
!***
!***                   The expression for U_4 requires the calculation of four 
!***                   independent propagators.
!***references
!***routines called
!***end prologue       propagator_real_time
!
  SUBROUTINE propagator(q)
  USE dvr_global
  USE dvrprop_global_it
  IMPLICIT NONE
  INTEGER                                  :: q, i, even_odd
  INTEGER                                  :: trips, ntrip
  REAL*8                                   :: tau
!
!
! Special case of only one region
!
  write(iout,1) prop_order, p_fac
!
  IF ( num_reg == 1) THEN
       write(iout,2)
       IF(prop_order==2) THEN
          tau = deltat
          call umat_reg(mat_reg%eigval_mat,             &
                        mat_reg%eigvec_mat,             &
                        mat_reg%exp_t_mat(:,:,1),       &
                        exp_tmp,tau,nfun_reg(1))
       ELSE IF(prop_order==4) THEN
          p_loc = p_fac
          trips = 2
          Write(iout,3) trips
          DO ntrip=1,trips
             write(iout,4) ntrip, p_loc
              tau = p_loc * deltat
              call umat_reg(mat_reg%eigval_mat,           &
                            mat_reg%eigvec_mat,           &
                            mat_reg%exp_t_mat(:,:,ntrip), &
                            exp_tmp,tau,nfun_reg(1))
!                                                                          
!             Modify p_loc for fourth order propagator so it is 
!             correct second time through last loop.            
!                                                                          
              p_loc = ( 1.d0 - 4.d0 * p_loc )
          END DO
       END IF
  ELSE
       p_loc = p_fac
       trips = 1
       IF ( prop_order == 4 ) THEN
            trips = 2
       END IF
       Write(iout,3) trips
!
       DO ntrip=1,trips
          write(iout,4) ntrip, p_loc
!
!         Construct the "outer" propagator at p_fac * delta_t/2
!
          tau = p_loc * deltat *.5d0
          DO i=1,num_reg,2
             call umat_reg(mat_reg(i)%eigval_mat,             &
                           mat_reg(i)%eigvec_mat,             &
                           mat_reg(i)%exp_t_mat(:,:,ntrip),   &
                           exp_tmp,tau,nfun_reg(i))
          END DO
!
!         Construct the "inner" propagator at p_loc * delta_t
!
          tau = p_loc * deltat
          DO i=2,num_reg,2
             call umat_reg(mat_reg(i)%eigval_mat,             &
                           mat_reg(i)%eigvec_mat,             &
                           mat_reg(i)%exp_t_mat(:,:,ntrip),     &
                           exp_tmp,tau,nfun_reg(i))
          END DO
!
!         Modify p_loc for fourth order propagator so its correct second time 
!         through last loop.
!
          p_loc = ( 1.d0 - 4.d0 * p_loc )
       END DO
  END IF
1 Format(/,10x,'Constructing the Propagators for order = ',i2, &
         /,10x,'Exponential factor                     = ',e15.8)
2 FORMAT(/,10x,'There is only one region')
3 FORMAT(/,10x,'Number of Propagator Passes            = ',i1)
4 FORMAT(/,10x,'Pass                                   = ',i1, &
         /,10x,'p_fac                                  = ',e15.8)
END SUBROUTINE propagator_real_time
!**********************************************************************
!**********************************************************************
END MODULE Propagator

