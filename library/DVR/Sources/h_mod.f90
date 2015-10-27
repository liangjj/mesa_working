!***********************************************************************
! FEDVR_Hamiltonian_Module
!**begin prologue     FEDVR_Hamiltonian_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Calculate the kinetic and potential energy matrix elements
!***                  in a FEDVR basis
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      FEDVR_Hamiltonian_Module
!***********************************************************************
!***********************************************************************
                           MODULE FEDVR_Hamiltonian_Module
                           USE FEDVR_Global
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE Kinetic_Energy                       
                       MODULE PROCEDURE KE,                             &
                                        Odd_KE,                         &
                                        Fourier_KE,                     &  
                                        Hermite_KE,                     &  
                                        Laguerre_KE  
                            END INTERFACE Kinetic_Energy
!
                            INTERFACE Matrix_Renormalization                                     
                       MODULE PROCEDURE Mat_Renormalization,            &
                                        Odd_Mat_Renormalization    
                            END INTERFACE Matrix_Renormalization                                     
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck KE.f
!***begin prologue     KE
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the parts of the regional matrices
!***                   not depending on the m quantum number.  since a different
!***                   form for DVR basis functions are needed for even and odd
!***                   m, two different kinetic energy operators need to be constructed.
!***                   Later, the even and odd m kinetic energy will be built from these.
!***                   Also, the formulas programmed assume that 1) the kinetic energy
!***                   is left in its standard form and no first derivative has been traansformed
!***                   away and 2) the volume element is included in the definition of the the matrix
!***                   element.  This means that if there is a singular point at the boundaries, the
!***                   integrand remains well behaved and it is not necessary to make that a quadrature
!***                   point.  This allows us to use Gauss-Radau rules for those elements instead of
!***                   Gauss-Lobatto and you do not have to remove the first basis function.  This was
!***                   first pointed out to me by Brett Esry.
!***references

!***routines called    iosys, util and mdutil
!***end prologue       Ke
  SUBROUTINE KE(grid,reg_mat,reg_poly)
  IMPLICIT NONE
  TYPE(coordinates)                    :: grid  
  TYPE(matrices), DIMENSION(:)         :: reg_mat  
  TYPE(functions), DIMENSION(:)        :: reg_poly  
  INTEGER                              :: ir
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
!
!
  DO ir = 1, nreg
     ALLOCATE( grid%reg_mat(ir)%tr( npt(ir), npt(ir) ) )
     grid%reg_mat(ir)%tr(:,:) = zero
     DO i = 1, npt(ir)
        DO j = 1, i
           DO k = 1, npt(ir)
              grid%reg_mat(ir)%tr(i,j) = grid%reg_mat(ir)%tr(i,j)          &
                                           -                               &
              grid%reg_pt_wt(ir)%qr_fac(k) * grid%reg_pt_wt(ir)%wtr(k)     &
                                           *                               &
              grid%reg_poly(ir)%dpr(k,i)   * grid%reg_poly(ir)%dpr(k,j) 
           END DO
           grid%reg_mat(ir)%tr(j,i) = grid%reg_mat(ir)%tr(i,j)
        END DO
     END DO
  END DO
  Call Matrix_Renormalization(grid,reg_mat)
!
END SUBROUTINE KE
!***********************************************************************
!***********************************************************************
!deck Odd_KE.f
!***begin prologue     Odd_KE
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the parts of the regional matrices
!***                   not depending on the m quantum number.  since a different
!***                   form for DVR basis functions are needed for even and odd
!***                   m, two different kinetic energy operators need to be constructed.
!***                   Later, the even and odd m kinetic energy will be built from these.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Odd_KE

  SUBROUTINE Odd_KE(grid,reg_mat_odd,reg_poly_odd)
  IMPLICIT NONE
  TYPE(coordinates)                    :: grid  
  TYPE(odd_functions), DIMENSION(:)    :: reg_poly_odd  
  TYPE(odd_matrices), DIMENSION(:)     :: reg_mat_odd  
  INTEGER                              :: ir
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
!
!
  DO ir = 1, nreg
     ALLOCATE( grid%reg_mat_odd(ir)%tr( npt(ir), npt(ir) ) )
     grid%reg_mat_odd(ir)%tr(:,:) = zero
     DO  i = 1, npt(ir)
         DO j = 1, i
            DO k = 1, npt(ir)
               grid%reg_mat_odd(ir)%tr(i,j) = grid%reg_mat_odd(ir)%tr(i,j)     &
                                          -                                    &
               grid%reg_pt_wt(ir)%qr_fac(k) * grid%reg_pt_wt(ir)%qr_fac(k)     &
                                            *                                  &
               grid%reg_pt_wt(ir)%wtr(k)    * grid%reg_poly_odd(ir)%dpr(k,i)   &
                                            * grid%reg_poly_odd(ir)%dpr(k,j)
            END DO
            grid%reg_mat_odd(ir)%tr(i,j) = grid%reg_mat_odd(ir)%tr(i,j)        &
                                      -                                        &
            pre_factor * ( grid%reg_pt_wt(ir)%qr(i)                            &
                                      *                                        &
                           grid%reg_pt_wt(ir)%qr_fac(i)                        &
                                      *                                        &
                           grid%reg_poly_odd(ir)%pr(i,i)                       &
                                      *                                        &
                           grid%reg_poly_odd(ir)%dpr(i,j)                      &
                                      +                                        &
                           grid%reg_pt_wt(ir)%qr(j)                            &
                                      *                                        &
                           grid%reg_pt_wt(ir)%qr_fac(j)                        &
                                      *                                        &
                           grid%reg_poly_odd(ir)%pr(j,j)                       &
                                      *                                        &
                           grid%reg_poly_odd(ir)%dpr(j,i) ) 
            grid%reg_mat_odd(ir)%tr(j,i) = grid%reg_mat_odd(ir)%tr(i,j)
         END DO
         grid%reg_mat_odd(ir)%tr(i,i) = grid%reg_mat_odd(ir)%tr(i,i)           &
                                     -                                         &
         grid%reg_pt_wt(ir)%qr(i)    * grid%reg_pt_wt(ir)%qr(i)                &
                                     *                                         &
         grid%reg_poly_odd(ir)%pr(i,i) * grid%reg_poly_odd(ir)%pr(i,i)
     END DO
  END DO
  Call Matrix_Renormalization(grid,reg_mat_odd)
!
END SUBROUTINE Odd_KE
!***********************************************************************
!***********************************************************************
!***deck Fourier_KE.f90
!***begin prologue     Fourier_KE
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)                                                               
!***keywords           kinetic energy, hamiltonian, fourier basis                                                    
!***author             schneider, barry (nsf)                                                                        
!***source                                                                                                           
!***purpose            generate matrix elements                                                                      
!***                   for a fourier expansion. only a single element                                                                    !***                   is considered. there is no integration by parts and
!***                   what is produced is the discretized second derivative
!***                   matrix. 
!***references                                                                                                       
  SUBROUTINE Fourier_KE(grid,reg_mat_fourier)
  IMPLICIT NONE
  TYPE(coordinates)                      :: grid  
  TYPE(fourier_matrices), DIMENSION(:)   :: reg_mat_fourier  
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: k
  INTEGER                                :: l
  INTEGER                                :: m
  INTEGER                                :: delim
  INTEGER                                :: pre
  INTEGER                                :: ifac
  REAL(idp)                              :: kii
  REAL(idp)                              :: kij
  REAL(idp)                              :: cosfac
  REAL(idp)                              :: sinfac_2
  CHARACTER (LEN=80)                     :: title
  ALLOCATE( grid%reg_mat_fourier(1) )
  ALLOCATE( grid%reg_mat_fourier(1)%tr( npt(1), npt(1) ) ) 
  j=( npt(1) - 1 )/2
  kii=-pi*pi*( npt(1)*npt(1) - one)/(three*box*box)
  kij=-two*pi*pi/(box*box)
  DO  i=1,npt(1)
      grid%reg_mat_fourier(1)%tr(i,i) = kii
  END DO
  k=-j
  DO i=1,npt(1)
     l=-j
     DO m=1,i-1
        delim=k-l
        cosfac=cos(pi*delim/npt(1))
        sinfac_2=sin(pi*delim/npt(1))
        sinfac_2=sinfac_2 * sinfac_2
        pre=1
        ifac=delim - 2*(delim/2)
        if(ifac == 1) THEN
           pre=-1
        END IF
        grid%reg_mat_fourier(1)%tr(i,m) = kij * pre * cosfac / sinfac_2
        l = l + 1
     END DO
     k = k + 1
  END DO
  DO i=1,npt(1)
     DO j=1,i
        grid%reg_mat_fourier(1)%tr(j,i) = grid%reg_mat_fourier(1)%tr(i,j)
     END DO
  END DO
1 Format(/,5x,'fourier weight = ',e15.8)
END SUBROUTINE Fourier_KE
  Call Pe_Fedvr(grid)
  grid%reg_mat_fourier(1)%tr(:,:) = dscale * grid%reg_mat_fourier(1)%tr(:,:)
  IF(prn(3)) THEN
     title='normalized kinetic energy matrix for fourier basis'
     CALL prntrm(title,grid%reg_mat_fourier(1)%tr,npt(1),npt(1),npt(1),npt(1),iout)
  END IF
  ALLOCATE( grid%reg_mat(1)%ham( npt(1), npt(1) ) ) 
  grid%reg_mat(1)%ham(:,:) =  grid%reg_mat(1)%tr(:,:)
  DO i=1,npt(1)
     grid%reg_mat(1)%ham(i,i) =  grid%reg_mat(1)%ham(i,i) + grid%reg_vec(1)%vr(i)
  END DO
  IF(prn(3)) THEN
     title='normalized hamiltonian matrix for fourier basis'
     CALL prntrm(title,grid%reg_mat(1)%ham,npt(1),npt(1),npt(1),npt(1),iout)
  END IF
  IF(.not.nodiag) THEN
     ALLOCATE( grid%reg_mat(1)%eigval(npt(1)),          &
               grid%reg_mat(1)%eigvec(npt(1),npt(1)),   &
               grid%reg_mat(1)%eigval_0(npt(1)),        &
               grid%reg_mat(1)%eigvec_0(npt(1),npt(1)), &
               grid%reg_mat(1)%work(5*npt(1)))
     call Diagonalize(grid%reg_mat(1)%tr,               &
                 grid%reg_mat(1)%eigval_0,              &
                 grid%reg_mat(1)%eigvec_0,              &
                 grid%reg_mat(1)%work,                  &
                 npt(1))
     call Diagonalize(grid%reg_mat(1)%ham,              &
                 grid%reg_mat(1)%eigval,                &
                 grid%reg_mat(1)%eigvec,                &
                 grid%reg_mat(1)%work,                  &
                 npt(1))
     DEALLOCATE(grid%reg_mat(1)%work)
  END IF
  CALL ebct(grid%reg_mat(1)%ham,                        &
            grid%reg_mat(1)%eigvec,                     &
            grid%reg_mat(1)%eigvec,                     &
            npt(1),npt(1),npt(1))
  IF(prn(3)) THEN
     title='Test of Normalization'
     CALL prntrm(title,grid%reg_mat(1)%ham,             &
                 npt(1),npt(1),npt(1),npt(1),iout)
  END IF
1 Format(/,5x,'fourier weight = ',e15.8)
END SUBROUTINE Fourier_KE
!***********************************************************************
!***********************************************************************
!deck Hermite_KE.f
!***begin prologue     Hermite_KE
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           kinetic energy with Hermite weight function
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements
!***                   for a hermite weight function.  only one region.
!***
!***description
!***references
!***routines called
!***end prologue       Hermite_KE
  SUBROUTINE Hermite_KE(grid,reg_mat_hermite)
  IMPLICIT NONE
  TYPE(coordinates)                      :: grid  
  TYPE(hermite_matrices), DIMENSION(:)   :: reg_mat_hermite  
  TYPE(functions), DIMENSION(:)          :: reg_poly  
  INTEGER                                :: i 
  CHARACTER (LEN=80)                     :: title
  ALLOCATE( grid%reg_mat_hermite(1)%tr( npt(1), npt(1) ) )
  grid%reg_mat_hermite(1)%tr(:,:) = zero
  DO  i=1,npt(1)
      grid%reg_mat_hermite(1)%tr(i,:)                                         &
                                 =                                            &
      grid%reg_mat_hermite(1)%tr(i,:)                                         &
                                 +                                            &
      grid%reg_poly(1)%pr(i,i)   * grid%reg_pt_wt(1)%wtr(i)                   &
                                 *                                            &
      ( grid%reg_poly(1)%ddpr(i,:) - two * grid%reg_pt_wt(1)%qr(i)            &
                                   * grid%reg_poly(1)%dpr(i,:) )
  END DO
  DO i=1,npt(1)
     grid%reg_mat_hermite(1)%tr(i,i)                                          &
                               =                                              &
     grid%reg_mat_hermite(1)%tr(i,i)                                          &
                               +                                              &
     grid%reg_poly(1)%pr(i,i)  * grid%reg_pt_wt(1)%wtr(i)                     &
                               *                                              &
     ( grid%reg_pt_wt(1)%qr(i) * grid%reg_pt_wt(1)%qr(i) - one )              &
                               * grid%reg_poly(1)%pr(i,i)
  END DO
END SUBROUTINE Hermite_KE
  Call Pe_Fedvr(grid)
  grid%reg_mat_hermite(1)%tr(:,:) = dscale * grid%reg_mat_hermite(1)%tr(:,:)
  IF(prn(3)) THEN
     title='normalized kinetic energy matrix for hermite DVR basis'
     CALL prntrm(title,grid%reg_mat_hermite(1)%tr,npt(1),npt(1),npt(1),npt(1),iout)
  END IF
  ALLOCATE( grid%reg_mat_hermite(1)%ham( npt(1), npt(1) ) ) 
  grid%reg_mat_hermite(1)%ham(:,:) =  grid%reg_mat_hermite(1)%tr(:,:)
  DO i=1,npt(1)
     grid%reg_mat_hermite(1)%ham(i,i) =  grid%reg_mat_hermite(1)%ham(i,i) + grid%reg_vec(1)%vr(i)
  END DO
  IF(prn(3)) THEN
     title='normalized hamiltonian matrix for hermite DVR basis'
     CALL prntrm(title,grid%reg_mat_hermite(1)%ham,npt(1),npt(1),npt(1),npt(1),iout)
  END IF
  IF(.not.nodiag) THEN
     ALLOCATE( grid%reg_mat(1)%eigval(npt(1)),          &
               grid%reg_mat(1)%eigvec(npt(1),npt(1)),   &
               grid%reg_mat(1)%work(5*npt(1)))
     call Diagonalize(grid%reg_mat(1)%ham,              &
                 grid%reg_mat(1)%eigval,                &
                 grid%reg_mat(1)%eigvec,                &
                 grid%reg_mat(1)%work,                  &
                 npt(1))
     DEALLOCATE(grid%reg_mat(1)%work)
  END IF
END SUBROUTINE Hermite_KE
!***********************************************************************
!**********************************************************************
!deck LaGuerre_KE.f
!***begin prologue     LaGuerre_KE
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           kinetic energy with Laguerre weight function
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements
!***                   for a laguerre weight function. only one region.
!***
!***description
!***references
!***routines called
!***end prologue       LaGuerre_KE
  SUBROUTINE LaGuerre_KE(grid,reg_mat_laguerre)
  IMPLICIT NONE
  TYPE(coordinates)                       :: grid  
  TYPE(laguerre_matrices), DIMENSION(:)   :: reg_mat_laguerre  
  TYPE(functions), DIMENSION(:)           :: reg_poly  
  REAL(idp)                               :: scale
  INTEGER                                 :: i
  CHARACTER (LEN=80)                      :: title
  grid%reg_mat_laguerre(1)%tr(:,:) = zero
  DO  i=1,npt(1)
      grid%reg_mat_laguerre(1)%tr(i,:)                                &
                                =                                     &
       grid%reg_mat_laguerre(1)%tr(i,:)                               &
                                +                                     &
       grid%reg_poly(1)%pr(i,i) * grid%reg_pt_wt(1)%wtr(i)            &
                                *                                     &
       ( grid%reg_poly(1)%ddpr(i,:) - grid%reg_poly(1)%dpr(i,:) )
  END DO
  DO i=1,npt(1)
     grid%reg_mat_laguerre(1)%tr(i,i)                                 &
                                 =                                    &
     grid%reg_mat_laguerre(1)%tr(i,i)                                 &
                                 +                                    &
     grid%reg_poly(1)%pr(i,i)    * grid%reg_pt_wt(1)%wtr(i) * quarter &
                                 *                                    &
                           grid%reg_poly(1)%pr(i,i)
  END DO
  Call Pe_Fedvr(grid)
  grid%reg_mat_laguerre(1)%tr(:,:) = dscale * grid%reg_mat_laguerre(1)%tr(:,:)
  IF(prn(3)) THEN
     title='normalized kinetic energy matrix for laguerre DVR basis'
     CALL prntrm(title,grid%reg_mat_laguerre(1)%tr,npt(1),npt(1),npt(1),npt(1),iout)
  END IF
  ALLOCATE( grid%reg_mat_laguerre(1)%ham( npt(1), npt(1) ) ) 
  grid%reg_mat_laguerre(1)%ham(:,:) =  grid%reg_mat_laguerre(1)%tr(:,:)
  DO i=1,npt(1)
     grid%reg_mat_laguerre(1)%ham(i,i) =  grid%reg_mat_laguerre(1)%ham(i,i) + grid%reg_vec(1)%vr(i)
  END DO
  IF(prn(3)) THEN
     title='normalized hamiltonian matrix for laguerre DVR basis'
     CALL prntrm(title,grid%reg_mat_laguerre(1)%ham,npt(1),npt(1),npt(1),npt(1),iout)
  END IF
  IF(.not.nodiag) THEN
     ALLOCATE( grid%reg_mat(1)%eigval(npt(1)),          &
               grid%reg_mat(1)%eigvec(npt(1),npt(1)),   &
               grid%reg_mat(1)%work(5*npt(1)))
     call Diagonalize(grid%reg_mat(1)%ham,              &
                 grid%reg_mat(1)%eigval,                &
                 grid%reg_mat(1)%eigvec,                &
                 grid%reg_mat(1)%work,                  &
                 npt(1))
     DEALLOCATE(grid%reg_mat(1)%work)
  END IF
END SUBROUTINE LaGuerre_KE
!***********************************************************************
!***********************************************************************
!deck Pe_Fedvr.f
!***begin prologue     Pe_Fedvr
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the parts of the regional matrices
!***                   not depending on the m quantum number.  since a different
!***                   form for DVR basis functions are needed for even and odd
!***                   m, two different kinetic energy operators need to be constructed.
!***                   Later, the even and odd m kinetic energy will be built from these.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Pe_Fedvr
  SUBROUTINE Pe_Fedvr(grid)
  IMPLICIT NONE
  TYPE (coordinates)        :: grid
  INTEGER                   :: i
!
!
  ALLOCATE( grid%reg_vec(1:nreg) )
  DO i = 1, nreg
     ALLOCATE( grid%reg_vec(i)%vr(1:npt(i)))
     Call vrmat(grid%reg_vec(i)%vr,grid%reg_pt_wt(i)%qr,grid%label,npt(i),i)
  END DO
!
END SUBROUTINE Pe_Fedvr
!***********************************************************************
!***********************************************************************
!deck vrmat.f
!***begin prologue     vrmat
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate potential energy matrix elements
!***references
!***routines called
!***end prologue       vrmat
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE vrmat(v,pt,label,nr,region)
  IMPLICIT NONE
  INTEGER                                :: nr
  REAL(idp), DIMENSION(nr)               :: v
  REAL(idp), DIMENSION(nr)               :: pt
  CHARACTER (LEN=*)                      :: label
  CHARACTER (LEN=3)                      :: itoc
  CHARACTER (LEN=16)                     :: datkey
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: region
  INTEGER                                :: i
  REAL(idp)                              :: escal
  REAL(idp)                              :: charg
  REAL(idp)                              :: depth
  REAL(idp)                              :: len
  REAL(idp)                              :: shift
  REAL(idp)                              :: omega
  REAL(idp)                              :: awell
  REAL(idp)                              :: fac
  REAL(idp)                              :: fpkey
  REAL(idp)                              :: scale
  CHARACTER (LEN=80)                     :: typ
  CHARACTER (LEN=80)                     :: chrkey
  CHARACTER (LEN=80)                     :: atom
  LOGICAL                                :: dollar
  LOGICAL                                :: logkey
  LOGICAL                                :: axis
  LOGICAL                                :: toau
  LOGICAL                                :: useau
  LOGICAL                                :: totrap
  LOGICAL                                :: prnt
  REAL(idp), DIMENSION(2)                :: amp
  REAL(idp), DIMENSION(2)                :: expnt
  REAL(idp)                              :: e_c
  INTEGER                                :: number
  INTEGER                                :: nwell
  INTEGER                                :: n_p
  INTEGER                                :: intkey
  INTEGER                                :: n_scale
  INTEGER                                :: lendat
  INTEGER                                :: lenth
  v=zero
  if(reuse) then
     datkey = '$'//label//'_v_reg_'//itoc(1)
  else
     datkey = '$'//label//'_v_reg_'//itoc(region)
  END IF
  lendat=lenth(datkey)
  IF( dollar( datkey,card,atom,inp ) )  THEN
      typ=chrkey(card,'potential','none',' ')
      scale=one
      write(iout,1) typ
!
!     calculate the potential
!
      IF(typ == 'none') then
         call none(v,pt,nr,prnt)
      ELSE IF(typ == 'well') then
         depth=fpkey(card,'well_depth',zero,' ')
         call vwell(v,depth,nr,prnt)
      ELSE IF(typ == 'exponential') then
         amp(1)=fpkey(card,'amplitude',-1.d0,' ')
         expnt(1)=fpkey(card,'exponent',expnt,' ')
         call vexp(v,pt,amp,expnt,nr,prnt)
      ELSE IF(typ == 'yukawa') then
         amp(1)=fpkey(card,'amplitude',-1.d0,' ')
         expnt(1)=fpkey(card,'exponent',expnt,' ')
         call vyukawa(v,pt,amp,expnt,nr,prnt)
      ELSE IF(typ == 'power_exponential') then
         amp(1)=fpkey(card,'amplitude',1.d0,' ')
         expnt(1)=fpkey(card,'exponent',expnt,' ')
         n_p=intkey(card,'power',0,' ')
         call v_pow_exp(v,pt,amp,expnt,n_p,nr,prnt)
      ELSE IF(typ == 'sum_exponential') then
         call fparr(card,'amplitudes',amp,2,' ')
         call fparr(card,'exponents',expnt,2,' ')
         call vexp_sum(v,pt,amp,expnt,nr,prnt)
      ELSE IF(typ == 'coulomb') then
         charg=fpkey(card,'charge',-one,' ')
         call vcoul(v,pt,charg,nr,prnt)
      ELSE IF(typ == 'eberlonium') then
         charg=fpkey(card,'charge',-one,' ')         
         n_p=intkey(card,'power',0,' ')
         amp(1)=fpkey(card,'a',1.d0,' ')
         amp(2)=fpkey(card,'b',1.d0,' ')
         call v_eberlonium(v,pt,charg,amp(1),amp(2),n_p,nr,prnt)
      ELSE IF(typ == 'inverse_r4') then
         call vir4(v,pt,nr,prnt)
      ELSE IF(typ == 'rounded_well') then
         nwell=intkey(card,'n_well',ten,' ')
         awell=fpkey(card,'a_well',14.d0,' ')
         call vrwell(v,pt,awell,nwell,nr,prnt)
      ELSE IF(typ == 'harmonic_oscillator') then
!
!        take the Cesium atom as model.
!         
         atom=chrkey(card,'atom','generic',' ')
         IF(atom == 'cs') then
            mass=fpkey(card,'mass',2.2d-25,' ')
            omega=fpkey(card,'omega',ten,' ')
         ELSE IF(atom == 'na') then
            mass=3.8176d-26
            omega=13.846d0
         ELSE IF(atom == 'na_3d') then
            mass=3.8176d-26
            omega=177.0d0
            axis=logkey(card,'short_axis',.false.,' ')
            IF(axis) then
              omega=sqrt(2.d0)*omega
            endIF
         ELSE IF(atom == 'generic') then
           mass=fpkey(card,'mass',massau,' ')
           omega=fpkey(card,'omega',1.d0/(two*pi*timau),' ')
         endIF
         omega=omega*two*pi
         fac=mass*omega*omega*half
         IF(units == 'atomic_units') then
            write(iout,*) '   converting to atomic units'
            omega=omega*timau
            mass=mass/massau
            fac=mass*omega*omega*half
            hbar=one
         endIF
         write(iout,2) mass, omega
         call vhmo(v,pt,fac,nr,prnt)
      ELSE IF(typ == 'anharmonic_oscillator') then
         call vanhmo(v,pt,nr,prnt)
      ELSE IF(typ == 'expres') then
         call fparr(card,'amplitude',amp,2,' ')
         call fparr(card,'exponent',expnt,2,' ')
         shift=fpkey(card,'exponent_shift',zero,' ')
         call vres(v,pt,amp,expnt,shift,nr,prnt)
      ELSE IF(typ == 'periodic') then
         n_scale=intkey(card,'n_i',10,' ')         
         e_c=fpkey(card,'e_c',.001d0,' ')         
         awell=n_scale/e_c
         call vperiod(v,pt,awell,nr,prnt)
      ELSE IF(typ == 'spheroidal') then
         IF ( label == 'eta') THEN
              v(:) = - R_ab * ( z_b - z_a ) * pt(:)
         ELSE IF ( label == 'xi' ) THEN
              v(:) = - R_ab * ( z_a + z_b ) * pt(:)
         END IF
      ELSE
         call lnkerr('error in potential')
      endIF
  END IF
  dscale = -.5d0*hbar*hbar/mass
  IF(units == 'atomic_units') then
     dscale = -.5d0
  END IF     
  IF(prn(3)) then
     title='potential matrix elements for region = '//itoc(region)
     call prntrm(title,v,nr,1,nr,1,iout)
  END IF
 1    format(/,1x,'potential type = ',a32)
 2    format(/,1x,'oscillator mass      = ',e15.8, &
             /,1x,'oscillator-frequency = ',e15.8)
 3    format(/,1x,'length scale = ',e15.8,1x,'energy scale = ',e15.8)
  END SUBROUTINE vrmat
!***********************************************************************
!***********************************************************************
!deck Mat_Renormalization.f
!***begin prologue     Mat_Renormalization
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the parts of the regional matrices
!***                   not depending on the m quantum number.  since a different
!***                   form for DVR basis functions are needed for even and odd
!***                   m, two different kinetic energy operators need to be constructed.
!***                   Later, the even and odd m kinetic energy will be built from these.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Mat_Renormalization

  SUBROUTINE Mat_Renormalization(grid,reg_mat)
  IMPLICIT NONE
  TYPE (coordinates)                  :: grid
  TYPE(matrices),  DIMENSION(:)       :: reg_mat
  REAL(idp)                           :: dum
  INTEGER                             :: i
!
!

  IF ( nreg == 1 ) THEN
!
!      One Region only.  No change needed except normalization.
!
       i = 1
       Call Re_KE (grid%reg_mat(i)%tr,                                        &
                   dum,                                                       &
                   dum,                                                       &
                   grid%reg_pt_wt(i)%inv_sqrt_wtr,                            &
                   npt(i),                                                    &
                   i )
  ELSE
       i = 1
!
!      Region one is a special case.  Only the last right hand element needs to be changed
!      and the elements renormalized.
!
       Call Re_KE (grid%reg_mat(i)%tr,                                        &
                   dum,                                                       &
                   grid%reg_mat(i+1)%tr(1,1),                                 &
                   grid%reg_pt_wt(i)%inv_sqrt_wtr,                            &
                   npt(i),                                                    &
                   i )
!
!      Now do the general case.  Both the left and right elemenents are changed.
! 
       DO i = 2 , nreg - 1
          Call Re_KE (grid%reg_mat(i)%tr,                                     &
                      grid%reg_mat(i-1)%tr(npt(i-1),npt(i-1)),                &
                      grid%reg_mat(i+1)%tr(1,1),                              &
                      grid%reg_pt_wt(i)%inv_sqrt_wtr,                         &
                      npt(i),                                                 &
                      i )
       END DO
       i = nreg
!
!      Special case of last region.  Only the left hand elements needs to be changed.
!
          Call Re_KE (grid%reg_mat(i)%tr,                                     &
                      grid%reg_mat(i-1)%tr(npt(i-1),npt(i-1)),                &
                      dum,                                                    &
                      grid%reg_pt_wt(i)%inv_sqrt_wtr,                         &
                      npt(i),                                                 &
                      i )
  END IF
END SUBROUTINE Mat_Renormalization
!***********************************************************************
!***********************************************************************
!deck Odd_Mat_Renormalization.f
!***begin prologue     Odd_Mat_Renormalization
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the parts of the regional matrices
!***                   not depending on the m quantum number.  since a different
!***                   form for DVR basis functions are needed for even and odd
!***                   m, two different kinetic energy operators need to be constructed.
!***                   Later, the even and odd m kinetic energy will be built from these.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Odd_Mat_Renormalization

  SUBROUTINE Odd_Mat_Renormalization(grid,reg_mat_odd)
  IMPLICIT NONE
  TYPE(coordinates)                   :: grid  
  TYPE(odd_matrices),  DIMENSION(:)   :: reg_mat_odd
  REAL(idp)                           :: dum
  INTEGER                             :: i
!
!

  IF ( nreg == 1 ) THEN
!
!      One Region only.  No change needed except normalization.
!
       i = 1
       Call Re_KE (grid%reg_mat_odd(i)%tr,                                    &
                   dum,                                                       &
                   dum,                                                       &
                   grid%reg_pt_wt(i)%inv_sqrt_wtr,                            &
                   npt(i),                                                    &
                   i )
  ELSE
       i = 1
!
!      Region one is a special case.  Only the last right hand element needs to be changed
!      and the elements renormalized.
!
       Call Re_KE (grid%reg_mat_odd(i)%tr,                                    &
                   dum,                                                       &
                   grid%reg_mat_odd(i+1)%tr(1,1),                             &
                   grid%reg_pt_wt(i)%inv_sqrt_wtr,                            &
                   npt(i),                                                    &
                   i )
!
!      Now do the general case.  Both the left and right elemenents are changed.
! 
       DO i = 2 , nreg - 1
          Call Re_KE (grid%reg_mat_odd(i)%tr,                                 &
                      grid%reg_mat_odd(i-1)%tr(npt(i-1),npt(i-1)),            &
                      grid%reg_mat_odd(i+1)%tr(1,1),                          &
                      grid%reg_pt_wt(i)%inv_sqrt_wtr,                         &
                      npt(i),                                                 &
                      i )
       END DO
       i = nreg
!
!      Special case of last region.  Only the left hand elements needs to be changed.
!
          Call Re_KE (grid%reg_mat_odd(i)%tr,                                 &
                      grid%reg_mat_odd(i-1)%tr(npt(i-1),npt(i-1)),            &
                      dum,                                                    &
                      grid%reg_pt_wt(i)%inv_sqrt_wtr,                         &
                      npt(i),                                                 &
                      i )
  END IF
END SUBROUTINE Odd_Mat_Renormalization

!***********************************************************************
!***********************************************************************
!deck Diagonalize.f
!***begin prologue     Diagonalize
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            diagonalization of DVR/FEM hamiltonian.
!***
!***references
!***routines called
!***end prologue       Diagonalize
  SUBROUTINE Diagonalize(h,eigv,eigvec,work,nphy)
  IMPLICIT NONE
  INTEGER                                :: nphy
  INTEGER                                :: info
  INTEGER                                :: i
  REAL(idp), DIMENSION(:,:)              :: h
  REAL(idp), DIMENSION(:,:)              :: eigvec
  REAL(idp), DIMENSION(:)                :: eigv
  REAL(idp), DIMENSION(:)                :: work
  CHARACTER (LEN=80)                     :: title
  eigvec(:,:)  = h(:,:)
  IF(prn(9)) THEN
     title='physical hamiltonian'
     call prntrm(title,h,nphy,nphy,nphy,nphy,iout)
  END IF
  CALL dsyev('v','l',nphy,h,nphy,eigv,work,5*nphy,info)
  IF(prn(9)) THEN
     title='eigenvalues'
     CALL prntrm(title,eigv,nphy,1,nphy,1,iout)
  END IF
  IF(prn(10)) THEN
     title='eigenvectors'
     CALL prntrm(title,eigvec,nphy,nphy,nphy,nphy,iout)
  END IF
END SUBROUTINE Diagonalize
!***********************************************************************
!***********************************************************************
!deck H0_FEDVR_Matrices.f
!***begin prologue     H0_FEDVR_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE H0_FEDVR_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (hamiltonian), DIMENSION(:,:), ALLOCATABLE  :: reg_ham
  INTEGER                                          :: ir
  INTEGER                                          :: i
  INTEGER                                          :: m
!
  DO m = 0, m_max, 2
     ALLOCATE( grid%reg_ham(1:nreg,m) )
     DO ir = 1, nreg
        ALLOCATE( grid%reg_ham(ir,m)%ham(1:npt(ir),1:npt(ir)) )
        grid%reg_ham(ir,m)%ham(:,:) = zero
        DO i = 1, npt(ir)
           grid%reg_ham(ir,m)%ham(i,i) =                            &
                            m * m *grid%reg_pt_wt(ir)%inv_qr_fac(i)
        END DO 
        grid%reg_ham(ir,m)%ham(:,:) = grid%reg_ham(ir,m)%ham(:,:)   &
                                             +                      &
                                      grid%reg_mat(ir)%tr(:,:)
     END DO
  END DO
  DO m = 1, m_max, 2
     ALLOCATE( grid%reg_ham(1:nreg,m) )
     DO ir = 1, nreg
        ALLOCATE( grid%reg_ham(ir,m)%ham(1:npt(ir),1:npt(ir)) )
        grid%reg_ham(ir,m)%ham(:,:) = zero
        DO i = 1, npt(ir)
           grid%reg_ham(ir,m)%ham(i,i) =                            &
                              m * m *grid%reg_pt_wt(ir)%inv_qr_fac(i)
        END DO 
        grid%reg_ham(ir,m)%ham(:,:) = grid%reg_ham(ir,m)%ham(:,:)   &
                                             +                      &
                                      grid%reg_mat_odd(ir)%tr(:,:)
     END DO
  END DO
END SUBROUTINE H0_FEDVR_Matrices
!***********************************************************************
!***********************************************************************
           END MODULE FEDVR_Hamiltonian_Module
!***********************************************************************
!***********************************************************************
