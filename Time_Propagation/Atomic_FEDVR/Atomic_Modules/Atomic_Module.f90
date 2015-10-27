!***********************************************************************
                           MODULE Atomic_Module
                           USE dvr_shared
                           USE dvr_global
                           USE dvrprop_global
                           IMPLICIT NONE
             INTEGER                       :: n_tri 
             CHARACTER(LEN=80)             :: ints

             REAL*8                        :: fourpi=12.566370614359172953850572d0
             CHARACTER(LEN=80), EXTERNAL   :: chrkey
             CHARACTER(LEN=3), EXTERNAL    :: itoc
             LOGICAL, EXTERNAL             :: dollar, logkey
             INTEGER, EXTERNAL             :: intkey
             REAL*8, EXTERNAl              :: fpkey
             REAL*8, DIMENSION(:,:,:,:),                         &
                     ALLOCATABLE           :: D_LM
             INTEGER, DIMENSION(:,:),                            &
                      ALLOCATABLE          :: index
             REAL*8, DIMENSION(:),                               &
                     ALLOCATABLE           :: rho
             REAL*8, DIMENSION(:),                               &
                     ALLOCATABLE           :: v
             REAL*8, DIMENSION(:,:),                             &
                     ALLOCATABLE           :: t
             REAL*8, DIMENSION(:,:),                             &
                     ALLOCATABLE           :: work
             INTEGER,DIMENSION(:),                               &
                     ALLOCATABLE           :: ipvt
             REAL*8, DIMENSION(:,:),                             &
                     ALLOCATABLE           :: h_one
             REAL*8, DIMENSION(:,:),                             &
                     ALLOCATABLE           :: h_two
             INTEGER                       :: number_of_centers
             TYPE Orbitals
               INTEGER                     :: l_orb_max
               INTEGER                     :: L_max
             END TYPE Orbitals
             TYPE (Orbitals)               :: Orbital
             TYPE Atoms
             REAL*8, DIMENSION(:),                               &
                     ALLOCATABLE           :: atomic_charges             
             REAL*8, DIMENSION(:,:),                             &
                     ALLOCATABLE           :: atomic_locations
             END TYPE Atoms
             TYPE (Atoms)                  :: Atom_Properties             
             TYPE Angular_Momentum
               INTEGER                     :: l_orb_max
               INTEGER                     :: L_max
               INTEGER                     :: M
               REAL*8, DIMENSION(:,:),                           &
                       ALLOCATABLE         :: 3_J             
             END TYPE Angular_Momentum
             TYPE (Angular_Momentum)       :: Momentum
!***********************************************************************
                           INTERFACE Atomic_Input
                       END INTERFACE Atomic_Input
                           INTERFACE Atomic_Basis
                       END INTERFACE Atomic_Basis
                          INTERFACE One_Electron_Integrals
                       END INTERFACE One_Electron_Integrals
                           INTERFACE Two_Electron_Integrals
                       END INTERFACE Two_Electron_Integrals
                           INTERFACE Onemat
                       END INTERFACE Onemat
                           INTERFACE Twomat
                       END INTERFACE Twomat
                           INTERFACE Poisson
                       END INTERFACE Poisson
                           INTERFACE Rhofun
                       END INTERFACE Rhofun
!***********************************************************************
!***********************************************************************
                           Contains
!***********************************************************************
!***********************************************************************
!deck Atomic_Input
!***begin prologue     Atomic_Input
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input, atomic
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for computing dvr basis sets,
!***                   kinetic energy, one-electron and two electron
!***                   radial integrals.
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Atomic_Input
!***********************************************************************
  SUBROUTINE Atomic_Input
!***********************************************************************
  IMPLICIT NONE
  INTEGER                                     :: l
  INTEGER                                     :: l_1, l_2, m_1, m_2
  INTEGER                                     :: l_3, m_3
  INTEGER                                     :: i, j
  INTEGER                                     :: l_sum, l_diff, count
  INTEGER                                     :: l_11, l_22, l_33
  CHARACTER(LEN=3)                            :: itoc
  REAL*8                                      :: factor, dfac
!
!     read the input by stopping at the keyword in the input stream.
!
  IF ( dollar('$Atomic_Data',card,cpass,inp) ) THEN
!
!
!     this is an atomic calculation
!
      number_of_centers=intkey(card,'number_of_centers',1,' ')
      ALLOCATE(Atom_Properties%atomic_charges(number_of_centers),         &
               Atom_Properties%atomic_locations(3,number_of_centers))
      Momentum%l_orb_max =                                                & 
                         intkey(card,'maximum-orbital-angular-momentum',  &
                                0,' ') 
      Momentum%M = intkey(card,'total M ',0,' ') 
      Momentum%L_max = Momentum%l_orb_max + Momentum%l_orb_max 
      DO i=1,number_of_centers        
         Call fparr(card,'atomic_location_'//itoc(i),3,                 &
                    Atom_Properties%atomic_locations(:,i), ' ')
         WRITE(iout,1) i, Atom_Properties%atomic_locations(:,i),        &
                          Atom_Properties%atomic_charges(i)
      END DO
      write(iout,2) Momentum%l_orb_max, Momentum%L_max, Momentum%M 
      ints=chrkey(card,'integrals','all-integrals',' ')
      IF(ints == 'poisson') then
         dentyp=chrkey(card,'type-density','exponential',' ')
      END IF
      units='atomic-units'
      typke='dvr'
      diag=.true.
      proj=.false.
      genpts=logkey(card,'automate-points',.false.,' ')
      proj=.false.
      WRITE(iout,3) ints, dentyp
      coord(1)='r'
  END IF
  DO l_1=0, Momentum%l_orb_max
     DO l_2=0,l_1
        l_sum=l_1+l_2
        l_diff=iabs( l_1 - l_2)
        ALLOCATE(Momentum%3_J(l_1,l_2,l_diff:l_sum)
        DO l=l_diff, l_sum
           3_J(l_1,l_2,l) = Three_J_Coef(l_1,0,l_2,0,l,0)           
           3_J(l_2,l_1,l) =3_J(l_1,l_2,l)
        END DO
     END DO      
  END DO                 
  count=0
  DO l_1=0, Momentum%l_orb_max
     count = count + l_1 +l_1 + 1
  END DO
  ALLOCATE(index(count,2))
  count=0
  DO l_1=0, Momentum%l_orb_max
     DO m_1=-l_1,l_1
        count=count+1
        index(count,1)=l_1
        index(count,2)=m_1
     END DO
  END DO 
  factor = 1.d0/sqrt(fourpi)
  D_LM(:,:,:,:)=0.d0
  DO i=1,count
     l_1=index(i,1)
     m_1=index(i,2)
     l_11 = l_1 + l_1 + 1 
     DO j=1,i
        l_2=index(j,1)
        m_2=index(j,2)
        l_22 = l_2 + l_2 + 1 
        l_sum=l_1+l_2
        l_diff=iabs( l_1 - l_2)
        m_3 = m_2 - m_1
        DO l_3=l_diff,l_sum
           l_33 = l_3 + l_3 + 1   
           dfac = l_11 * l_22 * l_33           
           D_LM(i,j,l_3,m_3) = sqrt ( dfac ) * factor           &
                               *                                &
                           c_tmp(l_1,l_2,l_3)                   & 
                               *                                &
                           Three_J_Coef(l_1,m_1,l_2,m_2,l_3,m_3)     
        END DO                
     END DO
  END DO
  DEALLOCATE(index,c_tmp)
1 FORMAT(/,20x,'generating atomic dvr functions and operators',   &
         /,35x,'basic data',                                      &
         /,5x,'atom center                      = ',i3,           &  
         /,5x,'atom position(x,y,z)             = ',3(e15.8,1x),  &  
         /,5x,'charge                           = ',e15.8)
2 FORMAT(/,5x,'maximum orbital angular momentum = ',i3,           &
         /,5x,'maximum total angular momentum   = ',i4,           &
         /,5x,'total M                          = ',i4)
3 Format(/20x,'integrals                        = ',a16,          &
         /,5x,'type density                     = ',a16)
4 FORMAT(/,10x'Allowed M Values-zero is yes, one is no')
5 FORMAT(/,5x,'Shell = ',i3,(/,5x,'Values =',20i1))
!***********************************************************************
  END Subroutine Atomic_Input
!***********************************************************************
!***********************************************************************
!deck Atomic_Basis
!***begin prologue     Atomic_Basis
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input, atomic
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for computing dvr basis sets,
!***                   kinetic energy, one-electron and two electron
!***                   radial integrals.
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Atomic_Basis
!***********************************************************************
  SUBROUTINE Atomic_Basis
!***********************************************************************
  IMPLICIT NONE
!
  ALLOCATE(grid(1),num_reg(1),nfun_reg(maxreg,1))
!
! The input routine must be called retaining the basis function at the
! last point.
!
  CALL dvr_input(nphy(1),nglobal(1),coord(1))
!
  IF(diag) then
     ALLOCATE(grid(1)%eigv_0(nphy(1)),                       &
              grid(1)%eigvec_0(nphy(1),nphy(1)),             &
              grid(1)%eigv(nphy(1)),                         &
              grid(1)%eigvec(nphy(1),nphy(1)),               &
              grid(1)%srf_0(nphy(1),2),                      &
              grid(1)%srf(nphy(1),2))
  END IF
!
  ALLOCATE(grid(1)%pt(nphy(1)),                              &
           grid(1)%wt(nphy(1)),                              &
           grid(1)%f(nphy(1),nphy(1)),                       &
           grid(1)%df(nphy(1),nphy(1)),                      &
           grid(1)%ddf(nphy(1),nphy(1)),                     &
           grid(1)%ke(nphy(1),nphy(1)),                      &
           grid(1)%p_mom(nphy(1),nphy(1)),                   &
           grid(1)%h(nphy(1),nphy(1)),                       &
           grid(1)%v(nphy(1)),grid(1)%srf_prm(2))
!
!     We calculate a radial basis and the matrix elements
!     needed to set up all one and two electron matrix elements
!     The ke matrix contains the pure radial part of the kinetic
!     energy.  The same basis is used for all (lm) values.
!     Angular momentum issues will be dealt with later.
!
  CALL dvr_basis(pt_0(1),                                    &
                 pt_n(1),                                    &
                 grid(1)%pt,                                 &
                 grid(1)%wt,                                 &
                 grid(1)%f,                                  &
                 grid(1)%df,                                 &
                 grid(1)%ddf,                                &
                 grid(1)%ke,                                 &
                 grid(1)%p_mom,                              &
                 grid(1)%eigv_0,                             &
                 grid(1)%eigvec_0,                           &
                 grid(1)%h,                                  &
                 grid(1)%eigv,                               &
                 grid(1)%eigvec,                             &
                 grid(1)%v,                                  &
                 grid(1)%srf_prm,                            &
                 grid(1)%srf_0,                              &
                 grid(1)%srf,                                &
                 coord(1),nphy(1),nglobal(1))
  DEALLOCATE(grid(1)%df,grid(1)%ddf,grid(1)%p_mom,           &
             grid(1)%h,grid(1)%srf_prm)
  row(1) = 2* nphy(1) -2
  IF(diag) then
     DEALLOCATE(grid(1)%srf_0,grid(1)%srf)
  END IF
  num_reg(1) = nreg
  IF(bcl == 0) then
     npt(1) = npt(1) - 1
  END IF
  IF(bcr == 0 ) then
     npt(nreg) = npt(nreg) - 1
  END IF
  nfun_reg(:,1) = npt
!
!     Generate the bare coulomb interaction
!
  CALL vcoul(grid(1)%v,grid(1)%pt,atomic_charge,nphy(1),.false.)
!
!     Generate 1/(r*r)  This will allow us to set up the one electron
!     matrix elements for any angular momentum.
!
  ALLOCATE(grid(1)%ang_pot(nphy(1)))
  CALL vir2(grid(1)%ang_pot,grid(1)%pt,nphy(1),.false.)
  CALL vscale(grid(1)%ang_pot,grid(1)%ang_pot,.5d0,nphy(1))
!
  n_tri=nphy(1)*(nphy(1) + 1 )/2
!***********************************************************************
  END Subroutine Atomic_Basis
!***********************************************************************
!***********************************************************************
!deck One_Electron_Integrals
!***begin prologue     One_Elecgtron_Integrals
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input, atomic
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for computing dvr basis sets,
!***                   kinetic energy, one-electron and two electron
!***                   radial integrals.
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       One_Electron_Integrals
!***********************************************************************
  SUBROUTINE One_Electron_Integrals
  IMPLICIT NONE
!
  ALLOCATE(h_one(n_tri,0:l_orb_max))
!
!     Here we take the pieces and set up all the one electron matrix
!     elements.  Upon exit, h_one will contain the kinetic plus
!     one electon potential integrals for all angular momenta
!     specified.
!
  CALL Onemat(h_one,grid(1)%ke,grid(1)%v,grid(1)%ang_pot,nphy(1))
!***********************************************************************
  END Subroutine One_Electron_Integrals
!***********************************************************************
!***********************************************************************
!deck Three_J_Coef
!***begin prologue     Three_J_Coef
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           spherical harmonic
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            Compute 3-J coefficient
!***                   
!***description        See below
!***references
!***routines called    
!***end prologue       Three_J_Coef
!***********************************************************************
  FUNCTION Three_J_Coef(L1,M1,L2,M2,L,M)
  IMPLICIT NONE
  INTEGER                                 :: L1, L2, M1, M2, L ,M
  INTEGER, DIMENSION(6)                   :: idum
  REAL*8                                  :: Three_J_Fac, Three_J_Coef
  REAL*8, External                        :: Z_3j
  idum(1) = L1 + L1 + 1
  idum(2) = M1 + M1 + 1
  idum(3) = L2 + L2 + 1
  idum(4) = M2 + M2 + 1
  idum(5) = L + L + 1
  idum(6) = M + M + 1
  Three_J_Fac=(-1)**( ( idum(1) - idum(3) + idum(6)-1)/2 ) * sqrt( DBLE(idum(5)) )
  Three_J_Coef = Three_J_Fac *                                                    &
                 Z_3j( idum(1),idum(2),idum(3),idum(4),idum(5),- idum(6) + 2 )
!
!***********************************************************************
  END FUNCTION Three_J_Coef
!***********************************************************************
!---------------------------------------------------------------------- 
  FUNCTION Z_3j (j1,m1,j2,m2,j3,m3) 
  Implicit NONE 
  Real*8                              :: Z_3j 
  Integer, intent(in)                 :: j1,m1,j2,m2,j3,m3 
  Integer                             :: i,i_max,k,kk,m,iz,iz_min,iz_max 
  Real*8                              :: x,y,z 
  Integer, DIMENSION(3)               :: a,b 
  Integer, DIMENSION(16)              :: J 
!-------------------------------------------------------------------- 
! 
!     determines the value of the 3j-symbols without direct using of 
!     factorials. The following expression for the 3j-symbols is used: 
!         (A.P.JUCYS, A.A.BANDZAITIS, 1977) 
! 
!     3j{j1,m1,j2,m2,j3,m3} = delta(m1+m2,m3) * (2j3+1)^1/2 * {j1,j2,j3} * 
!       sqrt[ (j1+m1)!*(j1-m1)!*(j2+m2)!*(j2-m2)!*(j3+m3)!*(j3-m3)! ] 
!                         SUM(z) {   (-1)^z  / 
!          [ z! *  (j1+j2-j3-z)! * (j1-m1-z)! * (j2-m2-z)! * 
!                  (j3-j2+m1+z)! * (j3-j1-m2+z)! ] } 
! 
!     where {a,b,c}=sqrt[ (a+b-c)! * (a-b+c)! * (b-a+c)! / (a+b+c+1)! ] 
! 
!     If we introduce the auxiliary values a(i) and b(i) 
!     (see below the text of program) then 
! 
!     3j =         (-1) ^ Sum[a(i)] 
!          sqrt{ Pr[ (b(j)-a(i))! ] / [ Sum (b(j)-a(i))+1 ] } 
!                  Sum(z) { (-1)^z  / 
!          [  Pr[ (z-a(i))! ]  * Pr[ (b(j)-z)! ]   ] } 
! 
!     (below the moments are used in (2J+1)-representation) 
! 
!-------------------------------------------------------------------- 
  Z_3j=0.d0 
  IF(M1+M2+M3-3 /= 0) RETURN    ! check of conservation rules 
  J(1)= J1+J2-J3-1 
  J(2)= J1-J2+J3-1 
  J(3)= J2-J1+J3-1 
  J(4)= J1+M1-2 
  J(5)= J1-M1 
  J(6)= J2-M2 
  J(7)= J2+M2-2 
  J(8)= J3+M3-2 
  J(9)= J3-M3 
  Do I=1,9 
     IF(J(i) < 0.or.mod(J(i),2) == 1) THEN
        RETURN 
     END IF 
  End do 
  a(1) = 0                         ! auxiliary values 
  a(2) = (j2-j3-m1+1)/2 
  a(3) = (j1-j3+m2-1)/2 
  b(1) = (j1+j2-j3-1)/2 
  b(2) = (j1-m1)/2 
  b(3) = (j2+m2-2)/2 
  IZ_min=MAX0(a(1),a(2),a(3))      ! limits of the sum 
  IZ_max=MIN0(b(1),b(2),b(3)) 
  IF(IZ_max.LT.IZ_min) Return 
     Do I=1,3                         ! constant factorial parameters 
        Do K=1,3 
           J(I+3*K-3)=b(i)-a(k) 
        End do 
     End do 
     J(10)=(j1+j2+j3-3)/2+1 
     Do I=1,3 
        J(I+10)=IZ_min-a(i)               ! initial factorial parameters 
        J(I+13)=b(i)-IZ_min               ! in the sum 
     End do 
     Z=0.d0 
     DO IZ=IZ_min,IZ_max                 ! summation 
        I_max=0                            ! max. factorial 
        Do I=1,16 
           IF (J(i) > I_max) THEN
              I_max=J(i) 
           END IF
        End do 
        Y=1.d0 
        DO I=2,I_max         ! estimation of one term in sum 
           K=0                 ! K - the extent of the integer I in term 
           DO M=1,9 
              IF(J(M) >= I) THEN
                 K=K+1 
              END IF
           End do 
           IF(J(10) >= I) THEN
              K=K-1 
           END IF
           DO M=11,16 
              IF(J(M) >= I) THEN
                 K=K-2 
              END IF
           End do 
           IF(K == 0) Cycle 
              X=DBLE(I)                   ! Y = Y * I ** K/2 
              KK=IABS(K)/2 
              IF(KK > 0 ) THEN 
                 DO M=1,KK 
                    IF(K > 0) THEN
                       Y=Y*X 
                    END IF
                    IF(K < 0) THEN
                       Y=Y/X 
                    END IF
                 END DO 
              END IF 
              IF(mod(K,2) == 1) THEN
                 Y=Y*SQRT(X) 
              END IF
              IF(mod(K,2) == -1) THEN
                 Y=Y/SQRT(X) 
              END IF
        End do 
        IF(mod(IZ,2) == 1) THEN
           Y=-Y 
        END IF
        Z=Z+Y 
        Do I=11,13                  ! new factorial parameters in sum 
           J(I)=J(I)+1 
        End do 
        DO I=14,16 
           J(I)=J(I)-1 
        End do 
     End do                       ! end of summation 
     K=a(1)+a(2)+a(3) 
     IF(mod(k,2) /= 0) THEN
        Z=-Z 
     END IF
     Z_3j=Z 
!***********************************************************************
      END FUNCTION Z_3j 
!***********************************************************************
!deck Two_Electron_Integrals
!***begin prologue     Two_Electron_Integrals
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input, atomic
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for computing dvr basis sets,
!***                   kinetic energy, one-electron and two electron
!***                   radial integrals.
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Two_Electron_Integrals
!***********************************************************************
  SUBROUTINE Two_Electron_Integrals
!***********************************************************************
  IMPLICIT NONE
  IF( ints == 'all-integrals') THEN
      ALLOCATE( h_two(n_tri,0:l_max) )
      CALL Twomat(h_two,grid(1)%ke,grid(1)%ang_pot,grid(1)%pt,           &
                  grid(1)%wt,pt_n(1),nphy(1))
  ELSE IF( ints == 'poisson-equation') THEN
      CALL Poisson(grid(1)%ke,grid(1)%ang_pot,grid(1)%pt,grid(1)%wt,      &
                   grid(1)%f,pt_n(1),nphy(1))
  END IF
!***********************************************************************
  END Subroutine Two_Electron_Integrals
!***********************************************************************
!deck Onemat
!***begin prologue     Onemat
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input, atomic
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            Compute one electron mtrix elements.
!***description        
!***references
!***routines called    
!***end prologue       Onemat
!***********************************************************************
  SUBROUTINE Onemat(h_one,ke,v,ang_pot,n_one)
!***********************************************************************
  IMPLICIT NONE
  Real*8, DIMENSION(:,:)                   :: h_one
  Real*8, DIMENSION (:,:)                  :: ke
  Real*8, DIMENSION (:)                    :: v, ang_pot
  Integer                                  :: n_one, l, i, j
  Integer                                  :: l_fac  , i_diag, i_count
  h_one(:,:) = 0.d0
  DO l=0,l_orb_max
     l_fac = l*(l+1)
     i_diag=0
     DO i=1,n_one
        i_diag = i_diag + i 
        h_one(i_diag,l) = h_one(i_diag,l) + v(i) + l_fac * ang_pot(i)
     END DO
     i_count = 0
     DO i=1,n_one
        DO j=1,i
           i_count = i_count + 1
           h_one(i_count,l) = h_one(i_count,l) + ke(i,j)
        END DO
     END DO
  END DO
!***********************************************************************
  END Subroutine Onemat
!***********************************************************************
!***********************************************************************
!deck Twomat.f
!***begin prologue     Twomat
!***date written       030131   (yymmdd)
!***revision date               (yymmdd)
!***keywords           two electron radial integrals
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate two electron radial integrals
!***                   in a dvr representation by solving the
!***                   poisson equation for the inner integral
!***                   and then applying quadrature. 

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Twomat
!***********************************************************************
  SUBROUTINE Twomat(v_two,ke,ang_pot,pt,wt,pt_n,n_one)
!***********************************************************************
  IMPLICIT NONE
  INTEGER                                :: n_one
  REAL*8,   DIMENSION(:)                 :: pt, wt, ang_pot
  REAL*8,   DIMENSION(:,:)               :: ke
  REAL*8,   DIMENSION(:,:)               :: v_two
  REAL*8                                 :: scale, pt_n, ptlst, lscale 
  INTEGER                                :: l, lval, twol
  INTEGER                                :: i, j, info, count
!
  scale=-2.d0
  ALLOCATE(t(n_one-1,n_one-1),work(n_one,5),ipvt(n_one))
  DO l=0,l_max
!
!    put the ke matrix into t, add the angular part and then scale it
!    to correspond to the definitions given in the notes.
!
     t(:,:) = ke(1:n_one-1,1:n_one-1)
     lval=l*(l+1)
     DO i=1,n_one-1
        t(i,i) = t(i,i) + lval * ang_pot(i)
     END DO
     t = scale*t
!
!    Invert the (n-1) submatrix of t to solve the Poisson equation
!    with homogeneous boundary conditions.
!
     CALL dsytrf('l',n_one-1,t,n_one-1,ipvt,work,5*(n_one-1),info)
     CALL dsytri('l',n_one-1,t,n_one-1,ipvt,work,info)
!
!    calculate the integrals
!
     lscale = fourpi / (2*l+1)
     work(:,1)=pt**l
     work(:,2)=1.d0/(pt*sqrt(wt))
     twol=2*l+1
     ptlst=1.d0/(pt_n**twol)
     count=0
     DO i=1,n_one-1
        DO j=1,i
           count=count+1
           v_two(count,l) = -lval*t(i,j)*work(i,2)*work(j,2) + &
                             work(i,1)*work(j,1)*ptlst
        END DO
     END DO
!
!    Do the last basis function which does not use the t matrix.
!
     DO j=1,n_one
        count=count+1
        v_two(count,l) = work(n_one,1)*work(j,1)*ptlst
     END DO
     v_two(:,l) = lscale * v_two(:,l)
  END DO
  DEALLOCATE(t,work,ipvt)
!***********************************************************************
  END SUBROUTINE Twomat
!**********************************************************************
!deck Poisson.f
!***begin prologue     Poisson
!***date written       030131   (yymmdd)
!***revision date               (yymmdd)
!***keywords           Poisson Equation
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Solve Poisson equation using dvr functions

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Poisson
!**********************************************************************
  SUBROUTINE Poisson(ke,ang_pot,pt,wt,f,pt_n,n_one)
!**********************************************************************
  IMPLICIT NONE
  INTEGER                                :: n_one
  REAL*8,   DIMENSION(:)                 :: pt, wt, ang_pot
  REAL*8,   DIMENSION(:,:)               :: ke, f
  REAL*8                                 :: scale, pt_n, ptlst, a 
  INTEGER                                :: l, l_val, twol
  INTEGER                                :: i, info
  CHARACTER(LEN=80)                      :: title
!
  ALLOCATE(rho(n_one),v(n_one))
  CALL rhofun(rho,pt)
  scale=-2.d0
  ALLOCATE(t(n_one-1,n_one-1),work(n_one,5),ipvt(n_one))
  DO l=0,l_max
!
!    put the second derivative sub-matrix into t and then scale it
!
     t(:,:)=ke(1:n_one-1,1:n_one-1)
     l_val = l*(l+1)
     DO i=1,n_one-1
        t(i,i) = t(i,i) + l_val * ang_pot(i)
     END DO
     t = t * scale
!
!    solve linear equations
!
     CALL dsytrf('l',n_one-1,t,n_one-1,ipvt,work,5*(n_one-1),info)
     CALL dsytrs('l',n_one-1,1,t,n_one-1,ipvt,rho,1,info)
!
!    calculate the potential
!
     twol=2*l+1
     work(:,1)=pt**(l+1)
     ptlst=1.d0/(pt_n**twol)
     a=0.d0
     DO i=1,n_one
        a = a + rho(i)*wt(i)*work(i,1)
     END DO
     a=a/ptlst
     v(n_one)= a*work(n_one,1)
     DO i=1,n_one-1
        v(i) = rho(i)*f(i,i) + a*work(i,1)
     END DO
     v = v / pt 
     title='potential'
     call prntrm(title,v,n_one,1,n_one,1,iout)
  END DO
  DEALLOCATE(t,work,ipvt)
  DEALLOCATE(rho,v)
!***********************************************************************
  END SUBROUTINE Poisson
!**********************************************************************
!deck Rhofun.f
!***begin prologue     Rhofun
!***date written       030131   (yymmdd)
!***revision date               (yymmdd)
!***keywords           Density
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Electron Density
!***references

!***routines called    iosys, util and mdutil
!***end prologue       Rhofun
!**********************************************************************
  SUBROUTINE Rhofun(rho,pt)
!**********************************************************************
  IMPLICIT NONE
  REAL*8,   DIMENSION(:)                 :: pt
  REAL*8,   DIMENSION(:)                 :: rho
! 
  rho = exp(-pt)
!***********************************************************************
  END SUBROUTINE Rhofun
!**********************************************************************
END MODULE Atomic_Module
