!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample calls to init_free.f90 :
! alpha=0 grid
      call init_free(0.d0,r0,weights0,V_ext0,Ekin0,sqrtw0,Pmax0,Rmax0)  
! alpha=1 grid
      call init_free(1.d0,r1,weights1,V_ext1,Ekin1,sqrtw1,Pmax1,Rmax1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Subroutine init_free(alpha,r,weights,V_ext,Ekin,sqrtw,Pmax,Rmax)

! Set up position grid on the DVR points for cylindrical box
! construct kinetic energy matrix 

   use BdGmod
   implicit none

   real*8,  INTENT(IN)                  :: alpha
   real*8,  INTENT(IN OUT)              :: r(ntot),sqrtw(ntot)
   real*8,  INTENT(IN OUT)              :: weights(ntot)    
   real*8,  INTENT(IN OUT)              :: V_ext(ntot) 
   real*8,  INTENT(IN OUT)              :: Ekin(ntot,ntot) 
   real*8 , dimension(ntot+1,0:alpha)   :: zeros
   real*8 , dimension(1:ntot)           :: zerolist 
   real*8,  INTENT(IN OUT)              :: Pmax,Rmax
   real*8 , dimension(1:ntot)           :: p,x
   real*8 , dimension(1:ntot)           :: norm                            
   real*8 , dimension(1:ntot,1:ntot)    :: B,InvB
   real*8 , dimension(1:ntot)           :: bess
real*8 :: maxoffdia, maxdia

! Find N+1 first roots of m'th order besselfunction
   call besselzero(ntot+1,int(alpha),zeros,.false.)

   do i=1,ntot
      zerolist(i)=zeros(i,alpha)
   end do

! Set up radial position and momentum grids
! ntot grid points, as the last point is fixed by boundary condition
   do i=1,ntot
      r(i)=zerolist(i)
      p(i)=zerolist(i)
   end do

   Pmax=zeros(ntot+1,alpha)/Rmax

   r(:)=r(:)/Pmax
   p(:)=p(:)/Rmax

! Quadrature weights
! Norm is J_{m+1}(zeros) 
   call bessel_m(zerolist(:),norm,ntot,int(alpha+1.0d0),.false.)

   sqrtw(:)=sqrt(2.0d0)/norm(:)/Pmax
   weights(:)=sqrtw(:)**2

! Symmetric transformation matrix
   B(:,:)=0.0d0
   do i=1,ntot
      x(:)=p(i)*r(:)
      call bessel_m(x(:),bess,ntot,int(alpha),.false.)
      do j=1,ntot
         B(i,j)=2.0d0*bess(j)/(norm(i)*norm(j)*Rmax*Pmax)
     end do 
   end do

   InvB=transpose(B)

! External potential
   V_ext(:)=0.0d0

!Set up kinetic energy matrix elements
   Ekin(:,:)=0.0d0
   do i=1,ntot
      Ekin(i,i)=0.5d0*p(i)**2
   end do
   Ekin(:,:)=matmul(Ekin,InvB)
   Ekin(:,:)=matmul(B,Ekin)

! Correct for the centrifugal potential in diagonalize.f90
! diagonalize.f90 sets up Hamiltonian which has explicit 
! centrifugal term.
   do i=1,ntot
      Ekin(i,i)=Ekin(i,i)-(alpha**2)/(2.0d0*r(i)**2)
   end do

end subroutine init_free


! ----------------------------------------------------------------------
   subroutine interpolate
! ----------------------------------------------------------------------
! Sets up interpolation matrices to go from 'odd' grid to 'even' grid
! Depending on the angular momentum value either u or v is transformed

   use BdGmod
   implicit none

   real*8 , dimension(1:ntot)        :: x
   real*8 , dimension(1:ntot)        :: bess
   real*8 , dimension(1:ntot,1:ntot) :: bessa,bessb

! Transformation matrix from 'odd' to 'even' grid
      LIP10(:,:)=0.0d0 
      do q=1,ntot
         x(:)=Pmax1*r1(q)*r0(:)/Rmax1
         call bessel_m(x(:),bess,ntot,1,.false.)
         bessa(:,q)=bess(:)
         x(:)=Pmax1*r1(q)*r1(:)/Rmax1
         call bessel_m(x(:),bess,ntot,1,.false.)
         bessb(:,q)=bess(:)
      end do   
      do i=1,ntot
         do j=1,ntot
            do q=1,ntot     
               LIP10(i,j)=LIP10(i,j)+(Pmax1/Rmax1)**2*weights1(q)*bessa(i,q)*bessb(j,q)
            end do
            LIP10(i,j)=LIP10(i,j)*weights1(j)
         end do 
      end do

! Transformation matrix from 'even' to 'odd' grid
      LIP01(:,:)=0.0d0
      do q=1,ntot
         x(:)=Pmax0*r0(q)*r1(:)/Rmax0
         call bessel_m(x(:),bess,ntot,0,.false.)
         bessa(:,q)=bess(:)
         x(:)=Pmax0*r0(q)*r0(:)/Rmax0
         call bessel_m(x(:),bess,ntot,0,.false.)
         bessb(:,q)=bess(:)
      end do  
      do i=1,ntot
         do j=1,ntot
            do q=1,ntot
               LIP01(i,j)=LIP01(i,j)+(Pmax0/Rmax0)**2*weights0(q)*bessa(i,q)*bessb(j,q)
            end do 
            LIP01(i,j)=LIP01(i,j)*weights0(j)
         end do 
      end do
    
! Transformation matrices
! Syntax: transform_u10 transforms u from the 'odd' grid to the 'even'
! grid, leaving v unchanged (expressed on the 'even' grid)
! Transformations like these are necessary when vortex=.true.
! as u/v have even/odd or odd/even angular momentum and hence are
! expressed on different grids in the Hamiltonian
! When vortex=.false. u and v are on the same grid and we must use
! transfomation matrices transform_01 and transform_10
   transform_u10(:,:)=0.0d0
   transform_v10(:,:)=0.0d0
   transform_u01(:,:)=0.0d0
   transform_v01(:,:)=0.0d0
   transform_01(:,:)=0.d0 
   transform_10(:,:)=0.d0
! 1st quadrant
   do i=1,ntot
      transform_v01(i,i)=1.0d0
      transform_v10(i,i)=1.0d0
      do j=1,ntot
         transform_u10(i,j)=LIP10(i,j)/sqrtw1(j)
         transform_u01(i,j)=LIP01(i,j)/sqrtw0(j)
         transform_10(i,j)=LIP10(i,j)/sqrtw1(j) 
         transform_01(i,j)=LIP01(i,j)/sqrtw0(j) 
      end do
   end do
! 2nd and 3rd quadrant are zero
! 4th quadrant
   do i=ntot+1,2*ntot 
      transform_u10(i,i)=1.0d0 
      transform_u01(i,i)=1.0d0
      do j=ntot+1,2*ntot
         transform_v10(i,j)=LIP10(i-ntot,j-ntot)/sqrtw1(j-ntot)
         transform_v01(i,j)=LIP01(i-ntot,j-ntot)/sqrtw0(j-ntot)
         transform_10(i,j)=LIP10(i-ntot,j-ntot)/sqrtw1(j-ntot)
         transform_01(i,j)=LIP01(i-ntot,j-ntot)/sqrtw0(j-ntot) 
      end do
   end do  
    
   end subroutine interpolate


! ----------------------------------------------------------------------
   subroutine diagonalize(alpha,vort,d0,d1,h0,h1,E_Fermi)
! ----------------------------------------------------------------------

!  Diagonalizes the Bogoliubov-de Gennes equations        
!  returns u, v and eigenvalues epsilon

   use BdGmod
   implicit none

   real*8                          :: alpha,beta
   logical                         :: vort
   real*8                          :: cent_m,cent_p
   real*8 , dimension(1:ntot)      :: x0,x1
   real*8 , dimension(1:ntot)      :: d0
   real*8 , dimension(1:ntot)      :: d1
   real*8 , dimension(1:ntot)      :: h0
   real*8 , dimension(1:ntot)      :: h1
   real*8                          :: E_Fermi
   integer*8 , dimension(1:2*ntot) :: sign 

! NB mu is set to be integer value  (m=mu-0.5)!!
! This is only important for vortex solutions
! where the ang mom of u and v differ by a factor of 1
! ie mu in 4th quadrant must be changed to m+1 !!!

! Use angular momentum zero grid for even angular momentum states
      if(mod(alpha,2.d0).eq.0) then
         grid='even'
! Use angular momentum one grid for odd angular momentum states 
      else
         grid='odd'
      end if

   sign(:)=0
 
   if(vort) then
      beta=alpha+1.d0
   else
      beta=alpha
   end if

   if(potential == 'harmonic') then
      x0(:)=r0(:)
      x1(:)=1.d0/r1(:)
   elseif(potential == 'free') then
      x0(:)=1.d0
      x1(:)=1.d0
   end if   

! Coefficient of centrifugal potential
   if(potential == 'free') then
      cent_m=alpha**2
      cent_p=beta**2
   elseif(potential == 'harmonic') then
      cent_m=(alpha**2-0.25d0)
      cent_p=(beta**2-0.25d0)
   end if

! set up Hamiltonian in Boguliubov 2 component formalism
   ham(:,:)=0.0d0
! 1st quadrant
   do i=1,ntot
      if(grid == 'even') then
         ham(i,i)=V_ext0(i)+h0(i)+cent_m/(2.d0*r0(i)**2)-E_Fermi+0.5d0*kz**2
      elseif(grid == 'odd') then
         ham(i,i)=V_ext1(i)+h1(i)+cent_m/(2.d0*r1(i)**2)-E_Fermi+0.5d0*kz**2
      end if
      do j=1,ntot
         if(grid == 'even') then
            ham(i,j)=ham(i,j)+Ekin0(i,j)
         elseif(grid == 'odd') then
            ham(i,j)=ham(i,j)+Ekin1(i,j)
         end if
      end do
   end do
! Put in interaction terms
! 2nd quadrant
   if(grid == 'even') then   
      do i=1,ntot 
         if(vort) then
            do j=ntot+1,2*ntot
               ham(i,j)=d0(i)*(sqrtw0(i)/sqrtw1(j-ntot))*x0(i)*LIP10(i,j-ntot) 
           end do
         else
            ham(i,i+ntot)=d0(i)
         end if
      end do
   elseif(grid == 'odd') then
      do i=1,ntot 
         if(vort) then
            do j=ntot+1,2*ntot
               ham(i,j)=d1(i)*(sqrtw1(i)/sqrtw0(j-ntot))*LIP01(i,j-ntot)*x1(i)
            end do
         else
            ham(i,i+ntot)=d1(i)
         end if 
      end do
   end if   
! 3rd quadrant
   if(grid == 'even') then
      do i=ntot+1,2*ntot
         if(vort) then 
            do j=1,ntot
               ham(i,j)=d1(i-ntot)*(sqrtw1(i-ntot)/sqrtw0(j))*LIP01(i-ntot,j)*x1(i-ntot)
            end do
         else
            ham(i,i-ntot)=d0(i-ntot)
         end if
      end do
   elseif(grid == 'odd') then
      do i=ntot+1,2*ntot 
         if(vort) then
            do j=1,ntot
               ham(i,j)=d0(i-ntot)*(sqrtw0(i-ntot)/sqrtw1(j))*x0(i-ntot)*LIP10(i-ntot,j)
            end do
         else
            ham(i,i-ntot)=d1(i-ntot)
         end if 
      end do
   end if
! 4th quadrant
   if(vort) then
      do i=ntot+1,2*ntot
         if(grid == 'odd') then
            ham(i,i)=V_ext0(i-ntot)+h0(i-ntot)+cent_p/(2.d0*r0(i-ntot)**2)-E_Fermi+0.5d0*kz**2
         elseif(grid == 'even') then 
            ham(i,i)=V_ext1(i-ntot)+h1(i-ntot)+cent_p/(2.d0*r1(i-ntot)**2)-E_Fermi+0.5d0*kz**2
         end if
         do j=ntot+1,2*ntot
            if(grid == 'odd') then      
               ham(i,j)=-ham(i,j)-Ekin0(i-ntot,j-ntot)
            elseif(grid == 'even') then
               ham(i,j)=-ham(i,j)-Ekin1(i-ntot,j-ntot)
            end if  
         end do
      end do 
   else
      do i=ntot+1,2*ntot
         if(grid == 'even') then
            ham(i,i)=V_ext0(i-ntot)+h0(i-ntot)+cent_m/(2.d0*r0(i-ntot)**2)-E_Fermi+0.5d0*kz**2
         elseif(grid == 'odd') then 
            ham(i,i)=V_ext1(i-ntot)+h1(i-ntot)+cent_m/(2.d0*r1(i-ntot)**2)-E_Fermi+0.5d0*kz**2
         end if
         do j=ntot+1,2*ntot
            if(grid == 'even') then      
               ham(i,j)=-ham(i,j)-Ekin0(i-ntot,j-ntot)
            elseif(grid == 'odd') then
               ham(i,j)=-ham(i,j)-Ekin1(i-ntot,j-ntot)
            end if  
         end do
      end do 
   end if 

! diagonalize to find eigenstates and eigenvectors
   call DSYEV('V','U',2*ntot,ham,LDA,W,WORK,LWORK,INFO)

! Keep positive energy solutions
   do i=1,2*ntot
      if(symmetry == 'cylindrical') then
         if(W(i).lt.0.d0) then 
            if(vort) then
               sign(i)=2
            end if
! Negative ang mom solutions with positive energy are given
! by positive ang mom solutions with negative energy  
          Energy(i)=-W(i)         
            do j=1,ntot
! E -> -E, u -> v*, v -> -u*
! u and v are real so no complex conjugation is needed
               gamma(j,i)=ham(j+ntot,i)        
               gamma(j+ntot,i)=-ham(j,i)
            end do 
         else
            if(vort) then
               sign(i)=1
            end if  
            Energy(i)=W(i)
            gamma(:,i)=ham(:,i)      
         end if
! NB in the case of no vortex for m=0 only positive eigenvalues contribute
! thus for m=0, Energy contains both positive and negative values,
! and later on we selctively use only Energy>=0
         if((.not.vort).and.(alpha.eq.0.d0)) then
            Energy(i)=W(i)
            gamma(:,i)=ham(:,i)     
         end if 
      elseif(symmetry == 'spherical') then
! NB negative eigenvalues are retained for the spherical case
         Energy(i)=W(i)
         gamma(:,i)=ham(:,i) 
      end if 
   end do

! Find quasi-particle states
   call find_u_v(vort,sign)

   end subroutine diagonalize

! ----------------------------------------------------------------------
   subroutine find_u_v(vort2,sign2)
! ----------------------------------------------------------------------

! Finds the quasi-particle functions u, v from the radial function
! solutions to the BdG equations
! Note that u,v has 2*ntot columns, as it contains the lowest states for m
! and for -m, all with E>0 (If E_F=0 there will be ntot of each)

   use BdGmod
   implicit none

   logical                         :: vort2
   integer*8 , dimension(1:2*ntot) :: sign2

! Interpolate 'odd' grid eigenvectors unto 'even' grid
! and vice-versa
! NB if vortex=.true. then u and v are expressed on different grids
   if(vort2) then
      do j=1,2*ntot
         if(((grid == 'even').and.(sign2(j).eq.1)).or.((grid == 'odd').and.(sign2(j).eq.2))) then
            gamma0(:,j)=matmul(transform_v10,gamma(:,j))
            gamma1(:,j)=matmul(transform_u01,gamma(:,j))
         elseif(((grid == 'odd').and.(sign2(j).eq.1)).or.((grid == 'even').and.(sign2(j).eq.2))) then
            gamma0(:,j)=matmul(transform_u10,gamma(:,j))
            gamma1(:,j)=matmul(transform_v01,gamma(:,j))
         end if
      end do
! If vortex=.false. then u and v are expressed on the same grid
   else
      if(grid == 'even') then
         gamma0(:,:)=gamma
         gamma1(:,:)=matmul(transform_01,gamma)
      elseif(grid == 'odd') then
         gamma0(:,:)=matmul(transform_10,gamma)
         gamma1(:,:)=gamma
      end if 
   end if


      do i=1,ntot
         u0(i,:)=gamma0(i,:)
         v0(i,:)=gamma0(i+ntot,:)
         u1(i,:)=gamma1(i,:)
         v1(i,:)=gamma1(i+ntot,:)
      end do
! Complete interpolation to get correct u,v 
      do j=1,2*ntot  
         if(grid == 'even') then
            if(vort2) then
               if(sign2(j).eq.1) then
                  u0(:,j)=u0(:,j)/sqrtw0(:)
                  v0(:,j)=v0(:,j)
                  u1(:,j)=u1(:,j)
                  v1(:,j)=v1(:,j)/sqrtw1(:)
                elseif(sign2(j).eq.2) then
                  u0(:,j)=u0(:,j)
                  v0(:,j)=v0(:,j)/sqrtw0(:)
                  u1(:,j)=u1(:,j)/sqrtw1(:)
                  v1(:,j)=v1(:,j)
                end if 
            else
               u0(:,j)=u0(:,j)/sqrtw0(:)
               v0(:,j)=v0(:,j)/sqrtw0(:)
               u1(:,j)=u1(:,j)
               v1(:,j)=v1(:,j)
            end if
         elseif(grid == 'odd') then
            if(vort2) then
               if(sign2(j).eq.1) then 
                  u0(:,j)=u0(:,j)
                  v0(:,j)=v0(:,j)/sqrtw0(:)
                  u1(:,j)=u1(:,j)/sqrtw1(:)
                  v1(:,j)=v1(:,j)
               elseif(sign2(j).eq.2) then 
                  u0(:,j)=u0(:,j)/sqrtw0(:)
                  v0(:,j)=v0(:,j)
                  u1(:,j)=u1(:,j)
                  v1(:,j)=v1(:,j)/sqrtw1(:)                      
               end if 
            else
               u0(:,j)=u0(:,j)
               v0(:,j)=v0(:,j)
               u1(:,j)=u1(:,j)/sqrtw1(:)
               v1(:,j)=v1(:,j)/sqrtw1(:)
            end if
         end if 
      end do      

      call normalize


   end subroutine find_u_v





!deck bessel_m

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-20  Time: 15:30:07

!***begin prologue     bessel_m
!***date written       010802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           generate bessel function of order m
!***author             Nygaard, Nicolai (NIST)
!***source             math
!***purpose            Find order m Bessel function
!***                   
!***                   NB m must be integer
!***

!***references         is adaption of subroutine bessel

!***routines called
!***end prologue       nwtrap

SUBROUTINE bessel_m(x,bess_m,np,m,prnt)

IMPLICIT INTEGER (a-z)
REAL*8, INTENT(IN)                       :: x(np)
REAL*8                                   :: j(np,0:m)
REAL*8, INTENT(OUT)                      :: bess_m(np) 
REAL*8                                   :: dj(np,0:m)
REAL*8                                   :: factrl
INTEGER, INTENT(IN)                      :: np
INTEGER, INTENT(IN)                      :: m
LOGICAL, INTENT(IN)                      :: prnt
INTEGER, PARAMETER :: dim=1000
COMMON /io/ inp, iout
REAL*8 zero, one, two, epsilon
REAL*8 half, pi, twonpi, fournpi, xinv
REAL*8 norm, large, renorm
REAL*8 besj0, besj1, j0, j1
REAL*8 jnpl1, jn, jnm1

CHARACTER (LEN=80) :: title


DATA zero, one, two / 0.d+00, 1.d+00, 2.d+00 /
DATA half, large, renorm / .5D0, 1.d+250, 1.d-250  /
DATA pi / 3.141592653589793238462643D+00 /

twonpi=two/pi
fournpi=two*twonpi

!----------------------------------------------------------------------c
!            upward recursion, downward recursion or series            c
!----------------------------------------------------------------------c
DO  i=1,np
  xinv=1.d0/x(i)
  j0 = besj0(x(i))
  j1 = besj1(x(i))
  j(i,0) = j0
  IF(m > 0) THEN
    j(i,1) = j1
  END IF
  IF(m > 1) THEN
    
!           test for up or down recursion
    
    nn=x(i)
    IF(nn > m-1) THEN
      
!                  recur up
      
      j(i,0)=j0
      j(i,1)=j1
      nbeg=1
      DO  n=2,m
        j(i,n)=(nbeg+nbeg)*j(i,n-1)*xinv - j(i,n-2)
        nbeg=nbeg+1
      END DO
      
!              series or downward recursion
      
    ELSE
      IF(ABS(x(i)) > one) THEN
!----------------------------------------------------------------------c
!              find the value of upper needed to accurately get        c
!              the needed n values                                     c
!----------------------------------------------------------------------c
        upper=m
        upper=msta1(x(i),200)
        IF(upper < m) THEN
          upper=m + m
        ELSE
          upper=msta2(x(i),m,15)
        END IF
        upper=MAX(upper,m+2)
        jnpl1=zero
        jn=renorm
        onelss=upper-1
        nbeg=onelss
        DO  n=onelss-1,0,-1
          jnm1 = ( nbeg + nbeg )*jn*xinv -  &
              jnpl1
          IF(ABS(jnm1) > large) THEN
            jnm1=jnm1*renorm
            jn=jn*renorm
          END IF
          IF(n <= m) THEN
            j(i,n) = jnm1
          END IF
          jnpl1=jn
          jn=jnm1
          nbeg = nbeg -1
        END DO
!----------------------------------------------------------------------c
!                 normalize the j                                      c
!----------------------------------------------------------------------c
        IF(ABS(j0) > ABS(j1)) THEN
          norm=j0/j(i,0)
        ELSE
          norm=j1/j(i,1)
        END IF
        DO  n=0,m
          j(i,n) = norm*j(i,n)
        END DO
      ELSE
        DO  n=2,m
          CALL bser(x(i),j(i,n),dj(i,n),n,0)
        END DO
      END IF
    END IF

  END IF
END DO

bess_m(:)=j(:,m)

IF (prnt) THEN
  title='regular bessel functions'
  CALL prntfm(title,j,np,m+1,np,m+1,iout)
END IF
RETURN
END SUBROUTINE bessel_m



!deck bessely_m

!***begin prologue     bessely_m
!***date written       010813   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           generate irregular bessel function of order m
!***author             Nygaard, Nicolai (NIST)
!***source             math
!***purpose            Find order m irregular Bessel function
!***                   
!***                   NB m must be integer
!***

!***references         is adaption of subroutine bessel

!***routines called
!***end prologue       nwtrap

SUBROUTINE bessely_m(x,bessy_m,np,m,prnt)

IMPLICIT INTEGER (a-z)
REAL*8, INTENT(IN)                       :: x(np)
REAL*8                                   :: y(np,0:m)
REAL*8, INTENT(OUT)                      :: bessy_m(np) 
REAL*8                                   :: dy(np,0:m)
REAL*8                                   :: factrl
INTEGER, INTENT(IN)                      :: np
INTEGER, INTENT(IN)                      :: m
LOGICAL, INTENT(IN)                      :: prnt
INTEGER, PARAMETER :: dim=1000
COMMON /io/ inp, iout
REAL*8 zero, one, two, epsilon
REAL*8 half, pi, twonpi, fournpi, xinv
REAL*8 norm, large, renorm
REAL*8 besy0, besy1, y0, y1

CHARACTER (LEN=80) :: title



DATA zero, one, two / 0.d+00, 1.d+00, 2.d+00 /
DATA half, large, renorm / .5D0, 1.d+250, 1.d-250  /
DATA pi / 3.141592653589793238462643D+00 /

twonpi=two/pi
fournpi=two*twonpi

DO  i=1,np
  xinv=1.d0/x(i)
  y0 = besy0(x(i))
  y1 = besy1(x(i))
  y(i,0) = y0
  IF(m > 0) THEN
    y(i,1) = y1
  END IF
 
!----------------------------------------------------------------------c
!           recur upward for y(i,n)                                    c
!----------------------------------------------------------------------c
! (This is always stable)

    nbeg=1
    DO  n=2,m
      y(i,n) = ( nbeg + nbeg )*y(i,n-1)*xinv - y(i,n-2)
      IF(ABS(y(i,n)) > large) THEN
        y(i,n)=large
      END IF
      nbeg=nbeg+1
    END DO

END DO

bessy_m(:)=y(:,m)

IF (prnt) THEN
  title='irregular bessel functions'
  CALL prntfm(title,y,np,m+1,np,m+1,iout)
END IF
RETURN
END SUBROUTINE bessely_m


subroutine besselzero(mzro,nmax,besszro,writeout)

!***begin prologue     besselzero
!***date written       010622   (yymmdd)
!***revision date      010622   (yymmdd)
!***author             Nygaard, Nicolai (NIST)
!***purpose            generate zeros of integer order bessel functions
!***description        uses bisection and Newton-Raphson to locate zeros
!***references

!***routines called    bessel, nwtrap, intrv
!***end prologue       

implicit none

real*8  :: right,left,convg,del
integer :: mzro,step,niter,nmax,ntot,cntzro,locz,lroot,number
integer :: i,n,nroot  
real*8 , dimension(:) , ALLOCATABLE            :: x,z,roots
real*8 , dimension(:,:) , ALLOCATABLE          :: j,dj,y,dy
real*8 , dimension(:), ALLOCATABLE             :: bess_n
real*8 , dimension(mzro,0:nmax) , INTENT(OUT)  :: besszro
logical                                        :: writeout

! Parameter values
right=(mzro+0.5d0*nmax)*4.0d0 ! Right boundary (see A&S 9.5.12)
left=1.d-06                   ! Left boundary
!mzro=50                      ! Number of zeros to be found
step=4000                     ! Number of steps
niter=100                     ! Maximum number of iterations
convg=1.d-08                  ! Convergence criterion
del=(right-left)/step         ! spacing between points
!nmax=50                      ! Maximum order of besselfunction
ntot=step+1                   ! Number of grid points 

ALLOCATE(x(ntot))
ALLOCATE(z(mzro+1))
ALLOCATE(j(ntot,0:nmax))
ALLOCATE(dj(ntot,0:nmax))
ALLOCATE(y(ntot,0:nmax))
ALLOCATE(dy(ntot,0:nmax))
ALLOCATE(bess_n(ntot))
ALLOCATE(roots(mzro))	

! Make the x grid
   x(1)=left
   do i=2,ntot
      x(i)=x(i-1)+del
   end do

! Find values of regular and irregular bessel functions 
! and their derivatives on the grid
   call bessel(x,j,dj,y,dy,ntot,nmax,.false.)

! Find the first mzro roots of the bessel functions up to order nmax
   do n=0,nmax
      bess_n(:)=j(:,n)

! Locate intervals where bessel changes sign 
      call intrv(bess_n,z,left,del,ntot,mzro,cntzro,.false.)

      locz=1
      lroot=1
! Locate zeros of bessel
      do nroot=1,cntzro 
         call nwtrap(z(locz),z(locz+1),roots(lroot),convg,niter,number,n,nmax,'j')
         lroot=lroot+1
         locz=locz+1
      end do
! Store the roots
      besszro(:,n)=roots(:)
   end do

if(writeout) then
   Open(unit=1,file='Besselzeros.dat')
      write(1,*) 'Zeros of Bessel Functions'
      write(1,*)
      write(1,100) 
      100 format(3x,'s',9x,'j_0,s',13x,'j_1,s',13x,'j_2,s') 
      write(1,*)
      do i=1,20
         write(1,200) i, besszro(i,0), besszro(i,1), besszro(i,2) 
      end do
      write(1,*)
      write(1,*)
      write(1,300) 
      300 format(3x,'s',9x,'j_3,s',13x,'j_4,s',13x,'j_5,s') 
      write(1,*)
      do i=1,20
         write(1,200) i, besszro(i,3), besszro(i,4), besszro(i,5) 
      end do
      write(1,*)
      write(1,*)
      write(1,400) 
      400 format(3x,'s',9x,'j_6,s',13x,'j_7,s',13x,'j_8,s') 
      write(1,*)
      do i=1,20
         write(1,200) i, besszro(i,6), besszro(i,7), besszro(i,8) 
      end do
      200 format(i4,2x,f16.8,2x,f16.8,2x,f16.8)
   Close(1)


   Open(unit=2,file='bessel')
   do i=1,ntot
      write(2,500) x(i), j(i,10), j(i,20), j(i,30)
   end do
   500 format(f10.6,1x,e18.6,1x,e18.6,1x,e18.6,1x)
   Close(2)
end if

DEALLOCATE(x)
DEALLOCATE(z)
DEALLOCATE(j)
DEALLOCATE(dj)
DEALLOCATE(y)
DEALLOCATE(dy)
DEALLOCATE(bess_n)
DEALLOCATE(roots)

END subroutine besselzero


!deck bessel

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-20  Time: 15:30:03

!***begin prologue     bessel
!***date written       xxxxxx   (yymmdd)
!***revision date      890422   (yymmdd)
!***keywords           m6004, link 6004, bessel, spline
!***author             schneider, barry (lanl)
!***source             m6004
!***purpose            generate bessel functions for integer values
!***description
!***references

!***routines called
!***end prologue       bessel

SUBROUTINE bessel(x,j,dj,y,dy,np,nmax,prnt)

IMPLICIT INTEGER (a-z)
REAL*8, INTENT(IN)                       :: x(np)
REAL*8, INTENT(OUT)                      :: j(np,0:nmax)
REAL*8, INTENT(OUT)                      :: dj(np,0:nmax)
REAL*8, INTENT(OUT)                      :: y(np,0:nmax)
REAL*8, INTENT(OUT)                      :: dy(np,0:nmax)
REAL*8                                   :: factrl
INTEGER, INTENT(IN)                      :: np
INTEGER, INTENT(IN)                      :: nmax
LOGICAL, INTENT(IN)                      :: prnt
INTEGER, PARAMETER :: dim=1000
COMMON /io/ inp, iout
REAL*8 zero, one, two, epsilon
REAL*8 half, pi, twonpi, fournpi, xinv
REAL*8 norm, large, renorm
REAL*8 besj0, besj1, besy0, besy1, j0, j1, y0, y1
REAL*8 jnpl1, jn, jnm1

CHARACTER (LEN=80) :: title



DATA zero, one, two / 0.d+00, 1.d+00, 2.d+00 /
DATA half, large, renorm / .5D0, 1.d+250, 1.d-250  /
DATA pi / 3.141592653589793238462643D+00 /

twonpi=two/pi
fournpi=two*twonpi
!----------------------------------------------------------------------c
!            upward recursion, downward recursion or series            c
!----------------------------------------------------------------------c
DO  i=1,np
  xinv=1.d0/x(i)
  j0 = besj0(x(i))
  j1 = besj1(x(i))
  y0 = besy0(x(i))
  y1 = besy1(x(i))
  j(i,0) = j0
  y(i,0) = y0
  dj(i,0) = - j1
  dy(i,0) = - y1
  IF(nmax > 0) THEN
    j(i,1) = j1
    y(i,1) = y1
    dj(i,1) = j(i,0) - j(i,1)*xinv
    dy(i,1)= y(i,0) - y(i,1)*xinv
  END IF
  IF(nmax > 1) THEN
    
!           test for up or down recursion
    
    nn=x(i)
    IF(nn > nmax-1) THEN
      
!                  recur up
      
      j(i,0)=j0
      j(i,1)=j1
      nbeg=1
      DO  n=2,nmax
        j(i,n)=(nbeg+nbeg)*j(i,n-1)*xinv - j(i,n-2)
        nbeg=nbeg+1
      END DO
      
!              series or downward recursion
      
    ELSE
      IF(ABS(x(i)) > one) THEN
!----------------------------------------------------------------------c
!              find the value of upper needed to accurately get        c
!              the needed n values                                     c
!----------------------------------------------------------------------c
        upper=nmax
        upper=msta1(x(i),200)
        IF(upper < nmax) THEN
          upper=nmax + nmax
        ELSE
          upper=msta2(x(i),nmax,15)
        END IF
        upper=MAX(upper,nmax+2)
        jnpl1=zero
        jn=renorm
        onelss=upper-1
        nbeg=onelss
        DO  n=onelss-1,0,-1
          jnm1 = ( nbeg + nbeg )*jn*xinv -  &
              jnpl1
          IF(ABS(jnm1) > large) THEN
            jnm1=jnm1*renorm
            jn=jn*renorm
          END IF
          IF(n <= nmax) THEN
            j(i,n) = jnm1
          END IF
          jnpl1=jn
          jn=jnm1
          nbeg = nbeg -1
        END DO
!----------------------------------------------------------------------c
!                 normalize the j                                      c
!----------------------------------------------------------------------c
        IF(ABS(j0) > ABS(j1)) THEN
          norm=j0/j(i,0)
        ELSE
          norm=j1/j(i,1)
        END IF
        DO  n=0,nmax
          j(i,n) = norm*j(i,n)
        END DO
      ELSE
        DO  n=2,nmax
          CALL bser(x(i),j(i,n),dj(i,n),n,0)
        END DO
      END IF
    END IF
!----------------------------------------------------------------------c
!           recur upward for y(i,n)                                    c
!----------------------------------------------------------------------c
    nbeg=1
    DO  n=2,nmax
      y(i,n) = ( nbeg + nbeg )*y(i,n-1)*xinv - y(i,n-2)
      IF(ABS(y(i,n)) > large) THEN
        y(i,n)=large
      END IF
      nbeg=nbeg+1
    END DO
!----------------------------------------------------------------------c
!               get derivatives                                        c
!----------------------------------------------------------------------c
    DO  n=2,nmax
      dj(i,n) = j(i,n-1) - n*j(i,n)*xinv
      dy(i,n) = y(i,n-1) - n*y(i,n)*xinv
    END DO
  END IF
END DO
IF (prnt) THEN
  title='regular bessel functions'
  CALL prntfm(title,j,np,nmax+1,np,nmax+1,iout)
  title='irrregular bessel functions'
  CALL prntfm(title,y,np,nmax+1,np,nmax+1,iout)
  title='derivative of regular bessel functions'
  CALL prntfm(title,dj,np,nmax+1,np,nmax+1,iout)
  title='derivative of irregular bessel functions'
  CALL prntfm(title,dy,np,nmax+1,np,nmax+1,iout)
END IF
RETURN
END SUBROUTINE bessel


!deck nwtrap

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-20  Time: 15:30:07

!***begin prologue     nwtrap
!***date written       930623   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           newton-raphson root finder
!***author             schneider, barry (nsf)
!***source             math
!***purpose            find zero of function using combination of bisection
!***                   and newton-raphson method
!***
!***

!***references

!***routines called
!***end prologue       nwtrap

SUBROUTINE nwtrap(x1,x2,rtsafe,convg,niter,nactul,l,ltop,TYPE)

IMPLICIT INTEGER (a-z)
REAL*8, INTENT(IN)                       :: x1
REAL*8, INTENT(IN)                       :: x2
REAL*8, INTENT(OUT)                      :: rtsafe
REAL*8, INTENT(IN)                       :: convg
INTEGER, INTENT(IN)                      :: niter
INTEGER, INTENT(OUT)                     :: nactul
REAL*8                                   :: f1(0:ltop)
REAL*8                                   :: f2(0:ltop)
REAL*8                                   :: f3(0:ltop)
REAL*8                                   :: f4(0:ltop)
INTEGER, INTENT(IN OUT)                  :: l
INTEGER, INTENT(IN OUT)                  :: ltop
CHARACTER (LEN=*), INTENT(IN)            :: TYPE
REAL*8  fl, df, fh, xl, xh, dxold
REAL*8 dx, f, temp


COMMON/io/inp, iout

CALL bessel(x1,f1,f2,f3,f4,1,ltop,.false.)
fl=f1(l)
IF(TYPE == 'jp') THEN
  fl=f2(l)
END IF
CALL bessel(x2,f1,f2,f3,f4,1,ltop,.false.)
fh=f1(l) 
IF(TYPE == 'jp') THEN
  fh=f2(l)
END IF

IF (fl*fh >= 0.d0) THEN
  WRITE (iout,1) x1, fl, x2, fh
  CALL lnkerr('quit')
END IF
IF (fl < 0.d0) THEN
  xl=x1
  xh=x2
ELSE
  xl=x2
  xh=x1
END IF
rtsafe=.5D0*(x1+x2)
dxold=ABS(x2-x1)
dx=dxold
CALL bessel(rtsafe,f1,f2,f3,f4,1,ltop,.false.)
f=f1(l)
df=f2(l)
!IF (TYPE == 'jp') THEN
!  f=f2(l)
!  CALL secder(f1(l),df,rtsafe,l,1)
!END IF
DO  iter=1,niter
  IF (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) >= 0.d0.OR.  &
        ABS(2.d0*f) > ABS(dxold*df)) THEN
    dxold=dx
    dx=.5D0*(xh-xl)
    rtsafe=xl+dx
    IF ( xl == rtsafe ) THEN
      nactul=iter
      RETURN
    END IF
  ELSE
    dxold=dx
    dx=f/df
    temp=rtsafe
    rtsafe=rtsafe-dx
    IF ( temp == rtsafe ) THEN
      nactul=iter
      RETURN
    END IF
  END IF
  IF (ABS(dx) < convg) THEN
    nactul=iter
    RETURN
  END IF
  CALL bessel(rtsafe,f1,f2,f3,f4,1,ltop,.false.)
  f=f1(l)
  df=f2(l)
!  IF (TYPE == 'jp') THEN
!    f=f2(l)
!    CALL secder(f1(l),df,rtsafe,l,1)
!  END IF
  IF (f < 0.d0) THEN
    xl=rtsafe
  ELSE
    xh=rtsafe
  END IF
END DO
WRITE (iout,2) niter
RETURN
1 FORMAT (/,5X,'root not bracketed. will quit.',/,5X,  &
    'xl = ',e15.8,1X,'fl = ',e15.8/,5X, 'xr = ',e15.8,1X,'fr = ',e15.8)
2 FORMAT (/,5X,'no convergence after ',i4,' iterations')
END SUBROUTINE nwtrap








!deck besj0

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-20  Time: 16:01:32

FUNCTION besj0 (x)
!***begin prologue  besj0
!***purpose  compute the bessel function of the first kind of order
!            zero.
!***library   slatec (fnlib)
!***category  c10a1
!***type      double precision
!***keywords  bessel function, first kind, fnlib, order zero,
!             special functions
!***author  fullerton, w., (lanl)
!***description

! besj0(x) calculates the double precision bessel function of
! the first kind of order zero for double precision argument x.

! series for bj0        on the interval  0.          to  1.60000e+01
!                                        with weighted error   4.39e-32
!                                         log weighted error  31.36
!                               significant figures required  31.21
!                                    decimal places required  32.00

!***references  (none)
!***routines called  r1mach, r9b0mp, csevl, initds
!***revision history  (yymmdd)
!   770701  date written
!   890531  changed all specific intrinsics to generic.  (wrb)
!   890531  revision date from version 3.2
!   891214  prologue converted to version 4.0 format.  (bab)
!***end prologue  besj0

REAL*8, INTENT(IN OUT)                   :: x
REAL*8  bj0cs(19), ampl, theta, xsml, y, r1mach, csevl
LOGICAL :: first
SAVE bj0cs, ntj0, xsml, first
DATA bj0cs(  1) / +.10025416196893913701073127264074d+0     /
DATA bj0cs(  2) / -.66522300776440513177678757831124d+0     /
DATA bj0cs(  3) / +.24898370349828131370460468726680d+0     /
DATA bj0cs(  4) / -.33252723170035769653884341503854d-1     /
DATA bj0cs(  5) / +.23114179304694015462904924117729d-2     /
DATA bj0cs(  6) / -.99112774199508092339048519336549d-4     /
DATA bj0cs(  7) / +.28916708643998808884733903747078d-5     /
DATA bj0cs(  8) / -.61210858663032635057818407481516d-7     /
DATA bj0cs(  9) / +.98386507938567841324768748636415d-9     /
DATA bj0cs( 10) / -.12423551597301765145515897006836d-10    /
DATA bj0cs( 11) / +.12654336302559045797915827210363d-12    /
DATA bj0cs( 12) / -.10619456495287244546914817512959d-14    /
DATA bj0cs( 13) / +.74706210758024567437098915584000d-17    /
DATA bj0cs( 14) / -.44697032274412780547627007999999d-19    /
DATA bj0cs( 15) / +.23024281584337436200523093333333d-21    /
DATA bj0cs( 16) / -.10319144794166698148522666666666d-23    /
DATA bj0cs( 17) / +.40608178274873322700800000000000d-26    /
DATA bj0cs( 18) / -.14143836005240913919999999999999d-28    /
DATA bj0cs( 19) / +.43910905496698880000000000000000d-31    /
DATA first /.true./
!***first executable statement  besj0
IF (first) THEN
  ntj0 = initds (bj0cs, 19, 0.1*REAL(r1mach(3)))
  xsml = SQRT(8.0D0*r1mach(3))
END IF
first = .false.

y = ABS(x)
IF (y > 4.0D0) GO TO 20

besj0 = 1.0D0
IF (y > xsml) besj0 = csevl (.125D0*y*y-1.d0, bj0cs, ntj0)
RETURN

20   CALL d9b0mp (y, ampl, theta)
besj0 = ampl * COS(theta)

RETURN
END FUNCTION besj0


!deck besj1

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-20  Time: 16:01:39

FUNCTION besj1 (x)
!***begin prologue  dbesj1
!***purpose  compute the bessel function of the first kind of order one.
!***library   slatec (fnlib)
!***category  c10a1
!***type      double precision (besj1-s, dbesj1-d)
!***keywords  bessel function, first kind, fnlib, order one,
!             special functions
!***author  fullerton, w., (lanl)
!***description

! dbesj1(x) calculates the double precision bessel function of the
! first kind of order one for double precision argument x.

! series for bj1        on the interval  0.          to  1.60000e+01
!                                        with weighted error   1.16e-33
!                                         log weighted error  32.93
!                               significant figures required  32.36
!                                    decimal places required  33.57

!***references  (none)
!***routines called  d1mach, d9b1mp, dcsevl, initds, xermsg
!***revision history  (yymmdd)
!   780601  date written
!   890531  changed all specific intrinsics to generic.  (wrb)
!   890531  revision date from version 3.2
!   891214  prologue converted to version 4.0 format.  (bab)
!   900315  calls to xerror changed to calls to xermsg.  (thj)
!   910401  corrected error in code which caused values to have the
!           wrong sign for arguments less than 4.0.  (wrb)
!***end prologue  dbesj1

REAL*8, INTENT(IN)                       :: x
REAL*8  bj1cs(19), ampl, theta, xsml, xmin, y, r1mach, csevl
LOGICAL :: first
SAVE bj1cs, ntj1, xsml, xmin, first
DATA bj1cs(  1) / -.117261415133327865606240574524003d+0    /
DATA bj1cs(  2) / -.253615218307906395623030884554698d+0    /
DATA bj1cs(  3) / +.501270809844695685053656363203743d-1    /
DATA bj1cs(  4) / -.463151480962508191842619728789772d-2    /
DATA bj1cs(  5) / +.247996229415914024539124064592364d-3    /
DATA bj1cs(  6) / -.867894868627882584521246435176416d-5    /
DATA bj1cs(  7) / +.214293917143793691502766250991292d-6    /
DATA bj1cs(  8) / -.393609307918317979229322764073061d-8    /
DATA bj1cs(  9) / +.559118231794688004018248059864032d-10   /
DATA bj1cs( 10) / -.632761640466139302477695274014880d-12   /
DATA bj1cs( 11) / +.584099161085724700326945563268266d-14   /
DATA bj1cs( 12) / -.448253381870125819039135059199999d-16   /
DATA bj1cs( 13) / +.290538449262502466306018688000000d-18   /
DATA bj1cs( 14) / -.161173219784144165412118186666666d-20   /
DATA bj1cs( 15) / +.773947881939274637298346666666666d-23   /
DATA bj1cs( 16) / -.324869378211199841143466666666666d-25   /
DATA bj1cs( 17) / +.120223767722741022720000000000000d-27   /
DATA bj1cs( 18) / -.395201221265134933333333333333333d-30   /
DATA bj1cs( 19) / +.116167808226645333333333333333333d-32   /
DATA first /.true./
!***first executable statement  dbesj1
IF (first) THEN
  ntj1 = initds (bj1cs, 19, 0.1*REAL(r1mach(3)))
  
  xsml = SQRT(8.0D0*r1mach(3))
  xmin = 2.0D0*r1mach(1)
END IF
first = .false.

y = ABS(x)
IF (y > 4.0D0) GO TO 20

besj1 = 0.0D0
IF (y == 0.0D0) RETURN
IF (y <= xmin) CALL lnkerr ('besj1:abs(x) j1 underflows')
IF (y > xmin) besj1 = 0.5D0*x
IF (y > xsml) besj1 = x*(.25D0 + csevl (.125D0*y*y-1.d0, bj1cs, ntj1) )
RETURN

20   CALL d9b1mp (y, ampl, theta)
besj1 = SIGN (ampl, x) * COS(theta)

RETURN
END FUNCTION besj1


!deck besy0

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-20  Time: 16:01:43

FUNCTION besy0 (x)
!***begin prologue  dbesy0
!***purpose  compute the bessel function of the second kind of order
!            zero.
!***library   slatec (fnlib)
!***category  c10a1
!***type      double precision (besy0-s, dbesy0-d)
!***keywords  bessel function, fnlib, order zero, second kind,
!             special functions
!***author  fullerton, w., (lanl)
!***description

! dbesy0(x) calculates the double precision bessel function of the
! second kind of order zero for double precision argument x.

! series for by0        on the interval  0.          to  1.60000e+01
!                                        with weighted error   8.14e-32
!                                         log weighted error  31.09
!                               significant figures required  30.31
!                                    decimal places required  31.73

!***references  (none)
!***routines called  d1mach, d9b0mp, dbesj0, dcsevl, initds, xermsg
!***revision history  (yymmdd)
!   770701  date written
!   890531  changed all specific intrinsics to generic.  (wrb)
!   890531  revision date from version 3.2
!   891214  prologue converted to version 4.0 format.  (bab)
!   900315  calls to xerror changed to calls to xermsg.  (thj)
!***end prologue  dbesy0

REAL*8, INTENT(IN)                       :: x
REAL*8  by0cs(19), ampl, theta, twodpi, xsml, y, r1mach, csevl, besj0
LOGICAL :: first
SAVE by0cs, twodpi, nty0, xsml, first
DATA by0cs(  1) / -.1127783939286557321793980546028d-1      /
DATA by0cs(  2) / -.1283452375604203460480884531838d+0      /
DATA by0cs(  3) / -.1043788479979424936581762276618d+0      /
DATA by0cs(  4) / +.2366274918396969540924159264613d-1      /
DATA by0cs(  5) / -.2090391647700486239196223950342d-2      /
DATA by0cs(  6) / +.1039754539390572520999246576381d-3      /
DATA by0cs(  7) / -.3369747162423972096718775345037d-5      /
DATA by0cs(  8) / +.7729384267670667158521367216371d-7      /
DATA by0cs(  9) / -.1324976772664259591443476068964d-8      /
DATA by0cs( 10) / +.1764823261540452792100389363158d-10     /
DATA by0cs( 11) / -.1881055071580196200602823012069d-12     /
DATA by0cs( 12) / +.1641865485366149502792237185749d-14     /
DATA by0cs( 13) / -.1195659438604606085745991006720d-16     /
DATA by0cs( 14) / +.7377296297440185842494112426666d-19     /
DATA by0cs( 15) / -.3906843476710437330740906666666d-21     /
DATA by0cs( 16) / +.1795503664436157949829120000000d-23     /
DATA by0cs( 17) / -.7229627125448010478933333333333d-26     /
DATA by0cs( 18) / +.2571727931635168597333333333333d-28     /
DATA by0cs( 19) / -.8141268814163694933333333333333d-31     /
DATA twodpi / 0.636619772367581343075535053490057d0 /
DATA first /.true./
!***first executable statement  besy0
IF (first) THEN
  nty0 = initds (by0cs, 19, 0.1*REAL(r1mach(3)))
  xsml = SQRT(4.0D0*r1mach(3))
END IF
first = .false.

IF (x <= 0.d0) CALL lnkerr ('dbesy0:x zero or negative')
IF (x > 4.0D0) GO TO 20

y = 0.d0
IF (x > xsml) y = x*x
besy0 = twodpi*LOG(0.5D0*x)*besj0(x) + .375D0 + csevl (  &
    .125D0*y-1.d0, by0cs, nty0)
RETURN

20   CALL d9b0mp (x, ampl, theta)
besy0 = ampl * SIN(theta)
RETURN

END FUNCTION besy0


!deck besy1

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-20  Time: 16:01:46

FUNCTION besy1 (x)
!***begin prologue  dbesy1
!***purpose  compute the bessel function of the second kind of order
!            one.
!***library   slatec (fnlib)
!***category  c10a1
!***type      double precision (besy1-s, dbesy1-d)
!***keywords  bessel function, fnlib, order one, second kind,
!             special functions
!***author  fullerton, w., (lanl)
!***description

! dbesy1(x) calculates the double precision bessel function of the
! second kind of order for double precision argument x.

! series for by1        on the interval  0.          to  1.60000e+01
!                                        with weighted error   8.65e-33
!                                         log weighted error  32.06
!                               significant figures required  32.17
!                                    decimal places required  32.71

!***references  (none)
!***routines called  d1mach, d9b1mp, dbesj1, dcsevl, initds, xermsg
!***revision history  (yymmdd)
!   770701  date written
!   890531  changed all specific intrinsics to generic.  (wrb)
!   890531  revision date from version 3.2
!   891214  prologue converted to version 4.0 format.  (bab)
!   900315  calls to xerror changed to calls to xermsg.  (thj)
!***end prologue  dbesy1

REAL*8, INTENT(IN)                       :: x
REAL*8  by1cs(20), ampl, theta, twodpi, xmin, xsml, y, r1mach, csevl, besj1
LOGICAL :: first
SAVE by1cs, twodpi, nty1, xmin, xsml, first
DATA by1cs(  1) / +.320804710061190862932352018628015d-1    /
DATA by1cs(  2) / +.126270789743350044953431725999727d+1    /
DATA by1cs(  3) / +.649996189992317500097490637314144d-2    /
DATA by1cs(  4) / -.893616452886050411653144160009712d-1    /
DATA by1cs(  5) / +.132508812217570954512375510370043d-1    /
DATA by1cs(  6) / -.897905911964835237753039508298105d-3    /
DATA by1cs(  7) / +.364736148795830678242287368165349d-4    /
DATA by1cs(  8) / -.100137438166600055549075523845295d-5    /
DATA by1cs(  9) / +.199453965739017397031159372421243d-7    /
DATA by1cs( 10) / -.302306560180338167284799332520743d-9    /
DATA by1cs( 11) / +.360987815694781196116252914242474d-11   /
DATA by1cs( 12) / -.348748829728758242414552947409066d-13   /
DATA by1cs( 13) / +.278387897155917665813507698517333d-15   /
DATA by1cs( 14) / -.186787096861948768766825352533333d-17   /
DATA by1cs( 15) / +.106853153391168259757070336000000d-19   /
DATA by1cs( 16) / -.527472195668448228943872000000000d-22   /
DATA by1cs( 17) / +.227019940315566414370133333333333d-24   /
DATA by1cs( 18) / -.859539035394523108693333333333333d-27   /
DATA by1cs( 19) / +.288540437983379456000000000000000d-29   /
DATA by1cs( 20) / -.864754113893717333333333333333333d-32   /
DATA twodpi / 0.636619772367581343075535053490057d0 /
DATA first /.true./
!***first executable statement  dbesy1
IF (first) THEN
  nty1 = initds (by1cs, 20, 0.1*REAL(r1mach(3)))
  
  xmin = 1.571D0 * EXP (MAX(LOG(r1mach(1)), -LOG(r1mach(2))) + 0.01D0)
  xsml = SQRT(4.0D0*r1mach(3))
END IF
first = .false.

IF (x <= 0.d0) CALL lnkerr ('besy1:x zero or negative')
IF (x > 4.0D0) GO TO 20

IF (x < xmin) CALL lnkerr ('besy1:x small y1 overflows')
y = 0.d0
IF (x > xsml) y = x*x
besy1 = twodpi * LOG(0.5D0*x)*besj1(x) + (0.5D0 +  &
    csevl (.125D0*y-1.d0, by1cs, nty1))/x
RETURN

20   CALL d9b1mp (x, ampl, theta)
besy1 = ampl * SIN(theta)
RETURN

END FUNCTION besy1



!deck bser

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-20  Time: 15:29:59

!***begin prologue     bser
!***date written       xxxxxx   (yymmdd)
!***revision date      890422   (yymmdd)
!***keywords           m6004, link 6004, bessel, spline
!***author             schneider, barry (lanl)
!***source             m6004
!***purpose            generate regular bessel functions for small argument
!***description        using series.
!***references

!***routines called
!***end prologue       bser

SUBROUTINE bser(x,j,dj,m,n)

IMPLICIT INTEGER (a-z)
REAL*8, INTENT(IN)                       :: x
REAL*8, INTENT(OUT)                      :: j
REAL*8, INTENT(OUT)                      :: dj
REAL*8                                   :: factrl
INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN)                      :: n
COMMON /io/ inp, iout

REAL*8 xfac, xfac2, pre, prex, term, sum, dsum
REAL*8 xfacm, zero, one, half, cnvrg
LOGICAL :: prnt
CHARACTER (LEN=80) :: title

DATA zero, one, half, cnvrg / 0.d0, 1.d0, .5D0, 1.d-20 /

mn=m-n
IF(mn < 0) THEN
  WRITE(iout,1)
  RETURN
END IF
xfac=x*half
xfacm=one/xfac
xfac2=xfac*xfac
pre=one/2**n
sum=zero
dsum=sum
IF(mn /= 0) THEN
  prex=xfac**mn
ELSE
  prex=one
END IF
imul=1
DO  i=0,100
  IF(prex <= cnvrg) THEN
    GO TO 20
  END IF
  term=imul*prex/(factrl(i)*factrl(m+i))
  sum=sum+term
  dsum=dsum+(m-n+i+i)*term*xfacm
  imul=-imul
  prex=prex*xfac2
END DO
WRITE(iout,2)
CALL lnkerr('quit')
20   j=pre*sum
dj=pre*dsum*half
RETURN
1    FORMAT(/,1X,'power of x must be ge 0')
2    FORMAT(/,1X,'no series convergence')
END SUBROUTINE bser


!deck d9b0mp

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-21  Time: 16:00:27

SUBROUTINE d9b0mp (x, ampl, theta)
!***begin prologue  d9b0mp
!***subsidiary
!***purpose  evaluate the modulus and phase for the j0 and y0 bessel
!            functions.
!***library   slatec (fnlib)
!***category  c10a1
!***type      double precision (d9b0mp-d)
!***keywords  bessel function, fnlib, modulus, phase, special functions
!***author  fullerton, w., (lanl)
!***description

! evaluate the modulus and phase for the bessel j0 and y0 functions.

! series for bm0        on the interval  1.56250e-02 to  6.25000e-02
!                                        with weighted error   4.40e-32
!                                         log weighted error  31.36
!                               significant figures required  30.02
!                                    decimal places required  32.14

! series for bth0       on the interval  0.          to  1.56250e-02
!                                        with weighted error   2.66e-32
!                                         log weighted error  31.57
!                               significant figures required  30.67
!                                    decimal places required  32.40

! series for bm02       on the interval  0.          to  1.56250e-02
!                                        with weighted error   4.72e-32
!                                         log weighted error  31.33
!                               significant figures required  30.00
!                                    decimal places required  32.13

! series for bt02       on the interval  1.56250e-02 to  6.25000e-02
!                                        with weighted error   2.99e-32
!                                         log weighted error  31.52
!                               significant figures required  30.61
!                                    decimal places required  32.32

!***references  (none)
!***routines called  d1mach, dcsevl, initds, xermsg
!***revision history  (yymmdd)
!   770701  date written
!   890531  changed all specific intrinsics to generic.  (wrb)
!   890531  revision date from version 3.2
!   891214  prologue converted to version 4.0 format.  (bab)
!   900315  calls to xerror changed to calls to xermsg.  (thj)
!   900720  routine changed from user-callable to subsidiary.  (wrb)
!   920618  removed space from variable names.  (rwc, wrb)
!***end prologue  d9b0mp

REAL*8, INTENT(IN)                       :: x
REAL*8, INTENT(OUT)                      :: ampl
REAL*8, INTENT(OUT)                      :: theta
REAL*8  bm0cs(37), bt02cs(39),  &
    bm02cs(40), bth0cs(44), xmax, pi4, z, r1mach, csevl
LOGICAL :: first
SAVE bm0cs, bth0cs, bm02cs, bt02cs, pi4, nbm0, nbt02,  &
    nbm02, nbth0, xmax, first
DATA bm0cs(  1) / +.9211656246827742712573767730182d-1      /
DATA bm0cs(  2) / -.1050590997271905102480716371755d-2      /
DATA bm0cs(  3) / +.1470159840768759754056392850952d-4      /
DATA bm0cs(  4) / -.5058557606038554223347929327702d-6      /
DATA bm0cs(  5) / +.2787254538632444176630356137881d-7      /
DATA bm0cs(  6) / -.2062363611780914802618841018973d-8      /
DATA bm0cs(  7) / +.1870214313138879675138172596261d-9      /
DATA bm0cs(  8) / -.1969330971135636200241730777825d-10     /
DATA bm0cs(  9) / +.2325973793999275444012508818052d-11     /
DATA bm0cs( 10) / -.3009520344938250272851224734482d-12     /
DATA bm0cs( 11) / +.4194521333850669181471206768646d-13     /
DATA bm0cs( 12) / -.6219449312188445825973267429564d-14     /
DATA bm0cs( 13) / +.9718260411336068469601765885269d-15     /
DATA bm0cs( 14) / -.1588478585701075207366635966937d-15     /
DATA bm0cs( 15) / +.2700072193671308890086217324458d-16     /
DATA bm0cs( 16) / -.4750092365234008992477504786773d-17     /
DATA bm0cs( 17) / +.8615128162604370873191703746560d-18     /
DATA bm0cs( 18) / -.1605608686956144815745602703359d-18     /
DATA bm0cs( 19) / +.3066513987314482975188539801599d-19     /
DATA bm0cs( 20) / -.5987764223193956430696505617066d-20     /
DATA bm0cs( 21) / +.1192971253748248306489069841066d-20     /
DATA bm0cs( 22) / -.2420969142044805489484682581333d-21     /
DATA bm0cs( 23) / +.4996751760510616453371002879999d-22     /
DATA bm0cs( 24) / -.1047493639351158510095040511999d-22     /
DATA bm0cs( 25) / +.2227786843797468101048183466666d-23     /
DATA bm0cs( 26) / -.4801813239398162862370542933333d-24     /
DATA bm0cs( 27) / +.1047962723470959956476996266666d-24     /
DATA bm0cs( 28) / -.2313858165678615325101260800000d-25     /
DATA bm0cs( 29) / +.5164823088462674211635199999999d-26     /
DATA bm0cs( 30) / -.1164691191850065389525401599999d-26     /
DATA bm0cs( 31) / +.2651788486043319282958336000000d-27     /
DATA bm0cs( 32) / -.6092559503825728497691306666666d-28     /
DATA bm0cs( 33) / +.1411804686144259308038826666666d-28     /
DATA bm0cs( 34) / -.3298094961231737245750613333333d-29     /
DATA bm0cs( 35) / +.7763931143074065031714133333333d-30     /
DATA bm0cs( 36) / -.1841031343661458478421333333333d-30     /
DATA bm0cs( 37) / +.4395880138594310737100799999999d-31     /
DATA bth0cs(  1) / -.24901780862128936717709793789967d+0     /
DATA bth0cs(  2) / +.48550299609623749241048615535485d-3     /
DATA bth0cs(  3) / -.54511837345017204950656273563505d-5     /
DATA bth0cs(  4) / +.13558673059405964054377445929903d-6     /
DATA bth0cs(  5) / -.55691398902227626227583218414920d-8     /
DATA bth0cs(  6) / +.32609031824994335304004205719468d-9     /
DATA bth0cs(  7) / -.24918807862461341125237903877993d-10    /
DATA bth0cs(  8) / +.23449377420882520554352413564891d-11    /
DATA bth0cs(  9) / -.26096534444310387762177574766136d-12    /
DATA bth0cs( 10) / +.33353140420097395105869955014923d-13    /
DATA bth0cs( 11) / -.47890000440572684646750770557409d-14    /
DATA bth0cs( 12) / +.75956178436192215972642568545248d-15    /
DATA bth0cs( 13) / -.13131556016891440382773397487633d-15    /
DATA bth0cs( 14) / +.24483618345240857495426820738355d-16    /
DATA bth0cs( 15) / -.48805729810618777683256761918331d-17    /
DATA bth0cs( 16) / +.10327285029786316149223756361204d-17    /
DATA bth0cs( 17) / -.23057633815057217157004744527025d-18    /
DATA bth0cs( 18) / +.54044443001892693993017108483765d-19    /
DATA bth0cs( 19) / -.13240695194366572724155032882385d-19    /
DATA bth0cs( 20) / +.33780795621371970203424792124722d-20    /
DATA bth0cs( 21) / -.89457629157111779003026926292299d-21    /
DATA bth0cs( 22) / +.24519906889219317090899908651405d-21    /
DATA bth0cs( 23) / -.69388422876866318680139933157657d-22    /
DATA bth0cs( 24) / +.20228278714890138392946303337791d-22    /
DATA bth0cs( 25) / -.60628500002335483105794195371764d-23    /
DATA bth0cs( 26) / +.18649748964037635381823788396270d-23    /
DATA bth0cs( 27) / -.58783732384849894560245036530867d-24    /
DATA bth0cs( 28) / +.18958591447999563485531179503513d-24    /
DATA bth0cs( 29) / -.62481979372258858959291620728565d-25    /
DATA bth0cs( 30) / +.21017901684551024686638633529074d-25    /
DATA bth0cs( 31) / -.72084300935209253690813933992446d-26    /
DATA bth0cs( 32) / +.25181363892474240867156405976746d-26    /
DATA bth0cs( 33) / -.89518042258785778806143945953643d-27    /
DATA bth0cs( 34) / +.32357237479762298533256235868587d-27    /
DATA bth0cs( 35) / -.11883010519855353657047144113796d-27    /
DATA bth0cs( 36) / +.44306286907358104820579231941731d-28    /
DATA bth0cs( 37) / -.16761009648834829495792010135681d-28    /
DATA bth0cs( 38) / +.64292946921207466972532393966088d-29    /
DATA bth0cs( 39) / -.24992261166978652421207213682763d-29    /
DATA bth0cs( 40) / +.98399794299521955672828260355318d-30    /
DATA bth0cs( 41) / -.39220375242408016397989131626158d-30    /
DATA bth0cs( 42) / +.15818107030056522138590618845692d-30    /
DATA bth0cs( 43) / -.64525506144890715944344098365426d-31    /
DATA bth0cs( 44) / +.26611111369199356137177018346367d-31    /
DATA bm02cs(  1) / +.9500415145228381369330861335560d-1      /
DATA bm02cs(  2) / -.3801864682365670991748081566851d-3      /
DATA bm02cs(  3) / +.2258339301031481192951829927224d-5      /
DATA bm02cs(  4) / -.3895725802372228764730621412605d-7      /
DATA bm02cs(  5) / +.1246886416512081697930990529725d-8      /
DATA bm02cs(  6) / -.6065949022102503779803835058387d-10     /
DATA bm02cs(  7) / +.4008461651421746991015275971045d-11     /
DATA bm02cs(  8) / -.3350998183398094218467298794574d-12     /
DATA bm02cs(  9) / +.3377119716517417367063264341996d-13     /
DATA bm02cs( 10) / -.3964585901635012700569356295823d-14     /
DATA bm02cs( 11) / +.5286111503883857217387939744735d-15     /
DATA bm02cs( 12) / -.7852519083450852313654640243493d-16     /
DATA bm02cs( 13) / +.1280300573386682201011634073449d-16     /
DATA bm02cs( 14) / -.2263996296391429776287099244884d-17     /
DATA bm02cs( 15) / +.4300496929656790388646410290477d-18     /
DATA bm02cs( 16) / -.8705749805132587079747535451455d-19     /
DATA bm02cs( 17) / +.1865862713962095141181442772050d-19     /
DATA bm02cs( 18) / -.4210482486093065457345086972301d-20     /
DATA bm02cs( 19) / +.9956676964228400991581627417842d-21     /
DATA bm02cs( 20) / -.2457357442805313359605921478547d-21     /
DATA bm02cs( 21) / +.6307692160762031568087353707059d-22     /
DATA bm02cs( 22) / -.1678773691440740142693331172388d-22     /
DATA bm02cs( 23) / +.4620259064673904433770878136087d-23     /
DATA bm02cs( 24) / -.1311782266860308732237693402496d-23     /
DATA bm02cs( 25) / +.3834087564116302827747922440276d-24     /
DATA bm02cs( 26) / -.1151459324077741271072613293576d-24     /
DATA bm02cs( 27) / +.3547210007523338523076971345213d-25     /
DATA bm02cs( 28) / -.1119218385815004646264355942176d-25     /
DATA bm02cs( 29) / +.3611879427629837831698404994257d-26     /
DATA bm02cs( 30) / -.1190687765913333150092641762463d-26     /
DATA bm02cs( 31) / +.4005094059403968131802476449536d-27     /
DATA bm02cs( 32) / -.1373169422452212390595193916017d-27     /
DATA bm02cs( 33) / +.4794199088742531585996491526437d-28     /
DATA bm02cs( 34) / -.1702965627624109584006994476452d-28     /
DATA bm02cs( 35) / +.6149512428936330071503575161324d-29     /
DATA bm02cs( 36) / -.2255766896581828349944300237242d-29     /
DATA bm02cs( 37) / +.8399707509294299486061658353200d-30     /
DATA bm02cs( 38) / -.3172997595562602355567423936152d-30     /
DATA bm02cs( 39) / +.1215205298881298554583333026514d-30     /
DATA bm02cs( 40) / -.4715852749754438693013210568045d-31     /
DATA bt02cs(  1) / -.24548295213424597462050467249324d+0     /
DATA bt02cs(  2) / +.12544121039084615780785331778299d-2     /
DATA bt02cs(  3) / -.31253950414871522854973446709571d-4     /
DATA bt02cs(  4) / +.14709778249940831164453426969314d-5     /
DATA bt02cs(  5) / -.99543488937950033643468850351158d-7     /
DATA bt02cs(  6) / +.85493166733203041247578711397751d-8     /
DATA bt02cs(  7) / -.86989759526554334557985512179192d-9     /
DATA bt02cs(  8) / +.10052099533559791084540101082153d-9     /
DATA bt02cs(  9) / -.12828230601708892903483623685544d-10    /
DATA bt02cs( 10) / +.17731700781805131705655750451023d-11    /
DATA bt02cs( 11) / -.26174574569485577488636284180925d-12    /
DATA bt02cs( 12) / +.40828351389972059621966481221103d-13    /
DATA bt02cs( 13) / -.66751668239742720054606749554261d-14    /
DATA bt02cs( 14) / +.11365761393071629448392469549951d-14    /
DATA bt02cs( 15) / -.20051189620647160250559266412117d-15    /
DATA bt02cs( 16) / +.36497978794766269635720591464106d-16    /
DATA bt02cs( 17) / -.68309637564582303169355843788800d-17    /
DATA bt02cs( 18) / +.13107583145670756620057104267946d-17    /
DATA bt02cs( 19) / -.25723363101850607778757130649599d-18    /
DATA bt02cs( 20) / +.51521657441863959925267780949333d-19    /
DATA bt02cs( 21) / -.10513017563758802637940741461333d-19    /
DATA bt02cs( 22) / +.21820381991194813847301084501333d-20    /
DATA bt02cs( 23) / -.46004701210362160577225905493333d-21    /
DATA bt02cs( 24) / +.98407006925466818520953651199999d-22    /
DATA bt02cs( 25) / -.21334038035728375844735986346666d-22    /
DATA bt02cs( 26) / +.46831036423973365296066286933333d-23    /
DATA bt02cs( 27) / -.10400213691985747236513382399999d-23    /
DATA bt02cs( 28) / +.23349105677301510051777740800000d-24    /
DATA bt02cs( 29) / -.52956825323318615788049749333333d-25    /
DATA bt02cs( 30) / +.12126341952959756829196287999999d-25    /
DATA bt02cs( 31) / -.28018897082289428760275626666666d-26    /
DATA bt02cs( 32) / +.65292678987012873342593706666666d-27    /
DATA bt02cs( 33) / -.15337980061873346427835733333333d-27    /
DATA bt02cs( 34) / +.36305884306364536682359466666666d-28    /
DATA bt02cs( 35) / -.86560755713629122479172266666666d-29    /
DATA bt02cs( 36) / +.20779909972536284571238399999999d-29    /
DATA bt02cs( 37) / -.50211170221417221674325333333333d-30    /
DATA bt02cs( 38) / +.12208360279441714184191999999999d-30    /
DATA bt02cs( 39) / -.29860056267039913454250666666666d-31    /
DATA pi4 / 0.785398163397448309615660845819876d0 /
DATA first /.true./
!***first executable statement  d9b0mp
IF (first) THEN
  eta = 0.1*REAL(r1mach(3))
  nbm0 = initds (bm0cs, 37, eta)
  nbt02 = initds (bt02cs, 39, eta)
  nbm02 = initds (bm02cs, 40, eta)
  nbth0 = initds (bth0cs, 44, eta)
  
  xmax = 1.0D0/r1mach(4)
END IF
first = .false.

IF (x < 4.d0) CALL lnkerr ('d9b0mp:x must be ge 4')

IF (x > 8.d0) GO TO 20
z = (128.d0/(x*x) - 5.d0)/3.d0
ampl = (.75D0 + csevl (z, bm0cs, nbm0))/SQRT(x)
theta = x - pi4 + csevl (z, bt02cs, nbt02)/x
RETURN

20   IF (x > xmax) CALL lnkerr ('d9b0mp:no precision x is big')

z = 128.d0/(x*x) - 1.d0
ampl = (.75D0 + csevl (z, bm02cs, nbm02))/SQRT(x)
theta = x - pi4 + csevl (z, bth0cs, nbth0)/x
RETURN

END SUBROUTINE d9b0mp



!deck d9b1mp

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-21  Time: 16:00:12

SUBROUTINE d9b1mp (x, ampl, theta)
!***begin prologue  d9b1mp
!***subsidiary
!***purpose  evaluate the modulus and phase for the j1 and y1 bessel
!            functions.
!***library   slatec (fnlib)
!***category  c10a1
!***type      double precision (d9b1mp-d)
!***keywords  bessel function, fnlib, modulus, phase, special functions
!***author  fullerton, w., (lanl)
!***description

! evaluate the modulus and phase for the bessel j1 and y1 functions.

! series for bm1        on the interval  1.56250e-02 to  6.25000e-02
!                                        with weighted error   4.91e-32
!                                         log weighted error  31.31
!                               significant figures required  30.04
!                                    decimal places required  32.09

! series for bt12       on the interval  1.56250e-02 to  6.25000e-02
!                                        with weighted error   3.33e-32
!                                         log weighted error  31.48
!                               significant figures required  31.05
!                                    decimal places required  32.27

! series for bm12       on the interval  0.          to  1.56250e-02
!                                        with weighted error   5.01e-32
!                                         log weighted error  31.30
!                               significant figures required  29.99
!                                    decimal places required  32.10

! series for bth1       on the interval  0.          to  1.56250e-02
!                                        with weighted error   2.82e-32
!                                         log weighted error  31.55
!                               significant figures required  31.12
!                                    decimal places required  32.37

!***see also  dbesj1, dbesy1
!***references  (none)
!***routines called  d1mach, dcsevl, initds, xermsg
!***revision history  (yymmdd)
!   770701  date written
!   890531  changed all specific intrinsics to generic.  (wrb)
!   890531  revision date from version 3.2
!   891214  prologue converted to version 4.0 format.  (bab)
!   900315  calls to xerror changed to calls to xermsg.  (thj)
!   900720  routine changed from user-callable to subsidiary.  (wrb)
!   920618  removed space from variable name and code restructured to
!           use if-then-else.  (rwc, wrb)
!***end prologue  d9b1mp

REAL*8, INTENT(IN)                       :: x
REAL*8, INTENT(OUT)                      :: ampl
REAL*8, INTENT(OUT)                      :: theta
REAL*8  bm1cs(37), bt12cs(39),  &
    bm12cs(40), bth1cs(44), xmax, pi4, z, r1mach, csevl
LOGICAL :: first
SAVE bm1cs, bt12cs, bth1cs, bm12cs, pi4, nbm1, nbt12,  &
    nbm12, nbth1, xmax, first
DATA bm1cs(  1) / +.1069845452618063014969985308538d+0      /
DATA bm1cs(  2) / +.3274915039715964900729055143445d-2      /
DATA bm1cs(  3) / -.2987783266831698592030445777938d-4      /
DATA bm1cs(  4) / +.8331237177991974531393222669023d-6      /
DATA bm1cs(  5) / -.4112665690302007304896381725498d-7      /
DATA bm1cs(  6) / +.2855344228789215220719757663161d-8      /
DATA bm1cs(  7) / -.2485408305415623878060026596055d-9      /
DATA bm1cs(  8) / +.2543393338072582442742484397174d-10     /
DATA bm1cs(  9) / -.2941045772822967523489750827909d-11     /
DATA bm1cs( 10) / +.3743392025493903309265056153626d-12     /
DATA bm1cs( 11) / -.5149118293821167218720548243527d-13     /
DATA bm1cs( 12) / +.7552535949865143908034040764199d-14     /
DATA bm1cs( 13) / -.1169409706828846444166290622464d-14     /
DATA bm1cs( 14) / +.1896562449434791571721824605060d-15     /
DATA bm1cs( 15) / -.3201955368693286420664775316394d-16     /
DATA bm1cs( 16) / +.5599548399316204114484169905493d-17     /
DATA bm1cs( 17) / -.1010215894730432443119390444544d-17     /
DATA bm1cs( 18) / +.1873844985727562983302042719573d-18     /
DATA bm1cs( 19) / -.3563537470328580219274301439999d-19     /
DATA bm1cs( 20) / +.6931283819971238330422763519999d-20     /
DATA bm1cs( 21) / -.1376059453406500152251408930133d-20     /
DATA bm1cs( 22) / +.2783430784107080220599779327999d-21     /
DATA bm1cs( 23) / -.5727595364320561689348669439999d-22     /
DATA bm1cs( 24) / +.1197361445918892672535756799999d-22     /
DATA bm1cs( 25) / -.2539928509891871976641440426666d-23     /
DATA bm1cs( 26) / +.5461378289657295973069619199999d-24     /
DATA bm1cs( 27) / -.1189211341773320288986289493333d-24     /
DATA bm1cs( 28) / +.2620150977340081594957824000000d-25     /
DATA bm1cs( 29) / -.5836810774255685901920938666666d-26     /
DATA bm1cs( 30) / +.1313743500080595773423615999999d-26     /
DATA bm1cs( 31) / -.2985814622510380355332778666666d-27     /
DATA bm1cs( 32) / +.6848390471334604937625599999999d-28     /
DATA bm1cs( 33) / -.1584401568222476721192960000000d-28     /
DATA bm1cs( 34) / +.3695641006570938054301013333333d-29     /
DATA bm1cs( 35) / -.8687115921144668243012266666666d-30     /
DATA bm1cs( 36) / +.2057080846158763462929066666666d-30     /
DATA bm1cs( 37) / -.4905225761116225518523733333333d-31     /
DATA bt12cs(  1) / +.73823860128742974662620839792764d+0     /
DATA bt12cs(  2) / -.33361113174483906384470147681189d-2     /
DATA bt12cs(  3) / +.61463454888046964698514899420186d-4     /
DATA bt12cs(  4) / -.24024585161602374264977635469568d-5     /
DATA bt12cs(  5) / +.14663555577509746153210591997204d-6     /
DATA bt12cs(  6) / -.11841917305589180567005147504983d-7     /
DATA bt12cs(  7) / +.11574198963919197052125466303055d-8     /
DATA bt12cs(  8) / -.13001161129439187449366007794571d-9     /
DATA bt12cs(  9) / +.16245391141361731937742166273667d-10    /
DATA bt12cs( 10) / -.22089636821403188752155441770128d-11    /
DATA bt12cs( 11) / +.32180304258553177090474358653778d-12    /
DATA bt12cs( 12) / -.49653147932768480785552021135381d-13    /
DATA bt12cs( 13) / +.80438900432847825985558882639317d-14    /
DATA bt12cs( 14) / -.13589121310161291384694712682282d-14    /
DATA bt12cs( 15) / +.23810504397147214869676529605973d-15    /
DATA bt12cs( 16) / -.43081466363849106724471241420799d-16    /
DATA bt12cs( 17) / +.80202544032771002434993512550400d-17    /
DATA bt12cs( 18) / -.15316310642462311864230027468799d-17    /
DATA bt12cs( 19) / +.29928606352715568924073040554666d-18    /
DATA bt12cs( 20) / -.59709964658085443393815636650666d-19    /
DATA bt12cs( 21) / +.12140289669415185024160852650666d-19    /
DATA bt12cs( 22) / -.25115114696612948901006977706666d-20    /
DATA bt12cs( 23) / +.52790567170328744850738380799999d-21    /
DATA bt12cs( 24) / -.11260509227550498324361161386666d-21    /
DATA bt12cs( 25) / +.24348277359576326659663462400000d-22    /
DATA bt12cs( 26) / -.53317261236931800130038442666666d-23    /
DATA bt12cs( 27) / +.11813615059707121039205990399999d-23    /
DATA bt12cs( 28) / -.26465368283353523514856789333333d-24    /
DATA bt12cs( 29) / +.59903394041361503945577813333333d-25    /
DATA bt12cs( 30) / -.13690854630829503109136383999999d-25    /
DATA bt12cs( 31) / +.31576790154380228326413653333333d-26    /
DATA bt12cs( 32) / -.73457915082084356491400533333333d-27    /
DATA bt12cs( 33) / +.17228081480722747930705920000000d-27    /
DATA bt12cs( 34) / -.40716907961286507941068800000000d-28    /
DATA bt12cs( 35) / +.96934745136779622700373333333333d-29    /
DATA bt12cs( 36) / -.23237636337765716765354666666666d-29    /
DATA bt12cs( 37) / +.56074510673522029406890666666666d-30    /
DATA bt12cs( 38) / -.13616465391539005860522666666666d-30    /
DATA bt12cs( 39) / +.33263109233894654388906666666666d-31    /
DATA bm12cs(  1) / +.9807979156233050027272093546937d-1      /
DATA bm12cs(  2) / +.1150961189504685306175483484602d-2      /
DATA bm12cs(  3) / -.4312482164338205409889358097732d-5      /
DATA bm12cs(  4) / +.5951839610088816307813029801832d-7      /
DATA bm12cs(  5) / -.1704844019826909857400701586478d-8      /
DATA bm12cs(  6) / +.7798265413611109508658173827401d-10     /
DATA bm12cs(  7) / -.4958986126766415809491754951865d-11     /
DATA bm12cs(  8) / +.4038432416421141516838202265144d-12     /
DATA bm12cs(  9) / -.3993046163725175445765483846645d-13     /
DATA bm12cs( 10) / +.4619886183118966494313342432775d-14     /
DATA bm12cs( 11) / -.6089208019095383301345472619333d-15     /
DATA bm12cs( 12) / +.8960930916433876482157048041249d-16     /
DATA bm12cs( 13) / -.1449629423942023122916518918925d-16     /
DATA bm12cs( 14) / +.2546463158537776056165149648068d-17     /
DATA bm12cs( 15) / -.4809472874647836444259263718620d-18     /
DATA bm12cs( 16) / +.9687684668292599049087275839124d-19     /
DATA bm12cs( 17) / -.2067213372277966023245038117551d-19     /
DATA bm12cs( 18) / +.4646651559150384731802767809590d-20     /
DATA bm12cs( 19) / -.1094966128848334138241351328339d-20     /
DATA bm12cs( 20) / +.2693892797288682860905707612785d-21     /
DATA bm12cs( 21) / -.6894992910930374477818970026857d-22     /
DATA bm12cs( 22) / +.1830268262752062909890668554740d-22     /
DATA bm12cs( 23) / -.5025064246351916428156113553224d-23     /
DATA bm12cs( 24) / +.1423545194454806039631693634194d-23     /
DATA bm12cs( 25) / -.4152191203616450388068886769801d-24     /
DATA bm12cs( 26) / +.1244609201503979325882330076547d-24     /
DATA bm12cs( 27) / -.3827336370569304299431918661286d-25     /
DATA bm12cs( 28) / +.1205591357815617535374723981835d-25     /
DATA bm12cs( 29) / -.3884536246376488076431859361124d-26     /
DATA bm12cs( 30) / +.1278689528720409721904895283461d-26     /
DATA bm12cs( 31) / -.4295146689447946272061936915912d-27     /
DATA bm12cs( 32) / +.1470689117829070886456802707983d-27     /
DATA bm12cs( 33) / -.5128315665106073128180374017796d-28     /
DATA bm12cs( 34) / +.1819509585471169385481437373286d-28     /
DATA bm12cs( 35) / -.6563031314841980867618635050373d-29     /
DATA bm12cs( 36) / +.2404898976919960653198914875834d-29     /
DATA bm12cs( 37) / -.8945966744690612473234958242979d-30     /
DATA bm12cs( 38) / +.3376085160657231026637148978240d-30     /
DATA bm12cs( 39) / -.1291791454620656360913099916966d-30     /
DATA bm12cs( 40) / +.5008634462958810520684951501254d-31     /
DATA bth1cs(  1) / +.74749957203587276055443483969695d+0     /
DATA bth1cs(  2) / -.12400777144651711252545777541384d-2     /
DATA bth1cs(  3) / +.99252442404424527376641497689592d-5     /
DATA bth1cs(  4) / -.20303690737159711052419375375608d-6     /
DATA bth1cs(  5) / +.75359617705690885712184017583629d-8     /
DATA bth1cs(  6) / -.41661612715343550107630023856228d-9     /
DATA bth1cs(  7) / +.30701618070834890481245102091216d-10    /
DATA bth1cs(  8) / -.28178499637605213992324008883924d-11    /
DATA bth1cs(  9) / +.30790696739040295476028146821647d-12    /
DATA bth1cs( 10) / -.38803300262803434112787347554781d-13    /
DATA bth1cs( 11) / +.55096039608630904934561726208562d-14    /
DATA bth1cs( 12) / -.86590060768383779940103398953994d-15    /
DATA bth1cs( 13) / +.14856049141536749003423689060683d-15    /
DATA bth1cs( 14) / -.27519529815904085805371212125009d-16    /
DATA bth1cs( 15) / +.54550796090481089625036223640923d-17    /
DATA bth1cs( 16) / -.11486534501983642749543631027177d-17    /
DATA bth1cs( 17) / +.25535213377973900223199052533522d-18    /
DATA bth1cs( 18) / -.59621490197413450395768287907849d-19    /
DATA bth1cs( 19) / +.14556622902372718620288302005833d-19    /
DATA bth1cs( 20) / -.37022185422450538201579776019593d-20    /
DATA bth1cs( 21) / +.97763074125345357664168434517924d-21    /
DATA bth1cs( 22) / -.26726821639668488468723775393052d-21    /
DATA bth1cs( 23) / +.75453300384983271794038190655764d-22    /
DATA bth1cs( 24) / -.21947899919802744897892383371647d-22    /
DATA bth1cs( 25) / +.65648394623955262178906999817493d-23    /
DATA bth1cs( 26) / -.20155604298370207570784076869519d-23    /
DATA bth1cs( 27) / +.63417768556776143492144667185670d-24    /
DATA bth1cs( 28) / -.20419277885337895634813769955591d-24    /
DATA bth1cs( 29) / +.67191464220720567486658980018551d-25    /
DATA bth1cs( 30) / -.22569079110207573595709003687336d-25    /
DATA bth1cs( 31) / +.77297719892989706370926959871929d-26    /
DATA bth1cs( 32) / -.26967444512294640913211424080920d-26    /
DATA bth1cs( 33) / +.95749344518502698072295521933627d-27    /
DATA bth1cs( 34) / -.34569168448890113000175680827627d-27    /
DATA bth1cs( 35) / +.12681234817398436504211986238374d-27    /
DATA bth1cs( 36) / -.47232536630722639860464993713445d-28    /
DATA bth1cs( 37) / +.17850008478186376177858619796417d-28    /
DATA bth1cs( 38) / -.68404361004510395406215223566746d-29    /
DATA bth1cs( 39) / +.26566028671720419358293422672212d-29    /
DATA bth1cs( 40) / -.10450402527914452917714161484670d-29    /
DATA bth1cs( 41) / +.41618290825377144306861917197064d-30    /
DATA bth1cs( 42) / -.16771639203643714856501347882887d-30    /
DATA bth1cs( 43) / +.68361997776664389173535928028528d-31    /
DATA bth1cs( 44) / -.28172247861233641166739574622810d-31    /
DATA pi4 / 0.785398163397448309615660845819876d0 /
DATA first /.true./
!***first executable statement  d9b1mp
IF (first) THEN
  eta = 0.1*REAL(r1mach(3))
  nbm1 = initds (bm1cs, 37, eta)
  nbt12 = initds (bt12cs, 39, eta)
  nbm12 = initds (bm12cs, 40, eta)
  nbth1 = initds (bth1cs, 44, eta)
  
  xmax = 1.0D0/r1mach(4)
END IF
first = .false.

IF (x < 4.0D0) THEN
  CALL lnkerr ('d9b1mp:x must be .ge. 4')
  ampl = 0.0D0
  theta = 0.0D0
ELSE IF (x <= 8.0D0) THEN
  z = (128.0D0/(x*x) - 5.0D0)/3.0D0
  ampl = (0.75D0 + csevl (z, bm1cs, nbm1))/SQRT(x)
  theta = x - 3.0D0*pi4 + csevl (z, bt12cs, nbt12)/x
ELSE
  IF (x > xmax) CALL lnkerr ('d9b1mp:no precision x too big')
  
  z = 128.0D0/(x*x) - 1.0D0
  ampl = (0.75D0 + csevl (z, bm12cs, nbm12))/SQRT(x)
  theta = x - 3.0D0*pi4 + csevl (z, bth1cs, nbth1)/x
END IF
RETURN
END SUBROUTINE d9b1mp



!deck envj

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-21  Time: 11:41:51

FUNCTION envj(n,x)


INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN)                       :: x
REAL*8  envj
envj=0.5D0*LOG10(6.28D0*n)-n*LOG10(1.36D0*x/n)
RETURN
END FUNCTION envj


!deck @(#)factrl.f 1.1  11/30/90

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-20  Time: 15:29:54

FUNCTION factrl(n)
!***begin prologue     factrl
!***date written       901115  (yymmdd)
!***revision date      yymmdd  (yymmdd)
!***keywords           factorial
!***author             martin, richard (lanl)
!***source             @(#)factrl.f 1.1   11/30/90
!***purpose            returns the factorial
!***description
!                      x=factrl(n)
!                        x= n * (n-1) * (n-2) ... *1

!***references
!***routines called    (none)
!***end prologue       factrl

IMPLICIT INTEGER(a-z)
INTEGER, INTENT(IN)                      :: n
REAL*8 factrl
REAL*8, PARAMETER :: one=1.0D+00


IF(n < 0) THEN
  CALL lnkerr('bad arguments to factorial function')
ELSE IF (n == 0) THEN
  factrl=one
ELSE
  factrl=one
  DO  i=n,2,-1
    factrl=factrl*(FLOAT(n))
  END DO
END IF


RETURN
END FUNCTION factrl


!deck initds

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-21  Time: 16:00:33

FUNCTION initds (os, nos, eta)
!***begin prologue  initds
!***purpose  determine the number of terms needed in an orthogonal
!            polynomial series so that it meets a specified accuracy.
!***library   slatec (fnlib)
!***category  c3a2
!***type      double precision (inits-s, initds-d)
!***keywords  chebyshev, fnlib, initialize, orthogonal polynomial,
!             orthogonal series, special functions
!***author  fullerton, w., (lanl)
!***description

!  initialize the orthogonal series, represented by the array os, so
!  that initds is the number of terms needed to insure the error is no
!  larger than eta.  ordinarily, eta will be chosen to be one-tenth
!  machine precision.

!             input arguments --
!   os     double precision array of nos coefficients in an orthogonal
!          series.
!   nos    number of coefficients in os.
!   eta    single precision scalar containing requested accuracy of
!          series.

!***references  (none)
!***routines called  xermsg
!***revision history  (yymmdd)
!   770601  date written
!   890531  changed all specific intrinsics to generic.  (wrb)
!   890831  modified array declarations.  (wrb)
!   891115  modified error message.  (wrb)
!   891115  revision date from version 3.2
!   891214  prologue converted to version 4.0 format.  (bab)
!   900315  calls to xerror changed to calls to xermsg.  (thj)
!***end prologue  initds


REAL*8, INTENT(IN OUT)                   :: os(*)
INTEGER, INTENT(IN)                      :: nos
REAL, INTENT(IN)                         :: eta

!***first executable statement  initds
IF (nos < 1) CALL lnkerr ('initds:no. coef. less than 1')

ERR = 0.
DO  ii = 1,nos
  i = nos + 1 - ii
  ERR = ERR + ABS(REAL(os(i)))
  IF (ERR > eta) GO TO 20
END DO

20 IF (i == nos) CALL lnkerr ('initds:cheby series too short')
initds = i

RETURN
END FUNCTION initds


!deck intrv.f

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-20  Time: 15:30:37

!***begin prologue     intrv
!***date written       930623   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           newton-raphson root finder
!***author             schneider, barry (nsf)
!***source             math
!***purpose            locate interval in which a function goes through zero
!***
!***                   Function is given in fun(n), zval(m) returns location 
!***                   of sign changes, cntzro returns number of sign changes
!***                   in interval starting at rl, ending at rl+del*n
!***                   m is maximum number of sign changes found 

!***references

!***routines called
!***end prologue       intrv

SUBROUTINE intrv(fun,zval,rl,del,n,m,cntzro,prnt)

IMPLICIT INTEGER (a-z)
REAL*8, INTENT(IN)                       :: fun(n)
REAL*8, INTENT(OUT)                      :: zval(m+1)
REAL*8, INTENT(IN)                       :: rl
REAL*8, INTENT(IN)                       :: del
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(OUT)                     :: cntzro
LOGICAL, INTENT(IN)                      :: prnt
REAL*8  val, x, prd
LOGICAL                                  :: zero  

COMMON/io/inp, iout


cntzro=0
k=1
zero=.true.
DO while(zero)
   x=rl+(k-1)*del
   zval(1)=x
   val=fun(k)
   k=k+1
   IF(abs(val).gt.0) then
      zero=.false.
   END IF
END DO

DO  i=k,n
  IF (cntzro == m) EXIT
  x=rl+(i-1)*del
  prd=val*fun(i)
  IF (prd < 0.d0) THEN
    cntzro=cntzro+1
    zval(cntzro+1)=x
    val=fun(i)
  END IF
END DO
!100 WRITE(iout,1) cntzro
!1 FORMAT(/,5X,'there are ',i3,' zeros in the input interval')
IF (prnt) THEN
  WRITE(iout,2) ( zval(i),i=1,cntzro+1)
END IF
2 FORMAT(/,5X,'the intervals for the zeros',(/,5(1X,e15.8)))
RETURN
END SUBROUTINE intrv


!deck msta1

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-21  Time: 11:33:38

FUNCTION msta1(x,mp)

!       ===================================================
!       purpose: determine the starting point for backward
!                recurrence such that the magnitude of
!                jn(x) at that point is about 10^(-mp)
!       input :  x     --- argument of jn(x)
!                mp    --- value of magnitude
!       output:  msta1 --- starting point
!       ===================================================


IMPLICIT INTEGER(a-z)
REAL*8, INTENT(IN OUT)                   :: x
INTEGER, INTENT(IN)                      :: mp
REAL*8  a0, envj, f0, f1, f

a0=ABS(x)
n0=INT(1.1*a0)+1
f0=envj(n0,a0)-mp
n1=n0+5
f1=envj(n1,a0)-mp
DO  it=1,20
  nn=n1-(n1-n0)/(1.0D0-f0/f1)
  f=envj(nn,a0)-mp
  IF(ABS(nn-n1) < 1) GO TO 20
  n0=n1
  f0=f1
  n1=nn
  f1=f
END DO
20     msta1=nn
RETURN
END FUNCTION msta1



!deck msta2

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-21  Time: 11:33:43

FUNCTION msta2(x,n,mp)

!       ===================================================
!       purpose: determine the starting point for backward
!                recurrence such that all jn(x) has mp
!                significant digits
!       input :  x  --- argument of jn(x)
!                n  --- order of jn(x)
!                mp --- significant digit
!       output:  msta2 --- starting point
!       ===================================================


IMPLICIT INTEGER (a-z)
REAL*8, INTENT(IN OUT)                   :: x
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: mp
REAL*8  a0, hmp, envj, ejn, obj, f0, f1, f

a0=DABS(x)
hmp=0.5D0*mp
ejn=envj(n,a0)
IF (ejn <= hmp) THEN
  obj=mp
  n0=INT(1.1*a0)
ELSE
  obj=hmp+ejn
  n0=n
END IF
f0=envj(n0,a0)-obj
n1=n0+5
f1=envj(n1,a0)-obj
DO  it=1,20
  nn=n1-(n1-n0)/(1.0D0-f0/f1)
  f=envj(nn,a0)-obj
  IF (ABS(nn-n1) < 1) GO TO 20
  n0=n1
  f0=f1
  n1=nn
  f1=f
END DO
20      msta2=nn+10
RETURN
END FUNCTION msta2


