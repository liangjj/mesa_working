!=======================================================================
  FUNCTION a(i,j,k)
!=======================================================================
!
!  determine the coefficient in the potential for electron i of
!  y^k(j,j)
!
!----------------------------------------------------------------------
!
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i,j,k
     REAL(KIND=8) :: a

      if (i > nclosd .and. j > nclosd) then
         istart = ijptr(i-nclosd,j-nclosd) + 1
         a = coef(istart + k/2)
      else if (i == j) then
         c = sum(i) - 1.d0
         if (k == 0) then
            a = c
         else
            a = -c*ca(l(i),k)
         end if
      else if (k.eq.0) then
         a = sum(j)
      else
         a = 0.d0
      end if
  END FUNCTION a
!=======================================================================
  SUBROUTINE add(c,k,i,j,first) 
!=======================================================================
!   add a slater integral to the data structure associated with the
!   energy expression
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k,i,j
    LOGICAL, INTENT(IN) :: first
    REAL(KIND=8), INTENT(IN) :: c

    ip = ijptr(i-nclosd,j-nclosd)
   
    if (first) then
       coef(ip+k/2+1) = c/sum(i) + coef(ip+k/2+1)
    else
       ip = ip + min(l(i),l(j)) +1 + (k-abs(l(i)-l(j)))/2 + 1
       coef(ip) = coef(ip) + c/sum(i)
    end if

  END FUNCTION add
!=======================================================================
  FUNCTION b(i,j,k) 
!=======================================================================
!
!   determine the coefficient of the y^k(i,j)p(j) term in the exchange
!   expression of electron i
!
!----------------------------------------------------------------------
!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,k
   REAL(KIND=8) :: b

 
   if (i == j) then
       b = 0.d0
   else if (i > nclosd .and. j > nclosd) then
 
    !  ll is the number of direct terms
    !  istart the beginning of the exchange terms
 
    ll = min(l(i),l(j)) + 1
    istart = ijptr(i-nclosd,j-nclosd) + 1 + ll
    kk = (k - abs(l(i)-l(j)))/2
    b = coef(istart + kk)
  else
    b = -sum(j)*cb(l(i),l(j),k)
  end if
  END FUNCTION b
!=======================================================================
        subroutine bxvpw(c,b,v,y,w)
!=======================================================================
!   Computes y = c* b * v + y   where b is a symmetric, banded matrix
!   and v, w are vectors
!
!   Written by C. F. Fischer
!------------------------------------------------------------------------
!
!
!   on entry
!   --------
!     nt      the leading dimension of arrays.
!     k       the number of diagonals
!     n       the order of the matrix
!     c,ic    coefficients
!     b       the symmetric, banded matrix in column storage mode
!     v       vector
!     w       working array
!
!   on exit
!   -------
!     y       y = c*B*v +y
!-----------------------------------------------------------------------
!
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i,j,k
     REAL(KIND=8), INTENT(IN), DIMENSION(nt,:) :: b
     REAL(KIND=8), INTENT(IN), DIMENSION(:) :: v,y
     REAL(KIND=8), INTENT(INOUT), DIMENSION(:) :: y

!
! ...   initialize the w array
!
	 w = 0.d0
!
!      .. contribution from sub-diagonals
!
       do i=1,k
         do j=k-i+1,n
	   w(j) = w(j)  +b(j,i)*v(j-k+i)
         end do
       end do
!
!      .. contribution from super-diagonals
!
       do i=1,k-1
         do j=1,n-k+i
	   w(j) = w(j) + b(j+k-i,i)*v(j+k-i)
         end do
       end do
!
       if ( c .ne. 1.d0) then
           y = y + c*w
       else
           y = y + w
       end if
    END SUBROUTINE bxvpw
      
!=======================================================================
   FUNCTION ca(l,k) 
!=======================================================================
! Computes the direct average interaction
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, k
    REAL(KIND=8) :: ca

      ca = rme(l,l,k)**2

    END FUNCTION ca

!=======================================================================
   FUNCTION cb(l,lp,k) 
!=======================================================================
!  Compute the average exchange interaction
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN): l,lp,k
    REAL(KIND=8):: cb

      cb = rme(l,lp,k)**2/(2*(2*l+1)*(2*lp+1))

    END FUNCTION cb
!=======================================================================
   FUNCTION ekin(i,ii,rel) 
!=======================================================================
!
!       returns the value of the integral of
!
!         (2/r)p (y p  + x )
!               j  i i    i
!
!   integrated with respect to r.
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,ii
    LOGICAL, INTENT(IN) :: rel
    REAL(KIND=8) :: ekin

    REAL(KIND=8), DIMENSION(ns,nt) :: yrm, xv, w

      call xch(i,rel,xv)
      call potl(i,rel,yrm)
      call bxvpw(nt,ks,ns,2.0d0,yrm,p(1,i),xv,w)
      ekin = SUM(p(:,ii)*xv)
      print *, 'Ekin', i,ii,ekin
   END FUNCTION ekin
!=======================================================================
  SUBROUTINE eptr(el, elsymb, iel) 
!=======================================================================
!
!   Determines the position of the electron in the electron list
!   Zero if not found.
!      
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    CHARACTER(LEN=*), DIMENSION(*), INTENT(IN) :: el
    CHARACTER(LEN=*), INTENT(IN) :: elsymb
    INTEGER, INTENT(OUT) :: iel

    INTEGER :: n, i
    n = size(el)
    iel = 0
    do  i=1,n
      if (el(i) .eq. elsymb ) then
        iel = i
        return
      endif
    end do
    end
!=======================================================================
   SUBROUTINE factrl(nfact)
!=======================================================================
!  Computes the log of the factorial function
!      gam(i) = log( gamma(i-1) ), where gamma(i) = factorial i-1
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nfact
    REAL(KIND=8), DIMENSION(*) :: gamma

      gamma = 1
      gam(1) = 0
      do  i=1,nfact-1
         gamma=i*gamma
         gam(i+1) = dlog(gamma)
      end do
      do i = nfact+1,(100)
         x = i-1
         gam(i) = gam(i-1) + dlog(x)
      end do
      
   END SUBROUTINE factrl
!=======================================================================
  FUNCTION grad(i,j,ier) 
!=======================================================================
!  The grad function computes the following directly
!          <P(J)^D + L(I)/R ^P(I)> WITH L(I) > L(J)
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j
    INTEGER, INTENT(OUT) :: ier
    REAL(KIND=8) :: grad

    LOGICAL :: ixch
    REAL(KIND=8), DIMENSION(ns,ks) :: g, a
      
    IF ( IABS(L(I) - L(J)) .NE. 1) THEN
      ier = 1
      grad = 0
    else
      LL = MAX0(L(I),L(J))
      i1 = i 
      i2 = j
      ixch = .false.
      IF ( L(I) < L(J) ) THEN
	i1 = j
	i2 = i
	ixch = .true.
      end if

      !  form the sums

      last = min0(max(i1),max(i2))
      grad = SUM( p(2:last,i1)*p(2:last,i2)*rm1(2:last,ks)
      grad = ll*grad
      do jp = 2,ks
        ! .. j' < jj
	lastjp = last-jp+1
        grad = grad  &       
	    +  SUM( g(2:lastjp,ks-jp+1)*                                     &
                          (p(2:lastjp,i2)*p(2:lastjp-jp+1,i1) -              &
	                   p(2:lastjp,i1)*p(2:lastjp-jp+1,i2)) )             &
	    + ll*SUM( rm1(2:lastjp,ks-jp+1)*                                 &
			  (p(2:lastjp,i2)*p(2:lastjp-jp+1,i1)+               &
                           p(2:lastjp,i1)*p(2:lastjp-jp+1,i2)) )
      end do
      IF ( ixch ) grad = -grad
    end if
    END FUNCTION grad
      integer function idmax(n,dx,incx)
!=======================================================================
  FUNCTION idmax(n,dx,incx) 
!=======================================================================
!     finds the index of element having largest value.
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) ::n, incx
    REAL(KIND=8), INTENT(IN) :: dx

    INTEGER :: idmax, i

    if( n < 1 ) then
      idmax = 0
    else if (n == 1) then
      idmax = 1
    else
      if (incx > 1) then
        ! code for increment not equal to 1
        ix = 1
        dmax = dx(1)
        ix = ix + incx
        do i = 2,n
         if(dx(ix) > dmax) then
           idmax = i
           dmax = dx(ix)
          end if
       	end do
      else
        !code for increment equal to 1
	idmax = 1
        dmax = dx(1)
        do i = 2,n
         if(dx(i) > dmax) then
           idmax = i
           dmax = dx(i)
	 end if
	end do
      end if
   END FUNCTION idmax
!=======================================================================
  FUNCTION lval(symbol) 
!=======================================================================
!   Convert symbol to l-value
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER :: lval
    CHARACTER, INTENT(IN) :: symbol
    CHARACTER*22 :: set = 'spdfghiklmnspdfghiklmn'

      locate = index(set,symbol)
      if ( locate .le. 11) then
            lval = locate - 1
         else
            lval = locate - 12
      endif

   END FUNCITON lval
!=======================================================================
  FUNCTION quad(a,b) 
!=======================================================================
!       Returns the value of <a, S b> where a and b are spline
!     expansion coefficients
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(ns):: a, b

      quad = 0
      do i = 2,ns
	quad = quad + a(i)*b(i)*sb(i,ks)
      end do
      do m = 2,ks
        ! i-m+1  < m
	do i = m+1,ns
	  quad = quad +(a(i)*b(i-m+1)+a(i-m+1)*b(i))*sb(i,ks-m+1)
        end do
      end do
    END FUNCTION quad
!=======================================================================
  SUBROUTINE reform(str1, str2) 
!=======================================================================
!   Reformat the configuration string from the free-format str1 to the 
!   fixed 8(1x,a3,1x,i4,1x) format
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: str1
    CHARACTER(LEN=80), INTENT(OUT) :: str2

    INTEGER :: i, is, js
!
    
    i = 0
    str2 = '  '
    is = 0
    DO
     js = index(str1(is+1:),'(')
     if (js .ne. 0) then
       if (js .gt. 8) then
         write(*,*) 'ERROR: too many shells for this code'
         STOP
       end if
       i = i+5
       str2(i-js+1:i) = str1(is+1:is+js)
       is = is + js
       js = index(str1(is+1:),')')
       i = i+5
       str2(i-js+1:i) = str1(is+1:is+js)
       is = is + js
     ELSE
       EXIT
     END IF
    END DO
  
  END SUBROUTINE reform

!=======================================================================
  SUBROUTINE reord(el, elc, nwf, ierr) 
!=======================================================================
!
!   Reorder the list of first appearance so that the FUNCTIONs (orbitals)
!   to be varied appear last in the list.
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    CHARACTER(LEN=*), DIMENSION(nwf), INTENT(INOUT) :: el
    CHARACTER(LEN=*), INTENT(IN) :: elc 
*
    call eptr(el, elc, i)
    if (i <> 0) then
      do j = i, nwf-1
        el(j) = el(j+1)
      end do
      el(nwf) = elc
      ierr = 0
    else
      ierr = 1
    end if
  END SUBROUTINE reord


!=======================================================================
  FUNCTION rme(l,lp,k) 
!=======================================================================
!  Evaluates the reduced matrix element (l//c(k)//lp)  -  see fano
!    and racah, irreducible tensorial sets, chap. 14, p. 81
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l,lp,k
    REAL(KIND=8) :: rme, qusqrt
    INTEGER :: i2g, ig, i1,i2,i3

      if (min(l,lp) == 0) then
         rme = 1
       else if ( k == 0) then
         rme = 2*l+1
         rme = dsqrt(rme)
       else
         i2g=l+lp+k
         ig=i2g/2
         if (i2g - 2*ig <> 0) then
             rme = 0
         else
            i1=ig-l
            i2=ig-lp
            i3=ig-k
            qusqrt=(2*l+1)*(2*lp+1)
            rme=dsqrt(qusqrt)*dexp((gam(2*i1+1)+gam(2*i2+1)+gam(2*i3+1)-
     :        gam(i2g+2))/2.d0 +gam(ig+1)-gam(i1+1)-gam(i2+1)-gam(i3+1))
         end if
      end if
   END FUNCTION rme
