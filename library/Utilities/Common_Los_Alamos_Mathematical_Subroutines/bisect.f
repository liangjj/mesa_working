*deck bisect
      subroutine bisect(n,eps1,d,e,e2,lb,ub,mm,m,w,ind,ierr,rv4,rv5)
! NOTE  REQUIRES FOLLOWING STURM-SEQUENCE SUBROUTINE THTW
! eigenvalues of symmetric tridiagonal matrix
! inputs:
! n matrix order  d(n) diagonals
! e(n) subdiagonals, e(1) arbitrary   e2(n) squares
! eps1 abs err tol for eigenvals, if <0 reset to - machp*norm
! lb, ub  lower and upper endpoints of interval to be searched
! mm upper bound on number of evals in (lb,ub)
! outputs:
! m number of evals in (lb,ub)
! w(mm) contains m evals in ascending order
! ind(mm) submatrix indices of evals
! ierr nonzero = 3*n+1 if m>mm, otherwise 0
! rv4(n) ,rv5(n) used to hold lower and upper bounds
      implicit none
      integer i,j,k,l,m,n,p,q,r,s,mm,m1,m2,tag,ierr
      real*8 d(n),e(n),e2(n),w(mm),rv4(n),rv5(n)
      real*8 u,v,lb,t1,t2,ub,xu,x0,x1,eps1,machep
      integer ind(mm)
	external thtw
        data machep/1.0d0/
	if (machep .eq. 1.0d0) then
	   do while (1.0d0 + machep .gt. 1.0d0)
	      machep = 0.5d0*machep
	   end do ! while
           machep = 2.0d0*machep
	end if  ! machep
	ierr = 0
	tag = 0
	t1 = lb
	t2 = ub
	e2(1) = 0.0d0
	do i = 2 ,n
	   if ( abs(e(i)) .le. machep * ( abs(d(i)) + abs(d(i-1))))
	1		e2(i) = 0.0d0
	end do ! i
	p = 1
	q = n
	call thtw(p,q,e,e2,d,ub,m)	! sturm seq at ub
	call thtw(p,q,e,e2,d,lb,s)	! sturm seq at lb
	m = m - s			! # evals in (lb,ub)
	if (m .gt. mm) then	! error exit, underestimate of number of ev
	   ierr = 3 * n + 1
	   lb = t1
	   ub = t2
	   return
	end if ! m
	q = 0
	r = 0
! establish and process submatrices, refining interval 
!					by gerschgorin bounds
	do while ( q .lt. n .and. r .ne. m )
	   tag = tag + 1
	   p = q + 1
	   xu = d(p)
	   x0 = d(p)
	   u = 0.0d0
	   v = 1.0d0
	   q = p - 1
	   do while ( v .ne. 0.0d0 .and. q .lt. n )
	      q = q + 1
              x1 = u
              u = 0.0d0
              v = 0.0d0
	      if ( q .ne. n ) then
                 u = abs(e(q+1))
                 v = e2(q+1)
	      end if ! q
              xu = min(d(q)-(x1+u),xu)
              x0 = max(d(q)+(x1+u),x0)
	   end do ! while v and q
           if (eps1 .le. 0.0d0) eps1 = - max(abs(xu),abs(x0)) * machep
           if (p .eq. q) then  ! check for isolated root in interval
              if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 1000
              m1 = p
              m2 = p
              rv5(p) = d(p)
	   else
              x1 = x1 * float(q-p+1)
              lb = max(t1,xu-x1)	! Gerschgorin refinement
              ub = min(t2,x0+x1)
	      call thtw(p,q,e,e2,d,lb,s)
              m1 = s + 1
	      call thtw(p,q,e,e2,d,ub,m2)
              if (m1 .gt. m2) go to 1000
! find roots by bisection
              x0 = ub
              do i = m1, m2
                 rv5(i) = ub
                 rv4(i) = lb
              end do ! i
! loop for k-th eigenvalue, starting at m2
              k = m2
	      do while ( k .ge. m1 )
	         xu = lb
	         i = k
	         do while ( i .ge. m1 .and. xu .ge. rv4(i) )
	            i = i - 1
	         end do ! while i
	         if ( i .ge. m1 ) xu = rv4(i)
		 x0 = min( x0 ,rv5(k) )
	         x1 = (xu + x0) * 0.5d0
		 do while ( x0-xu .gt. 
	1		abs(eps1) + 2*machep*(abs(xu)+abs(x0)))
	            call thtw(p,q,e,e2,d,x1,s)
! refine intervals
	            if (s .lt. k) then
                       xu = x1
                       if (s .lt. m1) then
                          rv4(m1) = x1
	               else		! s ge m1
                          rv4(s+1) = x1
                          rv5(s) = min( x1 ,rv5(s) )
	               end if ! s
	            else		! s ge k
                       x0 = x1
	            end if ! s
                    x1 = (xu + x0) * 0.5d0
	         end do	! while s2  ,   k-th eigenvalue found
	         rv5(k) = x1
                 k = k - 1
	      end do ! while k
	   end if ! p
! order eigenvalues and submatrix associations
           s = r
           r = r + m2 - m1 + 1
           j = 1
           k = m1
	   do l = 1 ,r
              if (j .le. s ) then
                 if (k .gt. m2) then
		    go to 1000
                 else
                    if (rv5(k) .lt. w(l)) then
	               do i = l + s - j ,l ,-1
                          w(i+1) = w(i)
                          ind(i+1) = ind(i)
                       end do
                       w(l) = rv5(k)
                       ind(l) = tag
                       k = k + 1
		    else
		       j = j + 1
		    end if ! rv5
	         end if ! k
	      else
                 w(l) = rv5(k)
                 ind(l) = tag
                 k = k + 1
	      end if ! j
	   end do ! q
1000	   continue
	end do ! while q and r
	lb = t1
	ub = t2
	return
        end
