
      subroutine ginvit(a,lda,n,c,ipvt,x,ierr)
      implicit real*8(a-h,o-z)
      dimension ipvt(lda),a(lda,*),c(*),x(*)
      data eps/1d-14/
*
*     Generalized inverse Iteration problem
*               Ac = 0
*     where
*             A = A* + k^2*B
*             A is dense matrix
*
*    lda  -  leading dimension of a
*    n    -  size of the matrix
*    c    -  initial guess of the eigenvector on input,
*            final solution on output
*    ipvt -  space for pivot information
*    x    -  working vector
*    ierr -  0 if converged, number of iterations if 30.
*
*    call for LU decomposition of A
*     print *, 'Entering ginvit with n,c =',n,(c(i),i=1,n)
      call dgefa(a,lda,n,ipvt,info)
CTSS  call sgefa(a,lda,n,ipvt,info)
CUNIC call sgefa(a,lda,n,ipvt,info)
      if(info.ne.0)then
	 write(0,*) 'Call to dgefa produced zero in position ',info
CTSS	 write(0,*) 'Call to sgefa produced zero in position ',info
CUNIC    write(0,*) 'Call to sgefa produced zero in position ',info
	 a(info, info) = 1.d-14
      endif

c   Inverse Iteration starts from here

      ierr = 0
      it = 0
      aold = 0
 200  it = it + 1
      call dgesl(a,lda,n,ipvt,c,0)
CTSS  call dgesl(a,lda,n,ipvt,c,0)
CUNIC call dgesl(a,lda,n,ipvt,c,0)
      anorm = 0.
      amaxx = 0.
      do 30 j =1,n
	 x(j) = c(j)
	 absx = dabs(x(j))
	 anorm = anorm + absx
	 if(absx.gt.amaxx) then
	    amaxx = absx
	    jmax = j
	 end if
 30   continue
C     write(6,*) it,' Maximum is x(',jmax,') = ',x(jmax)
      del = anorm - aold
      aold = anorm
      err = dabs(del/anorm)
      scale = 1.d0/x(jmax)
      do 40 i = 1,n
	 x(i) = x(i)*scale
	 c(i) = x(i)
 40   continue
      if(err.gt.eps .and. it .lt. 30) then
c        do 15 i = 1,n
c            x(i) = 0d0
c            do 17 j = 1,n
c               x(i) = x(i) + b(i,j)*c(j)
c17          continue
c15      continue
c        do 18 i = 1,n
c           c(i) = x(i)
c18      continue
         go to 200
      else
         do 50 i  = 1,n
           c(i)   = x(i)
 50      continue
	 if (it .eq. 30) then
	    ierr = 30
	 end if
      endif
      return
      end
