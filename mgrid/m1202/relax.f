*deck relax
c***begin prologue     relax
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            red/black gauss-seidel iteration
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       relax
      subroutine relax(u,rhs,hsq,ifac,n)
c
      implicit integer (a-z)
      real*8 u, rhs, hsq, ifac
      character*(*) type
      character*80 title
      dimension u(n,n), rhs(n,n)
      common/io/inp, iout
c      title='input rhs'
c      call prntrm(title,rhs,n,n,n,n,iout)
c      title='initial solution'
c      call prntrm(title,u,n,n,n,n,iout)
c     the points are labelled as red and black.  the red points depend for
c     relaxation on the black points only.  so we update these using our
c     black initial guess.  once we have those, we use the updated red values
c     to get our final updated black points.  
      ii=2
      jj=3
      do 10 j=2,n-1
         do 20 i=ii,n-1,2
            u(i,j)=ifac*( u(i+1,j) + u(i-1,j) + u(i,j+1) 
     1                             + u(i,j-1)
     2                             - hsq*rhs(i,j) )
   20    continue
         itmp=ii
         ii=jj
         jj=itmp
   10 continue
      ii=3
      jj=2
      do 30 j=2,n-1
         do 40 i=ii,n-1,2
            u(i,j)=ifac*( u(i+1,j) + u(i-1,j) + u(i,j+1) 
     1                             + u(i,j-1)
     2                             - hsq*rhs(i,j) )
   40    continue
         itmp=ii
         ii=jj
         jj=itmp
   30 continue                                   
c      title='relaxing'
c      call prntrm(title,u,n,n,n,n,iout)  
      return
      end





