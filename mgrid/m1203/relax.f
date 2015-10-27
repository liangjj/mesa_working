*deck relax
c***begin prologue     relax
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            gauss-seidel iteration
c***description        we consider planes which are stacked on top of one
c***                   another.  each plane consists of alternating points of
c**                    two colors.  we label plane one as a red/black plane 
c***                   and the one above (below) it as a blue/green plane.
c***                   if we consider red/green points it is easy to see that 
c***                   second differencing leads to coupling to only 
c***                   blue/black. so we do all of those first.  
c***                   then we return to finish the blue/black using the 
c***                   updated red/green.             
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
      dimension u(n,n,n), rhs(n,n,n)
      common/io/inp, iout
c     update the red points.  they depend only on their neighboring black
c     points in the same plane and the blue points in the plane above.      
      do 10 k=2,n-1,2
         ii=2
         jj=3
         do 20 j=2,n-1
            do 30 i=ii,n-1,2
               u(i,j,k)=ifac*( u(i+1,j,k) + u(i-1,j,k) + 
     1                         u(i,j+1,k) + u(i,j-1,k) +
     2                         u(i,j,k+1) + u(i,j,k-1) - 
     3                              hsq*rhs(i,j,k) )
 30         continue
            itmp=ii
            ii=jj
            jj=itmp 
 20      continue
 10   continue
c     update the green points.  same situation as previous loops.
      do 40 k=3,n-1,2
         ii=3
         jj=2
         do 50 j=2,n-1
            do 60 i=ii,n-1,2
               u(i,j,k)=ifac*( u(i+1,j,k) + u(i-1,j,k) + 
     1                         u(i,j+1,k) + u(i,j-1,k) +
     2                         u(i,j,k+1) + u(i,j,k-1) - 
     3                              hsq*rhs(i,j,k) )
   60       continue
            itmp=ii
            ii=jj
            jj=itmp
   50    continue
   40 continue
c     now that we have updated all of the red and green points we return to
c     to the black and blue points.
c     first we do the black points
      do 70 k=2,n-1,2
         ii=3
         jj=2
         do 80 j=2,n-1
            do 90 i=ii,n-1,2
               u(i,j,k)=ifac*( u(i+1,j,k) + u(i-1,j,k) + 
     1                         u(i,j+1,k) + u(i,j-1,k) +
     2                         u(i,j,k+1) + u(i,j,k-1) - 
     3                              hsq*rhs(i,j,k) )
 90         continue
            itmp=ii
            ii=jj
            jj=itmp 
 80      continue
 70   continue
c     now finish the job with the blue points
      do 100 k=3,n-1,2
         ii=2
         jj=3
         do 200 j=2,n-1
            do 300 i=ii,n-1,2
               u(i,j,k)=ifac*( u(i+1,j,k) + u(i-1,j,k) + 
     1                         u(i,j+1,k) + u(i,j-1,k) +
     2                         u(i,j,k+1) + u(i,j,k-1) - 
     3                              hsq*rhs(i,j,k) )
  300       continue
            itmp=ii
            ii=jj
            jj=itmp 
  200    continue
  100 continue                         
      return
      end





