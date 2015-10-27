*deck @(#)symlin.f	5.1  11/6/94
      subroutine symlin(am,nv,vec,ising)
      implicit none
c
      integer nv,ising
      real*8 am(1),vec(1)
c
      integer i,j,ijmul
      real*8 zero
c
      parameter (zero=0.0d+00)
c
      call syminv(am,nv,ising)
c     --- if matrix is nonsingular
      if (ising.eq.0) then
         do 30 i=1,nv
            am(i)=zero
            do 20 j=1,nv
               if (j.ge.i) then
                  ijmul=i*nv+j-(i*(i-1))/2
               else
                  ijmul=j*nv+i-(j*(j-1))/2
               endif
               am(i)=am(i)+am(ijmul)*vec(j)
   20       continue
   30    continue
      endif
c
c
      return
      end
