*deck @(#)bend.f	5.1  11/6/94
      subroutine bend(noint,i,j,k,b,ib,c)
c***begin prologue     bend.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)bend.f	5.1   11/6/94
c***purpose            
c***description
c        adapted from the normal coordinate analysis program of
c        schachtschneider, shell development .
c***references
c
c***routines called
c
c***end prologue       bend.f
      implicit none
c     --- input variables -----
      integer noint,i,j,k
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer ib(4,2)
      real*8 b(3,4,2),c(*)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer iaind,jaind,kaind,m
      real*8 rji(3),rjk(3),eji(3),ejk(3)
      real*8 zero,one,djisq,djksq,dji,djk,dotj,sinj,test,fuzz,abs
      parameter (zero=0.0d+00,one=1.0d+00,fuzz=1.0d-06)
c
c
      iaind=3*(i-1)
      jaind=3*(j-1)
      kaind=3*(k-1)
      ib(1,noint)=i
      ib(2,noint)=j
      ib(3,noint)=k
      djisq=zero
      djksq=zero
      do 10 m=1,3
         rji(m)=c(m+iaind)-c(m+jaind)
         rjk(m)=c(m+kaind)-c(m+jaind)
         djisq =djisq + rji(m)**2
         djksq =djksq + rjk(m)**2
   10 continue
      dji =sqrt(djisq)
      djk =sqrt(djksq)
c
      dotj = zero
      do 20 m=1,3
         eji(m)=rji(m)/dji
         ejk(m)=rjk(m)/djk
         dotj=dotj+eji(m)*ejk(m)
   20 continue
c
c     --- be careful about 180 degree angles.
      test=dotj*dotj-one
      if(abs(test).le.fuzz) then
         sinj=zero
         call rzero(b(1,1,noint),9)
      else
         sinj =sqrt(one-dotj**2)
         do 30 m=1,3
            b(m,3,noint)=(   (dotj*ejk(m)-eji(m)))/(djk*sinj)
            b(m,1,noint)=(   (dotj*eji(m)-ejk(m)))/(dji*sinj)
            b(m,2,noint)=-b(m,1,noint)-b(m,3,noint)
   30    continue
      endif
c
c
      return
      end
