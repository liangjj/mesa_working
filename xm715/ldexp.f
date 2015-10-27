*deck %W%  %G%
      subroutine ldexp(a,asave,ex,nprim,nij,i1,i2,j1,j2,cen1,cen3,
     $                 icen,jcen,c,t1,xyza,xyzam1,xyzam3,nat,ni,indx)
c***begin prologue     ldexp.f
c***date written       840808  
c***revision date      11/6/94      
c
c   17 november,1985   pws at lanl
c      modifying to add indx array
c***keywords           
c***author             saxe, paul
c***source             %W%   %G%
c***purpose            calculate as many intermediate arrays as
c    possible with the locations and exponents of just two atoms.
c    the arrays are then used in the computation of the 
c    two-electron integrals
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       ldexp.f
      implicit none
c     --- input variables -----
      integer nprim,nij,i1,i2,j1,j2
      integer nat,ni
      integer cen1,cen3,icen,jcen
c     --- input arrays (unmodified) ---
      real*8 ex(nprim),c(3,nat)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 a(nij),asave(nij),xyza(nij,3)
      real*8 xyzam1(nij,3),xyzam3(nij,3),t1(ni)
      integer indx(nij,2)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,j,ij,coord
      real*8 rijsq,scal1,scal2,scal3,scal4
c
c     --- form a = alpha(i) + alpha(j) and
c         asave = alpha(i) * alpha(j) * r(ij)**2 / a
c
      rijsq=(c(1,icen)-c(1,jcen))**2+(c(2,icen)-c(2,jcen))**2+
     $      (c(3,icen)-c(3,jcen))**2
      ij=0
      do 20 j=j1,j2
         scal1=ex(j)
         scal2=scal1*rijsq
         do 10 i=i1,i2
            ij=ij+1
            a(ij)=ex(i)+scal1
            asave(ij)=ex(i)*scal2/a(ij)
            indx(ij,1)=i-i1+1
            indx(ij,2)=j-j1+1
   10    continue
   20 continue
c
c     --- form xa, ya and za = [alpha(i)*x(i) + alpha(j)*x(j)] / a ---
c
      do 60 coord=1,3
         scal2=c(coord,icen)
         scal3=c(coord,cen1)
         scal4=c(coord,cen3)
         do 30 i=i1,i2
            t1(i-i1+1)=scal2*ex(i)
   30    continue
         ij=0
         do 50 j=j1,j2
            scal1=c(coord,jcen)*ex(j)
            do 40 i=i1,i2
               ij=ij+1
               xyza(ij,coord)=(scal1+t1(i-i1+1))/a(ij)
               xyzam1(ij,coord)=(xyza(ij,coord)-scal3)*a(ij)
               xyzam3(ij,coord)=(xyza(ij,coord)-scal4)*a(ij)
   40       continue
   50    continue
   60 continue
c
c
      return
      end
