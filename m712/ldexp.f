*deck @(#)ldexp.f	5.1  11/6/94
      subroutine ldexp(a,asave,ex,nprim,nij,i1,i2,j1,j2,cen1,cen3,
     #                 icen,jcen,c,t1,xyza,xyzam1,xyzam3,nat,ni,indx)
c
c***purpose: calculate as many intermediate arrays as possible
c    with the locations and exponents of just two atoms. the
c    arrays are then used in the computation of the two-electron
c    integrals
c
c paul saxe                       8 august 1984              lanl
c
c 17 november 1985   modified by pws at lanl to add indx array.
c
      implicit integer (a-z)
c
      real*8 a(nij),asave(nij),ex(nprim),c(3,nat),xyza(nij,3)
      real*8 xyzam1(nij,3),xyzam3(nij,3),t1(ni)
      real*8 rijsq,scal1,scal2,scal3,scal4
      integer indx(nij,2)
c
c     ----- timing -----
c
c
c     ----- form a = alpha(i) + alpha(j) and
c            asave = alpha(i) * alpha(j) * r(ij)**2 / a
c
      rijsq=(c(1,icen)-c(1,jcen))**2+(c(2,icen)-c(2,jcen))**2+
     #      (c(3,icen)-c(3,jcen))**2
      ij=0
      do 2 j=j1,j2
         scal1=ex(j)
         scal2=scal1*rijsq
         do 1 i=i1,i2
            ij=ij+1
            a(ij)=ex(i)+scal1
            asave(ij)=ex(i)*scal2/a(ij)
            indx(ij,1)=i-i1+1
            indx(ij,2)=j-j1+1
    1    continue
    2 continue
c
c     ----- form xa, ya and za = [alpha(i)*x(i) + alpha(j)*x(j)] / a -----
c
      do 6 coord=1,3
         scal2=c(coord,icen)
         scal3=c(coord,cen1)
         scal4=c(coord,cen3)
         do 3 i=i1,i2
            t1(i-i1+1)=scal2*ex(i)
    3    continue
         ij=0
         do 5 j=j1,j2
            scal1=c(coord,jcen)*ex(j)
            do 4 i=i1,i2
               ij=ij+1
               xyza(ij,coord)=(scal1+t1(i-i1+1))/a(ij)
               xyzam1(ij,coord)=(xyza(ij,coord)-scal3)*a(ij)
               xyzam3(ij,coord)=(xyza(ij,coord)-scal4)*a(ij)
    4       continue
    5    continue
    6 continue
c
c     ----- timing -----
c
c
c
      return
      end
