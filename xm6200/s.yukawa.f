h31541
s 00042/00000/00000
d D 1.1 94/02/16 20:35:14 mesa 1 0
c date and time created 94/02/16 20:35:14 by mesa
e
u
U
f e 0
t
T
I 1
*deck %W%  %G%
c***begin prologue     yukawa
c***date written       930502   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           yukawa, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            calculate yukawa potential
c***references         
c
c***routines called
c***end prologue       yukawa
      subroutine yukawa (yuk,grid,eta,cen,nr,nthet,nphi,ncen,
     1                   nang,nonsep)
      implicit integer (a-z)
      real*8 yuk, grid, eta, cen, dist
      logical nonsep
      dimension yuk(*), grid(3,*), eta(ncen), cen(3,ncen)
      common /io/ inp, iout
      count=0
      if (nonsep) then
          nprd=nang
      else
          nprd=nthet*nphi
      endif
      do 10 i=1,nr
         do 20 j=1,nprd
            count=count+1
            do 30 nc=1,ncen
               dist=sqrt(  (grid(1,count)-cen(1,nc))*
     1                     (grid(1,count)-cen(1,nc)) +
     2                     (grid(2,count)-cen(2,nc))*
     3                     (grid(2,count)-cen(2,nc)) +
     4                     (grid(3,count)-cen(3,nc))*
     5                     (grid(3,count)-cen(3,nc)) )
               yuk(count)=yuk(count) + 
     1                                exp(-eta(nc)*dist)/dist
   30       continue
   20    continue
   10 continue
      return
      end
E 1
