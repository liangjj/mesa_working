*deck %W%  %G%
      subroutine otsmbk(nsym,ml,m1,m2,m3,m4,n)
C
C***Begin prologue
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue
C
      implicit real*8(a-h,o-z)
c
      common /io/ inp,iout
c
c-----------------------------------------------------------------------
c
c --- description     this    routine computes the serial number for
c                      the integral block labeled ml, m1, m2, m3, m4
c                      in a canonical integral list.  note the m values
c                      for cinfv group begins with 0.
c
c-----------------------------------------------------------------------
c
      n = 0
      if (ml .le. 0) go to 500
      do 100 m = 1, ml
         im = nsym - m / 2
 100  n = n + (im * im + im) / 2
c
 500  km = (ml + 1) / 2
      im = m1 - km
      jm = m3 - km
      n = n + (im * im + im) / 2 + jm + 1
      write(iout,99) n
 99   format('  otsmbk   ndex = ',i6)
      return
      end
