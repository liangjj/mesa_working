*deck %W%  %G%
      subroutine eigout(a,b,nad,nbd)
C
C***Begin prologue     eigout
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
C***End prologue       eigout
C
      implicit real*8 (a-h,o-z)
c
      dimension a(nad,nbd),b(2)
      character*1 line
c
      common /io/ inp,iout
c
      data line /'-'/
c
 1    format(2x,5(7x,i5))
 2    format(2x,126a1)
 3    format(2x,i2,2x,5f12.7)
 4    format(/,' eval  ',5f12.7)
 44   format(/,' occ   ',5f12.7)
 5    format(/)
c
      ii=0
      jj=0
      m=nad
      n=nbd
 200  ii=ii+1
      jj=jj+1
      kk=5*jj
      nn=n
      if(n.gt.kk) nn=kk
      ll=(nn-ii+1)
      write(iout,1) (i,i=ii,nn)
      write(iout,4) (b(j),j=ii,nn)
      write(iout,2) (line,i=1,ll*5+4)
      do 101 i=1,m
 101  write(iout,3) i,(a(i,j),j=ii,nn)
      if(n.le.kk) go to 201
      write(iout,5)
      ii=kk
      go to 200
 201  return
      end
