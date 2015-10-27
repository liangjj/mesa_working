*deck @(#)project.f	5.1  11/6/94
      subroutine project(b,mdim,mvec,c0,npvec,scr,s)
      implicit real*8(a-h,o-z)
      dimension b(mdim,mvec),c0(mdim,npvec),scr(mvec,npvec)
      dimension s(mdim,mvec)
      common /io/ inp,iout
c
c
       call ebtc(scr,b,c0,mvec,mdim,npvec)
c
       call ebct(s,c0,scr,mdim,npvec,mvec)
c
       do 1 i=1,mvec
          do 2 j=1,mdim
             b(j,i)=b(j,i)-s(j,i)
  2       continue
  1    continue
c
       return
       end
