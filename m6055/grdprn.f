*deck @(#)grdprn.f	1.1 9/8/91
c***begin prologue     grdprn
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, print grid
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print co-ordinate grid
c***references         none
c
c***routines called    none
c***end prologue       grdprn
      subroutine grdprn(grid,npnts,reg)
      implicit integer (a-z)
      common /io/ inp, iout
      real *8 grid
      dimension grid(4,npnts)
      write (iout,10) reg
      write (iout,20)
      do 30 grpt=1,npnts
         write (iout,40) grid(1,grpt), grid(2,grpt), grid(3,grpt)
   30 continue
   40 format(5x,3f10.6)
      return
   10 format (/,5x,'grid points for region',1x,i4)
   20 format(/,5x,'     x     ','     y     ','     z     ')
      end
