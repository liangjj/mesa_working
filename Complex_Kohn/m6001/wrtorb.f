*deck @(#)wrtorb.f	1.1 9/7/91
c***begin prologue     wrtorb
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           wrtorb, link 6001, orbital file
c***author             schneider, barry (lanl)
c***source             m6001
c***purpose            orbital print
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       wrtorb
      subroutine wrtorb(fmo,reg,first,last,nfmax,npnts)
      implicit integer (a-z)
      real *8 fmo
      common /io/ inp,iout
      dimension fmo(npnts,nfmax)
      write (iout,100) reg, npnts
      concnt=0
      do 10 i=first,last
         concnt=concnt+1
         write (iout,200) i
         write(iout,300) (fmo(j,i),j=1,npnts)
   10 continue
  100 format (/,5x,'integration region',1x,i3,1x,'no. points',1x,i6)
  200 format (/,5x,'contracted function',1x,i4)
  300 format ( (/,5x,5(e15.8,1x)) )
      return
      end 
