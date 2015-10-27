*deck reight
      subroutine reight(wt,stp)
c***begin prologue     reight
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            newton cotes eight point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       reight
c
      implicit integer (a-z)
      real*8 wt, stp, d
      dimension wt(8,7)
      data d / 120960.d0 /
      common /io/ inp, iout
      wt(1,1)=36799.d0/d
      wt(2,1)=139849.d0/d
      wt(3,1)=-121797.d0/d
      wt(4,1)=123133.d0/d
      wt(5,1)=-88547.d0/d
      wt(6,1)=41499.d0/d
      wt(7,1)=-11351.d0/d
      wt(8,1)=1375.d0/d
      wt(1,2)=-1375.d0/d
      wt(2,2)=47799.d0/d
      wt(3,2)=101349.d0/d
      wt(4,2)=-44797.d0/d
      wt(5,2)=26883.d0/d
      wt(6,2)=-11547.d0/d
      wt(7,2)=2999.d0/d
      wt(8,2)=-351.d0/d
      wt(1,3)=351.d0/d
      wt(2,3)=-4183.d0/d
      wt(3,3)=57627.d0/d
      wt(4,3)=81693.d0/d
      wt(5,3)=-20227.d0/d
      wt(6,3)=7227.d0/d
      wt(7,3)=-1719.d0/d
      wt(8,3)=191.d0/d
      wt(1,4)=-191.d0/d
      wt(2,4)=1879.d0/d
      wt(3,4)=-9531.d0/d
      wt(4,4)=68323.d0/d
      wt(5,4)=68323.d0/d
      wt(6,4)=-9531.d0/d
      wt(7,4)=1879.d0/d
      wt(8,4)=-191.d0/d
      wt(1,5)=191.d0/d
      wt(2,5)=-1719.d0/d
      wt(3,5)=7227.d0/d
      wt(4,5)=-20227.d0/d
      wt(5,5)=81693.d0/d
      wt(6,5)=57627.d0/d
      wt(7,5)=-4183.d0/d
      wt(8,5)=351.d0/d
      wt(1,6)=-351.d0/d
      wt(2,6)=2999.d0/d
      wt(3,6)=-11547.d0/d
      wt(4,6)=26883.d0/d
      wt(5,6)=-44797.d0/d
      wt(6,6)=101349.d0/d
      wt(7,6)=47799.d0/d
      wt(8,6)=-1375.d0/d
      wt(1,7)=1375.d0/d
      wt(2,7)=-11351.d0/d
      wt(3,7)=41499.d0/d
      wt(4,7)=-88547.d0/d
      wt(5,7)=123133.d0/d
      wt(6,7)=-121797.d0/d
      wt(7,7)=139849.d0/d
      wt(8,7)=36799.d0/d
      do 10 i=1,8
         do 20 j=1,7
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
      return
      end















