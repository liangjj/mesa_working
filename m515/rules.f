*deck @(#)rules.f	5.1  11/28/95
      subroutine rtwo(wt,stp)
c***begin prologue     rules.f
c***date written       920324   (yymmdd)
c***revision date      11/6/94
c***keywords
c***author             schneider, barry (nsf)
c***source             @(#)rules.f	5.1   11/28/95
c***purpose            newton cotes for two point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rules.f
c
      implicit none
      real*8 stp
      real*8 wt(2,1)
      integer inp,iout
c
      common /io/ inp, iout
c
c
      wt(1,1)=stp*.5d0
      wt(2,1)=wt(1,1) 
c
c
      return
      end
      subroutine rthree(wt,stp)
c***begin prologue     rules.f
c***date written       920324   (yymmdd)
c***revision date      11/6/94
c***keywords
c***author             schneider, barry (nsf)
c***source             @(#)rules.f	5.1   11/28/95
c***purpose            newton cotes for three point quadrature
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rules.f
c
      implicit none
      real*8 stp
      real*8 wt(3,2)
      real*8 d
      integer i,j
      integer inp,iout
c
      data d / 12.d0 /
      common /io/ inp, iout
c
      wt(1,1)=5.d0/d
      wt(2,1)=8.d0/d
      wt(3,1)=-1.d0/d
      wt(1,2)=wt(3,1)
      wt(2,2)=wt(2,1)
      wt(3,2)=wt(1,1)
      do 10 i=1,3
         do 20 j=1,2
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
c
c
      return
      end
      subroutine rfour(wt,stp)
c***begin prologue     rules.f
c***date written       920324   (yymmdd)
c***revision date      11/6/94
c***keywords
c***author             schneider, barry (nsf)
c***source             @(#)rules.f	5.1   11/28/95
c***purpose            newton cotes for four point quadrature
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rules.f
c
      implicit none
      real*8 stp
      real*8 wt(4,3)
      real*8 d
      integer i,j
      integer inp,iout
c

      data d / 24.d0 /
      common /io/ inp, iout
c
c
      wt(1,1)=9.d0/d
      wt(2,1)=19.d0/d
      wt(3,1)=-5.d0/d
      wt(4,1)=1.d0/d
      wt(1,2)=-1.d0/d
      wt(2,2)=13.d0/d
      wt(3,2)=13.d0/d
      wt(4,2)=-1.d0/d
      wt(1,3)=1.d0/d
      wt(2,3)=-5.d0/d
      wt(3,3)=19.d0/d
      wt(4,3)=9.d0/d
      do 10 i=1,4
         do 20 j=1,3
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
c
c
      return
      end
      subroutine rfive(wt,stp)
c***begin prologue     rules.f
c***date written       920324   (yymmdd)
c***revision date      11/6/94
c***keywords
c***author             schneider, barry (nsf)
c***source             @(#)rules.f	5.1   11/28/95
c***purpose            newton cotes five point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rules.f
c
      implicit none
      real*8 stp
      real*8 wt(5,4)
      real*8 d
      integer i,j
      integer inp,iout
c
      data d / 720.d0 /
      common /io/ inp, iout
c
c
      wt(1,1)=251.d0/d
      wt(2,1)=646.d0/d
      wt(3,1)=-264.d0/d
      wt(4,1)=106.d0/d
      wt(5,1)=-19.d0/d
      wt(1,2)=-19.d0/d
      wt(2,2)=346.d0/d
      wt(3,2)=456.d0/d
      wt(4,2)=-74.d0/d
      wt(5,2)=11.d0/d
      wt(1,3)=11.d0/d
      wt(2,3)=-74.d0/d
      wt(3,3)=456.d0/d
      wt(4,3)=346.d0/d
      wt(5,3)=-19.d0/d
      wt(1,4)=-19.d0/d
      wt(2,4)=106.d0/d
      wt(3,4)=-264.d0/d
      wt(4,4)=646.d0/d
      wt(5,4)=251.d0/d
      do 10 i=1,5
         do 20 j=1,4
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
c
c
      return
      end
      subroutine rsix(wt,stp)
c***begin prologue     rules.f
c***date written       920324   (yymmdd)
c***revision date      11/6/94
c***keywords
c***author             schneider, barry (nsf)
c***source             @(#)rules.f	5.1   11/28/95
c***purpose            newton cotes six point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rules.f
      implicit none
      real*8 stp
      real*8 wt(6,5)
      real*8 d
      integer i,j
      integer inp,iout
c
      data d / 1440.d0 /
      common /io/ inp, iout
c
      wt(1,1)=475.d0/d
      wt(2,1)=1427.d0/d
      wt(3,1)=-798.d0/d
      wt(4,1)=482.d0/d
      wt(5,1)=-173.d0/d
      wt(6,1)=27.d0/d
      wt(1,2)=-27.d0/d
      wt(2,2)=637.d0/d
      wt(3,2)=1022.d0/d
      wt(4,2)=-258.d0/d
      wt(5,2)=77.d0/d
      wt(6,2)=-11.d0/d
      wt(1,3)=11.d0/d
      wt(2,3)=-93.d0/d
      wt(3,3)=802.d0/d
      wt(4,3)=802.d0/d
      wt(5,3)=-93.d0/d
      wt(6,3)=11.d0/d
      wt(1,4)=-11.d0/d
      wt(2,4)=77.d0/d
      wt(3,4)=-258.d0/d
      wt(4,4)=1022.d0/d
      wt(5,4)=637.d0/d
      wt(6,4)=-27.d0/d
      wt(1,5)=27.d0/d
      wt(2,5)=-173.d0/d
      wt(3,5)=482.d0/d
      wt(4,5)=-798.d0/d
      wt(5,5)=1427.d0/d
      wt(6,5)=475.d0/d
      do 10 i=1,6
         do 20 j=1,5
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
c
c
      return
      end
      subroutine rseven(wt,stp)
c***begin prologue     rules.f
c***date written       920324   (yymmdd)
c***revision date      11/6/94
c***keywords
c***author             schneider, barry (nsf)
c***source             @(#)rules.f	5.1   11/28/95
c***purpose            newton cotes seven point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rules.f
      implicit none
      real*8 stp
      real*8 wt(7,6)
      real*8 d
      integer i,j
      integer inp,iout
c

      data d / 60480.d0 /
      common /io/ inp, iout
c
c
      wt(1,1)=19087.d0/d
      wt(2,1)=65112.d0/d
      wt(3,1)=-46461.d0/d
      wt(4,1)=37504.d0/d
      wt(5,1)=-20211.d0/d
      wt(6,1)=6312.d0/d
      wt(7,1)=-863.d0/d
      wt(1,2)=-863.d0/d
      wt(2,2)=25128.d0/d
      wt(3,2)=46989.d0/d
      wt(4,2)=-16256.d0/d
      wt(5,2)=7299.d0/d
      wt(6,2)=-2088.d0/d
      wt(7,2)=271.d0/d
      wt(1,3)=271.d0/d
      wt(2,3)=-2760.d0/d
      wt(3,3)=30819.d0/d
      wt(4,3)=37504.d0/d
      wt(5,3)=-6771.d0/d
      wt(6,3)=1608.d0/d
      wt(7,3)=-191.d0/d
      wt(1,4)=-191.d0/d
      wt(2,4)=1608.d0/d
      wt(3,4)=-6771.d0/d
      wt(4,4)=37504.d0/d
      wt(5,4)=30819.d0/d
      wt(6,4)=-2760.d0/d
      wt(7,4)=271.d0/d
      wt(1,5)=271.d0/d
      wt(2,5)=-2088.d0/d
      wt(3,5)=7299.d0/d
      wt(4,5)=-16256.d0/d
      wt(5,5)=46989.d0/d
      wt(6,5)=25128.d0/d
      wt(7,5)=-863.d0/d
      wt(1,6)=-863.d0/d
      wt(2,6)=6312.d0/d
      wt(3,6)=-20211.d0/d
      wt(4,6)=37504.d0/d
      wt(5,6)=-46461.d0/d
      wt(6,6)=65112.d0/d
      wt(7,6)=19087.d0/d
      do 10 i=1,7
         do 20 j=1,6
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
c
c
      return
      end
      subroutine reight(wt,stp)
c***begin prologue     rules.f
c***date written       920324   (yymmdd)
c***revision date      11/6/94
c***keywords
c***author             schneider, barry (nsf)
c***source             @(#)rules.f	5.1   11/28/95
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
c***end prologue       rules.f
      implicit none
      real*8 stp
      real*8 wt(8,7)
      real*8 d
      integer i,j
      integer inp,iout
c
      data d / 120960.d0 /
      common /io/ inp, iout
c
c
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
c
c
      return
      end
      subroutine rnine(wt,stp)
c***begin prologue     rules.f
c***date written       920324   (yymmdd)
c***revision date      11/6/94
c***keywords
c***author             schneider, barry (nsf)
c***source             @(#)rules.f	5.1   11/28/95
c***purpose            newton cotes nine point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rules.f
      implicit none
      real*8 stp
      real*8 wt(9,8)
      real*8 d
      integer i,j
      integer inp,iout
c
      data d / 3628800.d0 /
      common /io/ inp, iout
c
c
      wt(1,1)=1070017.d0/d
      wt(2,1)=4467094.d0/d
      wt(3,1)=-4604594.d0/d
      wt(4,1)=5595358.d0/d
      wt(5,1)=-5033120.d0/d
      wt(6,1)=3146338.d0/d
      wt(7,1)=-1291214.d0/d
      wt(8,1)=312874.d0/d
      wt(9,1)=-33953.d0/d
      wt(1,2)=-33953.d0/d
      wt(2,2)=1375594.d0/d
      wt(3,2)=3244786.d0/d
      wt(4,2)=-1752542.d0/d
      wt(5,2)=1317280.d0/d
      wt(6,2)=-755042.d0/d
      wt(7,2)=294286.d0/d
      wt(8,2)=-68906.d0/d
      wt(9,2)=7297.d0/d
      wt(1,3)=7297.d0/d
      wt(2,3)=-99626.d0/d
      wt(3,3)=1638286.d0/d
      wt(4,3)=2631838.d0/d
      wt(5,3)=-833120.d0/d
      wt(6,3)=397858.d0/d
      wt(7,3)=-142094.d0/d
      wt(8,3)=31594.d0/d
      wt(9,3)=-3233.d0/d
      wt(1,4)=-3233.d0/d
      wt(2,4)=36394.d0/d
      wt(3,4)=-216014.d0/d
      wt(4,4)=1909858.d0/d
      wt(5,4)=2224480.d0/d
      wt(6,4)=-425762.d0/d
      wt(7,4)=126286.d0/d
      wt(8,4)=-25706.d0/d
      wt(9,4)=2497.d0/d
      wt(1,5)=2497.d0/d
      wt(2,5)=-25706.d0/d
      wt(3,5)=126286.d0/d
      wt(4,5)=-425762.d0/d
      wt(5,5)=2224480.d0/d
      wt(6,5)=1909858.d0/d
      wt(7,5)=-216014.d0/d
      wt(8,5)=36394.d0/d
      wt(9,5)=-3233.d0/d
      wt(1,6)=-3233.d0/d
      wt(2,6)=31594.d0/d
      wt(3,6)=-142094.d0/d
      wt(4,6)=397858.d0/d
      wt(5,6)=-833120.d0/d
      wt(6,6)=2631838.d0/d
      wt(7,6)=1638286.d0/d
      wt(8,6)=-99626.d0/d
      wt(9,6)=7297.d0/d 
      wt(9,7)=-33953.d0/d
      wt(8,7)=1375594.d0/d
      wt(7,7)=3244786.d0/d
      wt(6,7)=-1752542.d0/d
      wt(5,7)=1317280.d0/d
      wt(4,7)=-755042.d0/d
      wt(3,7)=294286.d0/d
      wt(2,7)=-68906.d0/d
      wt(1,7)=7297.d0/d     
      wt(9,8)=1070017.d0/d
      wt(8,8)=4467094.d0/d
      wt(7,8)=-4604594.d0/d
      wt(6,8)=5595358.d0/d
      wt(5,8)=-5033120.d0/d
      wt(4,8)=3146338.d0/d
      wt(3,8)=-1291214.d0/d
      wt(2,8)=312874.d0/d
      wt(1,8)=-33953.d0/d
      do 10 i=1,9
         do 20 j=1,8
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
c
c
      return
      end
