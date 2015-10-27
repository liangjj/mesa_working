*deck rnine
      subroutine rnine(wt,stp)
c***begin prologue     rnine
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
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
c***end prologue       rnine
c
      implicit integer (a-z)
      real*8 wt, stp, d
      dimension wt(9,8)
      data d / 3628800.d0 /
      common /io/ inp, iout
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
      return
      end















