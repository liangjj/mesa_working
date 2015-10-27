*deck @(#)generf.f	1.1 9/7/91
c***begin prologue     generf
c***date written                (yymmdd)
c***revision date      890417   (yymmdd)
c***keywords           generf, link 6003, error function
c***author             unknown
c***source             m6002
c***purpose            special functions for gaussian integrals
c*** 
c
c***references       
c
c***routines called    fmch
c***end prologue       generf
      subroutine generf
      parameter (maxt=8 , maxr=350)
      implicit real *8 (a-h,o-z)
      common/store/str0(maxr),str1(maxr),str2(maxr),str3(maxr),
     1             str4(maxr),str5(maxr),str6(maxr),str7(maxr),
     2             str8(maxr),str9(maxr),str10(maxr)
      x1 = 1.d0
      x3 = 1.d0/3.d0
      x5 = 1.d0/5.d0
      x7 = 1.d0/7.d0
      x9 = 1.d0/9.d0
      x11= 1.d0/11.d0
      x13= 1.d0/13.d0
      x15= 1.d0/15.d0
      x17= 1.d0/17.d0
      x19= 1.d0/19.d0
      t=0.d0
      do 779 i=1,maxr
         y= exp(-t)
         u=2.d0*t
         str10(i) = fmch(10,t,y)
         str9(i)=(u*str10(i)+y)*x19
         str8(i)=(u*str9(i)+y)*x17
         str7(i)=(u*str8(i)+y)*x15
         str6(i)=(u*str7(i)+y)*x13
         str5(i)=(u*str6(i)+y)*x11
         str4(i)=(u*str5(i)+y)*x9
         str3(i)=(u*str4(i)+y)*x7
         str2(i)=(u*str3(i)+y)*x5
         str1(i)=(u*str2(i)+y)*x3
         str0(i)=(u*str1(i)+y)
         t=t+0.1d0
  779 continue
      return
      end 
