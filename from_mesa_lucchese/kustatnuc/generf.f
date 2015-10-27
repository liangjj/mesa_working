      subroutine generf (maxtyp,maxrng)
      implicit real*8 (a-h,o-z)
      common/store/str0(350),str1(350),str2(350),str3(350),str4(350),
     1str5(350),str6(350),str7(350),str8(350),str9(350),str10(350)
      x1 = 1.e0
      x3 = 1.e0/3.e0
      x5 = 1.e0/5.e0
      x7 = 1.e0/7.e0
      x9 = 1.e0/9.e0
      x11= 1.e0/11.e0
      x13= 1.e0/13.e0
      x15= 1.e0/15.e0
      x17= 1.e0/17.e0
      x19= 1.e0/19.e0
      maxrng=350
      t=0.e0
      do 779 i=1,maxrng
      y= exp(-t)
      u=2.e0*t
      str10(i) = fmch(10,t,y)
      str9(i)=(u*str10(i)+y)*x19
      str8(i)=(u*str9(i)+y)*x17
  766 str7(i)=(u*str8(i)+y)*x15
      str6(i)=(u*str7(i)+y)*x13
      str5(i)=(u*str6(i)+y)*x11
      str4(i)=(u*str5(i)+y)*x9
      str3(i)=(u*str4(i)+y)*x7
      str2(i)=(u*str3(i)+y)*x5
      str1(i)=(u*str2(i)+y)*x3
      str0(i)=(u*str1(i)+y)
  779 t=t+0.1e0
      return
      end
