*deck @(#)derfc.f	5.1 11/6/94
      function derfc(z)
      implicit double precision(a-h,o-z)
      common zah(31),den(31)
      if(z-4.0d0)10,100,100
   10 derfc=z*1.128379167095512d0
      d=derfc
      tz2=z*z*2.0d0
      do 20 n=2,60
      d=d*tz2/dfloat(n+n-1)
   20 derfc=derfc+d
      derfc=1.0d0-derfc*dexp(-z*z)
      return
  100 if(z-13.0)102,102,101
  101 derfc=0.0d0
      return
  102 top=zah(1)
      bot=den(1)
      do 110 n=2,31
      top=top*z+zah(n)
  110 bot=bot*z+den(n)
      derfc=top/bot*0.564189583547756d0*dexp(-z*z)
      return
      end
