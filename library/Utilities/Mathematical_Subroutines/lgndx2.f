*deck lgndx2
c***begin prologue     lgndx2
c***date written       861117   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           gauss-legendre, quadrature
c***author             collins lee (lanl)
c***source             mylib
c***purpose            quadrature
c***description        points and weights for gauss-legendre quadrature
c***                   for 1, 2, 3 points
c***references         none
c
c***routines called    none
      subroutine lgndx2 (n,j,wtg,rtg)
      implicit real *8 (a-h,o-z)
*
      data x11, w11 / .5d0 , 1.d0/
      data x21, w21 /0.5773502692d0,1.d0/
      data x31, x32 /0.7745966692d0,0.0d0/
      data w31, w32 /0.55555556d0,0.88888888d0/
      data half, one /0.5d0,1.0d0/
*
      if (n.eq.1) then
          rtg=x11
          wtg=w11
      else
      if (n.eq.3) go to 20
*
      wtg=w21
      if (j.eq.2) go to 10
      rtg=-x21
      go to 60
   10 rtg=x21
      go to 60
   20 go to (30,40,50), j
   30 wtg=w31
      rtg=-x31
      go to 60
   40 wtg=w32
      rtg=x32
      go to 60
   50 wtg=w31
      rtg=x31
   60 wtg=wtg*half
      rtg=half*(rtg+one)
      endif
      return
      end
