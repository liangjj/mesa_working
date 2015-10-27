c----------------------------------------------------------------------c
c                   ricatti-bessel or bessel function                  c
c                           programs                                   c
c----------------------------------------------------------------------c
*deck rbes
      subroutine rbes (type,l,z,b,bp,bpp,y,yp,ypp)
      implicit real*8 (a-h,o-z)
      character *(*) type
      call bffgh (z,l,b,bp,bpp,y,yp,ypp,z*z,2,1.d0)
c----------------------------------------------------------------------c
c                make ricatti-bessels if you want                      c
c----------------------------------------------------------------------c   
      if (type.eq.'ricatti-bessel') then
          bpp=z*bpp+2.d+00*bp
          ypp=z*ypp+2.d+00*yp
          bp=z*bp+b
          yp=z*yp+y
          b=z*b
          y=z*y
      endif
      return
      end
