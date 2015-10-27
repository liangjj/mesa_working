*deck honv
      subroutine honv(vout,mat,vin,n)
      implicit integer(a-z)
      real*8 vout, mat, vin
      call ebc(vout,mat,vin,n,n,1)
      return
      end
