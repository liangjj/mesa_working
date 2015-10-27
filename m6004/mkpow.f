      subroutine mkpow(ipow,lmax)
      implicit integer (a-z)
      dimension ipow(0:lmax)
      data zero, one/ 0, 1/
      ipow(0)=zero
      do 10 i=1,lmax
         ipow(i)=one
   10 continue
      return
      end
