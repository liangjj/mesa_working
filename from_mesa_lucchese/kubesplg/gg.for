      function gg(x)
      implicit real*8 (a-h,o-z)
      real*8 mu
      parameter (nsmall=4015,lbig=#maxltop)
      common /integ/ xstart,ystart,ypstart,al,mu,h,xmax
      common /energy/ e,eta,ak
      common/spl/rd,r(nsmall),w(nsmall,0:lbig),cs(nsmall),rd26,l
      klo=(x-xstart)/rd
      klo=klo+1
      a=(r(klo+1)-x)/rd
      b=(x-r(klo))/rd
      gg=a*w(klo,l)+b*w(klo+1,l)+(a*(a*a-1.)*cs(klo)+b*
     1 (b*b-1.)*cs(klo+1))*rd26
      return
      end
