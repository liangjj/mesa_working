*deck qlplfn
      function qlplfn(eta,l,fact)
c***begin prologue     qlplfn
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***description        q /p  function for construction of irregular positive
c                       l  l energy coulomb functions   
c***references          nbs mathematical handbook (abromowitz and stegun)
c
c***routines called    cpsi(complex psi function:clams)
c
c***end prologue       qlplfn
c
      implicit integer (a-z)
      real*8 qlplfn, eta, fact, zero, one, two, eulerc, etasq
      real*8 sum1, sum2, rlfun, plfun
      complex*16 etafac, eye, cpsi
      dimension fact(0:100)
      common/io/inp,iout
      data zero, one, two /0.d0, 1.d0, 2.d0/
      data eulerc /.577215664901532860606512d0/
      data eye / (0.d0,1.d0) /
      etafac=one+eta*eye
      if (l.eq.0) then
          qlplfn = -one + real( cpsi(etafac) ) + two*eulerc
      else
          etasq=eta*eta
          sum1=zero
          do 10 n=1,l
             sum1=sum1+n/(n*n+etasq)
   10     continue          
          sum2=zero
          do 20 n=1,l+l+1
             sum2=sum2+one/n
   20     continue
          qlplfn = sum1-sum2 + real( cpsi(etafac) ) +two*eulerc +
     1                       rlfun(eta,l,fact)/plfun(eta,l)
      endif          
      return
      end














