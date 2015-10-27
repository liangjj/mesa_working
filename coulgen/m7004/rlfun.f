*deck rlfun
      function rlfun(eta,l,fact)
c***begin prologue     rlfun
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***description        r function for construction of irregular positive
c                       l energy coulomb functions   
c***references          nbs mathematical handbook (abromowitz and stegun)
c
c***routines called
c
c***end prologue       rlfun
c
      implicit integer (a-z)
      real*8 rlfun, eta, fact, deno, one, two, eulerc
      complex*16 cometa, etafac, num, crl, eye
      dimension fact(0:100)
      common/io/inp,iout
      data zero, one, two /0.d0, 1.d0, 2.d0/
      data eulerc /.577215664901532860606512d0/
      data eye / (0.d0,1.d0) /
      if (l.eq.0) then
          rlfun=zero 
      else
          cometa=eta*eye
          ntrms=l+l+1
          lfac=ntrms
          etafac=cometa-l
          num=one
          crl=num/lfac
          do 10 n=2,ntrms
             lfac=lfac-1
             num=num*two*etafac
             deno=lfac*fact(n-1)
             crl=crl+num/deno
             etafac=etafac+one
   10     continue
          rlfun=imag(crl)
          test=l-2*(l/2)
          if (test.eq.0) then
              rlfun=-rlfun
          endif
          rlfun=rlfun/fact(l+l)
      endif          
      return
      end
