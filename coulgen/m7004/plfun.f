*deck plfun
      function plfun(eta,lval)
c***begin prologue     plfun
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***description        coulomb p  function
c                               l
c***references         NBS handbook
c
c***routines called
c
c***end prologue       plfun
c
      implicit integer (a-z)
      real*8 plfun, eta, etasq, one, two, three, four
      common/io/inp,iout
      data one, two, three, four/1.d0,2.d0,3.d0,4.d0/
      plfun=two*eta
      if (lval.gt.0) then
          etasq=eta*eta
          plfun=plfun*(one+etasq)/three
          if (lval.gt.1) then
              do 10 l=2,lval
                 twoel=l+l
                 twoelp=twoel+1
                 twoelm=twoel-1 
                 plfun=plfun*four*(l*l+etasq)
     1                           / (twoelp*twoel*twoel*twoelm)
   10         continue
          endif
      endif              
      return
      end


