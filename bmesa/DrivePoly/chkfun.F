*deck chkfun.f 
c***begin prologue     chkfun
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             chkfun
c***purpose            expansion of functions and derivatives in coordinate functions
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       chkfun
      subroutine chkfun(cp,dcp,ddcp,pt,wt,fex,dfex,ddfex,
     1                  fapp,dfapp,ddfapp,c,work,type,parity,n)
c
      implicit integer (a-z)
      character*(*) type, parity
      real*8 cp, dcp, ddcp, pt, wt, fex, dfex, ddfex 
      real*8 fapp, dfapp, ddfapp, c, work
      dimension cp(n,n), dcp(n,n), ddcp(n,n), pt(n), wt(n)
      dimension fex(n), dfex(n), ddfex(n)
      dimension fapp(n), dfapp(n), ddfapp(n), c(n), work(n,n)
      common/io/inp, iout
c
      if(type.eq.'x') then
         do 10 i=1,n
            fex(i)=pt(i)
            dfex(i)=1.d0
            ddfex(i)=0.d0
 10      continue
      elseif(type.eq.'xx') then
         do 20 i=1,n
            fex(i)=pt(i)*pt(i)
            dfex(i)=2.d0*pt(i)
            ddfex(i)=2.d0
 20      continue
      elseif(type.eq.'xxx') then
         do 30 i=1,n
            fex(i)=pt(i)*pt(i)*pt(i)
            dfex(i)=3.d0*pt(i)*pt(i) 
            ddfex(i)=6.d0*pt(i)
 30      continue   
      elseif(type.eq.'xxxx') then
         do 40 i=1,n
            fex(i)=pt(i)*pt(i)*pt(i)*pt(i)
            dfex(i)=4.d0*pt(i)*pt(i)*pt(i) 
            ddfex(i)=12.d0*pt(i)*pt(i)
 40      continue   
      elseif(type.eq.'negative-exponential') then
         do 50 i=1,n
            fex(i) = exp(-pt(i))
            dfex(i)=-fex(i)
            ddfex(i)=fex(i)
 50      continue
      elseif(type.eq.'positive-exponential') then
         do 60 i=1,n
            fex(i) = exp(pt(i))
            dfex(i)=fex(i)
            ddfex(i)=fex(i)
 60      continue
      elseif(type.eq.'negative-gaussian') then
         do 70 i=1,n
            fex(i)=exp(-pt(i)*pt(i))
            dfex(i)=-2.d0*pt(i)*fex(i)
            ddfex(i)=-2.d0*pt(i)*dfex(i) - 2.d0*fex(i)
 70      continue
      elseif(type.eq.'positive-gaussian') then
         do 80 i=1,n
            fex(i)=exp(pt(i)*pt(i))
            dfex(i)=2.d0*pt(i)*fex(i)
            ddfex(i)=2.d0*pt(i)*dfex(i) + 2.d0*fex(i)
 80      continue
      elseif(type.eq.'sine') then
         do 90 i=1,n
            fex(i)=sin(pt(i))
            dfex(i)=cos(pt(i))
            ddfex(i)=-fex(i)
 90      continue
      elseif(type.eq.'cosine') then
         do 100 i=1,n
            fex(i)=cos(pt(i))
            dfex(i)=-sin(pt(i))
            ddfex(i)=-fex(i)
 100        continue
      elseif(type.eq.'one') then
         do 200 i=1,n
	    fex(i) = 1.d0
	    dfex(i) = 0.d0
	    ddfex(i) = 0.d0
 200     continue
      else	       
         call lnkerr('error')
      endif   
      do 300 i=1,n
         c(i) = cp(i,i)*wt(i)*fex(i)
 300  continue
      call ebc(fapp,cp,c,n,n,1)
      if(parity.ne.'none') then
         do 400 i=1,n
            do 410 j=1,n
               work(j,i) = 2.d0*pt(j)*dcp(j,i)
 410        continue
 400     continue   
         call ebc(dfapp,work,c,n,n,1)
         do 500 i=1,n
            do 510 j=1,n
               work(j,i) = 2.d0*dcp(j,i) + 4.d0*pt(j)*pt(j)*ddcp(j,i) 
 510        continue
 500     continue
         call ebc(ddfapp,work,c,n,n,1)   
      else
         call ebc(dfapp,dcp,c,n,n,1)
         call ebc(ddfapp,ddcp,c,n,n,1)
      endif
      write(iout,1) type  
      do 600 i=1,n
         write(iout,2) pt(i), fex(i), fapp(i)
 600  continue
      write(iout,3) type
      do 700 i=1,n
         write(iout,2) pt(i), dfex(i), dfapp(i)  
 700  continue
      write(iout,4) type
      do 800 i=1,n
         write(iout,2) pt(i), ddfex(i), ddfapp(i)   
 800  continue   
      return
 1    format(/,'comparison of exact and approximate function = ',a24)
 2    format(2x,'point = ',e15.8,1x,'exact = ',e15.8,1x,
     1          'approximate = ',e15.8)
 3    format(/,'comparison of exact and approximate derivative = ',a24)
 4    format(/,'comparison of exact and approximate second derivative '
     1         '= ',a24)
      end


















