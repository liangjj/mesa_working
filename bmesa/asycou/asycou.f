*deck asycou.f 
c***begin prologue     asycou
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           coulomb functions
c***                   
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            driver for asymptotic package for
c***                   coulomb functions.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       asycou
      program asycou
c
      implicit integer (a-z)
      character*4096 ops
      logical dollar, logkey
      logical prnt
      real*8 charge, energy, fpkey, k, eta, r, rho 
      real*8 z, fl, gl, dfl, dgl, wron, fl0, gl0, dfl0, dgl0, wron0
      common/io/inp, iout      
      pointer (pz,z(1))
c
      call drum
      write(iout,1)
      call iosys ('read character options from rwf',-1,0,0,ops)
      prn=logkey(ops,'print=m6297=expansion-coefficients',.false.,' ')
      charge=fpkey(ops,'charge',-1.d0,' ')
      energy=fpkey(ops,'energy',1.d0,' ')
      l=intkey(ops,'angular-momentum',0,' ')
      ntrms=intkey(ops,'number-of-terms-in-series',50,' ')
      k=sqrt(2.d0*energy)
      eta=charge/k
      rho=abs(eta/.1d0)
      r=rho/k
      write(iout,2) energy, eta
      a=1
      b=a+ntrms+1
      need=wpadti(b+ntrms+1)
      call memory(need,pz,ngot,'z',0)
      call lrcoef(z(a),z(b),eta,l,ntrms,prn)       
      call asymp(fl,gl,dfl,dgl,z(a),z(b),rho,eta,l,ntrms,wron,
     1           fl0,gl0,dfl0,dgl0,wron0)
      dfl=k*dfl
      dgl=k*dgl
      wron=wron*k
      write(iout,3) r, rho
      write(iout,4) fl, dfl
      write(iout,5) gl, dgl
      write(iout,6) wron
      dfl0=dfl0*k
      dgl0=dgl0*k
      wron0=wron0*k
      write(iout,7) fl0, dfl0
      write(iout,8) gl0, dgl0
      write(iout,9) wron0
      call memory(-ngot,pz,idum,'z',idum)
      call chainx(0)               
      stop
 1    format(/,20x,'asymptotic coulomb package')      
 2    format(/,1x,'energy = ',e15.8,1x,'eta = ',e15.8)
 3    format(/,1x,'r = ',e15.8,1x,'rho =',e15.8)
 4    format(/,1x,'F  = ',e15.8,
     1       /,1x,'DF = ',e15.8)
 5    format(/,1x,'I  = ',e15.8,
     1       /,1x,'DI = ',e15.8)
 6    format(/,10x,'WRON = ',e15.8)
 7    format(/,1x,'F0  = ',e15.8,
     1       /,1x,'DF0 = ',e15.8)
 8    format(/,1x,'I0  = ',e15.8,
     1       /,1x,'DI0 = ',e15.8)
 9    format(/,10x,'WRON0 = ',e15.8)
      end
