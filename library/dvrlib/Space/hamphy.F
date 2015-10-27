*deck hamphy.f
c***begin prologue     hamphy
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            physical hamiltonian and potential.
c***                   
c***references         
c
c***routines called    
c***end prologue       hamphy
      subroutine hamphy(qin,wtin,pin,dpin,ddpin,hamin,vin,
     1                  q0,qout,wtout,pout,dpout,ddpout,hamout,vout,
     2                  srf,m,n,start)
      implicit integer (a-z)
      real*8 qin, wtin, pin, dpin, ddpin, hamin, vin
      real*8 q0, qout, wtout, pout, dpout, ddpout
      real*8 hamout, vout, srf
      character*80 title
      dimension qin(m), wtin(m), pin(m,m), dpin(m,m), ddpin(m,m)
      dimension hamin(m,m), vin(m) 
      dimension qout(n), wtout(n), pout(n,n), dpout(n,n), ddpout(n,n)
      dimension hamout(n,n), vout(n)
      dimension srf(2)
      common/io/inp, iout
      q0=qin(1)
      call copy(qin(start),qout,n)
      call copy(wtin(start),wtout,n)      
      call mmove(pin(start,start),pout,n,n,m,n)
      call mmove(dpin(start,start),dpout,n,n,m,n)
      call mmove(ddpin(start,start),ddpout,n,n,m,n)
      call mmove(hamin(start,start),hamout,n,n,m,n)
      call copy(vin(start),vout,n)      
      srf(1)=pin(1,1)
      srf(2)=pin(m,m)
      return
      end       



