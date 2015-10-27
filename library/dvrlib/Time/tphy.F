*deck tphy.f
c***begin prologue     tphy
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            physical hamiltonian.
c***                   
c***references         
c
c***routines called    
c***end prologue       tphy
      subroutine tphy(tin,twtin,ptin,dptin,ddptin,hin,
     1                t0out,tout,twtout,ptout,dptout,ddptout,hout,
     2                m,n,prn)
      implicit integer (a-z)
      real*8 tin, twtin, ptin, dptin, ddptin, hin
      real*8 t0out, tout, twtout, ptout, dptout, ddptout, hout
      character*80 title
      logical prn
      dimension tin(m), twtin(m), ptin(m,m), dptin(m,m), ddptin(m,m)
      dimension hin(m,m) 
      dimension tout(n), twtout(n), ptout(n,n), dptout(n,n)
      dimension ddptout(n,n), hout(n,n)
      dimension srf(2)
      common/io/inp, iout
      t0out=tin(1)
      call copy(tin(2),tout,n)
      call copy(twtin(2),twtout,n)      
      call mmove(ptin(2,2),ptout,n,n,m,n)
      call mmove(dptin(2,2),dptout,n,n,m,n)
      call mmove(ddptin(2,2),ddptout,n,n,m,n)
      call mmove(hin(2,2),hout,n,n,m,n)
      if(prn) then
         title='initial grid'
         call prntrm(title,tin,m,1,m,1,iout)
         title='initial weights'
         call prntrm(title,twtin,m,1,m,1,iout)
         title='initial functions'
         call prntrm(title,ptin,m,m,m,m,iout)
         title='initial first derivatives'
         call prntrm(title,dptin,m,m,m,m,iout)
         title='initial second derivatives'
         call prntrm(title,ddptin,m,m,m,m,iout)
         title='initial hamiltonian'
         call prntrm(title,hin,m,m,m,m,iout)      
         title='final grid'
         call prntrm(title,tout,n,1,n,1,iout)
         title='final weights'
         call prntrm(title,twtout,n,1,n,1,iout)
         title='final functions'
         call prntrm(title,ptout,n,n,n,n,iout)
         title='final first derivatives'
         call prntrm(title,dptout,n,n,n,n,iout)
         title='final second derivatives'
         call prntrm(title,ddptout,n,n,n,n,iout)
         title='final hamiltonian'
         call prntrm(title,hout,n,n,n,n,iout)
      endif
      return
      end       




