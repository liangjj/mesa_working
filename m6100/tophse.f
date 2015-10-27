*deck @(#)tophse.f
c***begin prologue     tophse
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           basis, link 6080
c***author             schneider, barry (lanl)
c***source             m6090
c***purpose            to form the kohn matrix, solve the equations and extract
c***                   the t matrix and phase,
c***description        matrix elements of ( h - E ) are assembled fron
c***                   bound-bound, bound-free and free-free integrals.
c***                   the same is done for the right hand side. the linear
c***                   equations are solved and the scattering information
c***                   extracted.
c***references       
c
c***routines called
c***end prologue       tophse
      subroutine tophse(hbbt,hfbt,hfft,energy,mtrx,rhs,ipvt,nbfn,prnt)
      implicit integer (a-z)
      real *8 hbbt, energy, phase
      complex *16 hfft, hfbt, mtrx, rhs, eye, tmat, tndel
      logical prnt
      character *80 title
      dimension hbbt(nbfn,nbfn), hfbt(2,nbfn), hfft(2,2)
      dimension mtrx(nbfn+1,nbfn+1), rhs(nbfn+1), ipvt(nbfn+1)
      common /io/ inp, iout
c
      eye=cmplx(0.d0,1.d0)
      call czero(mtrx,(nbfn+1)*(nbfn+1))
      call czero(rhs,nbfn+1)
      do 10 i=1,nbfn
         do 20 j=1,nbfn
            mtrx(i,j)=hbbt(i,j)
   20    continue
   10 continue   
      do 30 i=1,nbfn
         mtrx(i,nbfn+1)=hfbt(2,i)
         mtrx(nbfn+1,i)=hfbt(2,i)
         rhs(i)=-hfbt(1,i)
   30 continue
      mtrx(nbfn+1,nbfn+1)=hfft(2,2)       
      rhs(nbfn+1)=-hfft(2,1)
      if (prnt) then
          title='complex kohn matrix'
          call prntcmn(title,mtrx,nbfn+1,nbfn+1,nbfn+1,nbfn+1,iout,'e')
          title='complex kohn rhs'
          call prntcmn(title,rhs,nbfn+1,1,nbfn+1,1,iout,'e')
      endif
      call cgefa(mtrx,nbfn+1,nbfn+1,ipvt,info)
      call cgesl(mtrx,nbfn+1,nbfn+1,ipvt,rhs,0)
      write(iout,100)
      tmat=rhs(nbfn+1)
      tndel=tmat/(1.d0+eye*tmat)
      write(iout,200) tmat, tndel
      phase=atan(real(tndel))
      write(iout,300) phase
      tmat=hfft(1,1) 
      do 40 i=1,nbfn
         tmat=tmat+hfbt(1,i)*rhs(i)
   40 continue
      tmat=tmat+hfft(2,1)*rhs(nbfn+1)
      tmat=-2.d0*tmat
      write(iout,400)     
      tndel=tmat/(1.d0+eye*tmat)
      write(iout,200) tmat, tndel
      phase=atan(real(tndel))
      write(iout,300) phase
      return
  100 format(//,25x,'kohn uncorrected results') 
  200 format(//,10x,'t matrix =(',e15.8,',',e15.8,')',/,10x,
     1              'tangent phase shift = (',e15.8,',',e15.8,')')
  300 format(//,10x,'phase shift = ',e15.8) 
  400 format(//,25x,'kohn corrected results') 
      end






