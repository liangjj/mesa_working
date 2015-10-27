*deck frmopt.f
c***begin prologue     frmopt
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           optical potential
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form the optical potential from the input static-exchange
c***                   piece plus the energy dependent, non-local dynamical 
c***                   potential.
c***description                           
c***                   
c***references         
c
c***routines called    
c***end prologue       frmopt
      subroutine frmopt(opt,xqp,hqp,header,energy,n,npvec,prnt)
      implicit integer (a-z)
      real*8 xqp, hqp, energy
      character*80 title
      character*(*) header
      character*16 fptoc
      logical prnt
      dimension xqp(n,npvec)
      dimension hqp(n,npvec)
      dimension opt(npvec,npvec)
      common/io/inp, iout 
c
c     read in the "static-exchange" part
c
      call iosys('read real '//header//' from hamiltonian '//
     #           'without rewinding',npvec*npvec,opt,0,' ')
      call apbtc(opt,hqp,xqp,npvec,n,npvec)
      if(prnt) then
         title='static-exchange plus optical potential energy = '
     #          //fptoc(energy)
         call prntrm(title,opt,npvec,npvec,npvec,npvec,iout)
      endif
c
c     overwrite the "static-exchange" hamiltonian
c     with the full matrix
c 
      call iosys('write real '//header//' to hamiltonian '//
     #           'without rewinding',npvec*npvec,opt,0,' ')                  
      return
      end       
