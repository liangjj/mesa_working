*deck diag
      subroutine diag(ham,bfbox,dbfbox,hfbox,dhfbox,eig,scr,
     1                ncon,type)
      implicit integer (a-z)
      real*8 ham, bfbox, dbfbox, hfbox, dhfbox, eig, scr
      character*(*) type
      dimension ham(ncon,ncon), bfbox(ncon), dbfbox(ncon)
      dimension eig(ncon), scr(*), hfbox(ncon), dhfbox(ncon)
      common /io/ inp, iout
      call tred2(ncon,ncon,ham,eig,scr,ham)
      call tql2(ncon,ncon,eig,scr,ham,ierr)
      if (type.eq.'s-wave') then
          write(iout,1) (eig(i),i=1,ncon)  
      else
          write(iout,2) (eig(i),i=1,ncon)           
      endif
      call ebtc(hfbox,ham,bfbox,ncon,ncon,1)
      write(iout,3) (hfbox(i),i=1,ncon)
      call ebtc(dhfbox,ham,dbfbox,ncon,ncon,1)
      write(iout,4) (dhfbox(i),i=1,ncon)      
 1    format(/,1x,'r matrix eigenvalues for s wave manifold',
     1             (/,5e15.8))
 2    format(/,1x,'r matrix eigenvalues for p wave manifold',
     1             (/,5e15.8))
 3    format(/,1x,'hamiltonian eigenfunction values at box',
     1             (/,5e15.8))
 4    format(/,1x,'derivative of hamiltonian eigenfunction values at box
     1',         (/,5e15.8))
      return
      end
