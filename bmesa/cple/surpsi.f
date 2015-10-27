*deck surpsi
c***begin prologue     surpsi
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           surface, functions
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate surface functions in r-matrix basis
c***description        
c***references       
c
c***routines called
c***end prologue       surpsi
      subroutine surpsi(ham,psi,dum,nclosd,nopen,nc,lenclo,n,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 ham, psi, dum
      character*80 title
      logical prnt
      dimension ham(n,n), psi(n), dum(n,nc), nopen(nc), nclosd(nc)
      if (prnt) then
          title='transformation matrix'
          call prntrm(title,ham,n,n,n,n,iout)
      endif
      do 10 i=1,n
         klwr=lenclo
         do 20 j=1,nc
            dum(i,j)=0.d0
            klwr=klwr+1
            kupr=klwr+nopen(j)-1
            do 30 k=klwr,kupr
               dum(i,j)=dum(i,j)+ham(k,i)*psi(k)
   30       continue
            klwr=kupr
   20    continue      
   10 continue
      call iosys ('write real "hamiltonian eigenfunction surface '//
     1            'projections" to lamdat',n*nc,dum,0,' ')        
      if (prnt) then
          title='surface projections'
          call prntrm(title,dum,n,nc,n,nc,iout)
      endif
      return
      end
