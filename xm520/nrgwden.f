*deck @(#)nrgwden.f	5.1 11/6/94
      subroutine nrgwden(c,eigval,d,nbf,nnp,nocc)
c***begin prologue     nrgwden.f
c***date written       840820  (yymmdd)
c***revision date      11/6/94
c
c***keywords           
c***author             russo, thomas (lanl)
c***source             nrgwden.f
c***purpose            make energy weighted density matrix
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       nrgwden.f
      implicit none
c     --- input variables -----
      integer nbf,nnp,nocc
c     --- input arrays (unmodified) ---
      real*8 c(nbf,nbf),eigval(nbf)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
      real*8 d(nnp)
c     --- local variables ---
      integer i,ij,mu,nu
      integer inp,iout
c
      common/io/inp,iout
c
c
      
      ij=0
      do 1 mu=1,nbf
         do 2 nu=1,mu
            ij=ij+1
            d(ij)=0.0
            do 3 i=1,nocc
               d(ij)=d(ij)+c(mu,i)*c(nu,i)*eigval(i)
 3          continue 
 2       continue 
 1    continue 

      call iosys('write real "dft nrg weighted density matrix" on rwf',
     $     nnp,d,0,' ')

      return
      end
