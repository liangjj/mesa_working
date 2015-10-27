*deck extphs
      subroutine extphs(f,r,l,eta,ak,bk,ntrms,prnt)
c***begin prologue     extphs
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            extract coulomb phase shift from asymptotic
c***                   wavefunction
c***description
c***references
c
c***routines called
c
c***end prologue       extphs
c
      implicit integer (a-z)
      real*8 f, r, rinv, eta, fl, gl, dfl, dgl, ak, bk, wron
      logical prnt
      dimension ak(*), bk(*)
      common/io/inp,iout
c**********************************************************************c
c                use the asymptotic routine at one point to compare    c
c                the numerically calculated phase shift obtained by    c
c                outward integration with the "exact" value.           c
c**********************************************************************c
      rinv=1.d0/r
      npt=1
      call asymp(fl,gl,dfl,dgl,ak,bk,r,rinv,eta,l,npt,ntrms,wron,prnt)
      write(iout,1)
      write(iout,2) r, f, fl
    1 format(/,5x,'approximate and exact asymptotic regular coulomb func
     1tion')   
    2 format(/,10x,'r = ',e15.8,1x,'fun = ',e15.8,1x,'exact = = ',e15.8)
      return
      end   





