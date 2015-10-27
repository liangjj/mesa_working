*deck @(#)optpot.f
c***begin prologue     optpot
c***date written       xxxxxx   (yymmdd)
c***revision date      920409   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (nsf)
c***source             m6060
c***purpose            optical potential construction
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       optpot
      subroutine optpot(ovlmbm,vbb,energy,eigval,nolam,
     1                  nmo,prnt)
      implicit integer (a-z)
      real *8 energy 
      complex *16  eigval, ovlmbm, vbb, cfac
      character *80 title
      logical prnt
      dimension ovlmbm(nolam,nmo), eigval(nolam)
      dimension vbb(nmo,nmo)
      common /io/ inp, iout
c----------------------------------------------------------------------c
c                  scale matrix elements by energy                     c
c                         denominator                                  c
c----------------------------------------------------------------------c
      do 10 lam=1,nolam
         cfac=1.d0/sqrt(energy-eigval(lam))
         do 20 j=1,nmo
            ovlmbm(lam,j)=ovlmbm(lam,j)*cfac
   20    continue
   10 continue
c----------------------------------------------------------------------c
c                  form optical potential sums                         c
c----------------------------------------------------------------------c
      call cebtc(vbb,ovlmbm,ovlmbm,nmo,nolam,nmo)
      if (prnt) then
          title='optical potential'
          call prntcm(title,vbb,nmo,nmo,nmo,nmo,iout)
      endif
      return
      end





