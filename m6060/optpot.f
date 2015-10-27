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
      subroutine optpot(ovplm,ovmlm,ovlmbm,ovpbm,ovmbm,vpp,vpm,vmm,
     1                  vbb,vbp,vbm,energy,eigval,nchan,maxlm,
     2                  nolam,nmo,ntri,north)
      implicit integer (a-z)
      real *8 energy 
      complex *16  eigval, ovplm, ovmlm, ovlmbm, ovpbm, ovmbm, cfac
      complex *16 vpp, vpm, vmm, vbb, vbp, vbm
      logical north
      dimension ovplm(1:maxlm,nolam,nchan), ovmlm(1:maxlm,nolam,nchan)
      dimension ovlmbm(nolam,nmo,nchan), ovpbm(1:maxlm,nmo,nchan)
      dimension ovmbm(1:maxlm,nmo,nchan), eigval(nolam)
      dimension vpp(1:maxlm,1:maxlm,ntri), vmm(1:maxlm,1:maxlm,ntri)
      dimension vpm(1:maxlm,1:maxlm,nchan,nchan), vbb(nmo,nmo,ntri)
      dimension vbp(nmo,1:maxlm,nchan,nchan)
      dimension vbm(nmo,1:maxlm,nchan,nchan)
      common /io/ inp, iout
c----------------------------------------------------------------------c
c          remove the bound state contribution from the integrals      c
c----------------------------------------------------------------------c
      if(.not.north) then
          do 10 i=1,nchan
             call cambct(ovplm(1,1,i),ovpbm(1,1,i),ovlmbm(1,1,i),maxlm,
     1                   nmo,nolam)
             call cambct(ovmlm(1,1,i),ovmbm(1,1,i),ovlmbm(1,1,i),maxlm,
     1                   nmo,nolam)
   10     continue
      endif
c----------------------------------------------------------------------c
c                  scale matrix elements by energy                     c
c                         denominator                                  c
c----------------------------------------------------------------------c
      do 20 lam=1,nolam
         cfac=1.d0/sqrt(.5d0*energy-eigval(lam))
         do 30 j=1,maxlm
            do 40 k=1,nchan
               ovplm(j,lam,k)=ovplm(j,lam,k)*cfac 
               ovmlm(j,lam,k)=ovmlm(j,lam,k)*cfac 
   40       continue
   30    continue
         do 50 j=1,nmo
            do 60 k=1,nchan
               ovlmbm(lam,j,k)=ovlmbm(lam,j,k)*cfac
   60       continue
   50    continue
   20 continue
c----------------------------------------------------------------------c
c                  form optical potential sums                         c
c----------------------------------------------------------------------c
      do 70 i=1,nchan
         do 80 j=1,nchan
            call cebct(vpm(1,1,i,j),ovplm(1,1,i),ovmlm(1,1,j),
     1                 maxlm,nolam,maxlm)
            call cebtct(vbp(1,1,i,j),ovlmbm(1,1,i),ovplm(1,1,j),nmo,
     1                  nolam,maxlm)
            call cebtct(vbm(1,1,i,j),ovlmbm(1,1,i),ovmlm(1,1,j),nmo,
     1                  nolam,maxlm)

   80    continue
   70 continue
      icnt=0
      do 100 i=1,nchan
         do 200 j=1,i     
            icnt=icnt+1
            call cebct(vpp(1,1,icnt),ovplm(1,1,i),ovplm(1,1,j),maxlm,
     1                 nolam,maxlm)
            call cebct(vmm(1,1,icnt),ovmlm(1,1,i),ovmlm(1,1,j),maxlm,
     1                 nolam,maxlm)
            call cebtc(vbb(1,1,icnt),ovlmbm(1,1,i),ovlmbm(1,1,j),nmo,
     1                 nolam,nmo)
  200    continue
  100 continue
      return
      end





