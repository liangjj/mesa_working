*deck @(#)mccga.f	1.2  7/30/91
      subroutine mccga(nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,lok,mix,g1,g2,
     $     nda1,lda1,nf35,nf36,
     $     cg,cr,icr,ncor,bufix,lbufso,
     $     rabcx,rabix,raibx,incor,sg,tg)
c
c***begin prologue     mccga
c***date written       871022   (yymmdd)
c***revision date      871116   (yymmdd)
c
c   16 november 1987   bhl at lanl
c    setting npqrs=1
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mccga.f	1.2   7/30/91
c
c***purpose
c
c***description
c
c
c-----------------------------------------------------------------------
c
c --- description     this program is the front-end for the
c                     programs that construct c*g.
c
c --- input
c
c     nsym            no of symmetries.
c     nbf(nsym)       number of basis functions.
c     nob(nsym)       number of orbitals.
c     nfob(nsym)      number of frozen orbitals
c     ncob(nsym)      number of core orbitals.
c     naob(nsym)      number of active orbitals.
c     cv(--)          orbitals.
c     npqrs           no of integral blocks.
c     locsym(nsym)    pointers into the arrays len and loc for the
c                     update vectors.
c     len(--)         no of nonzero componenets of update vectors.
c     lok(--)         location of update vectors in arrays mix and cg.
c     mix(--)         indices of nonzero components of update vectors.
c     g1(--)          one-electron density matrix elements.
c     g2(--)          two-electron density matrix elements.
c     nda1            direct access dataset containing the core
c                     fock-matrix.
c     lda1            address of the f-matrix in nda1.
c     nf35            fortran number of dataset containing the
c                     two electron integrals needed in the construction
c                     of the updated integrals.
c     nf36            fortran number of dataset containing the
c                     integrals needed to make core contributions
c                     to the updated integrals.
c
c --- output
c
c     cg(--)          c*g, packed according to len and mix
c
c --- working storage
c
c     cr(ncor)
c
c-----------------------------------------------------------------------
c
c
c***references
c
c***routines called    (none)
c
c***end prologue       mccga
c
      implicit real*8(a-h,o-z)
      integer wpadti
cc
cc
      real*8 bufix(lbufso),sg(*),tg(*)
      real*8 rabcx(*),rabix(*),raibx(*)
      dimension nbf(2),nob(2),ncob(2),naob(2),cv(2)
      dimension locsym(2),len(2),lok(2),mix(2),g1(2),g2(2)
      dimension cg(2),nfob(2)
      dimension cr(2),icr(2)
c
      common /io/ inp,iout
c-----------------------------------c
c      make storage allocation      c
c-----------------------------------c
      npqrs=1
      npq = 0
      maxoo = 0
      maxo = 0
      maxb = 0
      ncgx = 0
      do 100 l = 1, nsym
         maxoo = max(ncob(l),naob(l),maxoo)
         maxo = max(nob(l),maxo)
         maxb = max(nbf(l),maxb)
         npq = npq + nbf(l) * (nbf(l) + 1) / 2
         ncgx = ncgx + (ncob(l) + naob(l)) * nob(l)
 100  continue
c
      lr = 1
      lcgx = lr + maxoo * maxo
      ipqrs = wpadti(lcgx + ncgx)
      lfpq = iadtwp(ipqrs + npqrs)
      ltpi = lfpq + npq
      last = wpadti(ltpi + maxb)
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(last,cr,ncor,'mccga',0)
c
cpsdyn      if (last .le. ncor) go to 500
cpsdyn      write (iout,9001) last, ncor
cpsdyn 9001 format(//' ****** mccga  insufficient core (need,given0 ',2i10)
cpsdyn      stop
c
      call rzero(cr(1),iadtwp(last))
c
      call mccgb(nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,lok,mix,g1,g2,
     $     nda1,lda1,nf35,nf36,cg,cr(lcgx),
     $     icr(ipqrs),cr(lr),cr(lfpq),cr(ltpi),bufix,lbufso,
     $     rabcx,rabix,raibx,incor,sg,tg)
      return
      end
