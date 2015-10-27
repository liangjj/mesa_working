*deck @(#)mcgvca.f	1.2  7/30/91
      subroutine mcgvca(nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,lok,mix,b,nint,
     $     nda1,lda1,nf35,nf36,nf37,ltrb,
     $     cr,icr,ncor,intape,buf,lbufso,
     $     rabcx,rabix,raibx,incor,sg)
c
c***begin prologue     mcgva
c***date written       871022   (yymmdd)
c***revision date      871116   (yymmdd)
c
c   16 november 1987   bhl at lanl
c    setting npqrs=1
c
c   10 november 1987   bhl at brl
c      changing core allocation to number of active, not occupied,
c      orbitals. ('naot').
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcgvca.f	1.2   7/30/91
c
c***purpose
c
c***description
c
c-----------------------------------------------------------------------
c
c --- description     this routine is the front-end for the
c                     routines that calculates the updated one-
c                     and two-electron integrals.
c
c --- input
c
c     nsym            no of symmetries.
c     nbf(nsym)       number of basis functions.
c     nob(nsym)       number of orbitals.
c     ncob(nsym)      number of core orbitals.
c     naob(nsym)      number of active orbitals.
c     cv(--)          orbitals.
c     npqrs           no of integral blocks.
c     nint            no of integrals.
c     locsym(nsym)    pointers into the arrays len and loc for the
c                     update vectors.
c     len(--)         no of nonzero componenets of update vectors.
c     lok(--)         location of update vectors in arrays mix and cm.
c     mix(--)         indices of nonzero components of update vectors.
c     b(--)           nonzero components of update vectors.
c     nda1            direct access dataset containing the core
c                     fock-matrix.
c     lda1            address of the f-matrix in nda1.
c     nf35            fortran number of dataset containing the
c                     two electron integrals needed in the construction
c                     of the updated integrals.
c     nf36            fortran number of dataset containing the
c                     integrals needed to make core contributions
c                     to the updated integrals.
c     nf37            fortran number of dataset containing the updated
c                     integrals.
c     ltrb            record length for nf37.
c
c --- working storage
c
c     cr(ncor)
c
c-----------------------------------------------------------------------
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cc
cc
      integer wpadti
      logical debug
      real*8 buf(lbufso),sg(*)
      real*8 rabcx(*),rabix(*),raibx(*)
      dimension nbf(2),nob(2),ncob(2),naob(2),cv(2)
      dimension locsym(2),len(2),lok(2),mix(2),b(2)
      dimension cr(2),icr(2),nfob(2)
c
      data debug/.false./
      common /io/ inp,iout
c
c
      if (debug) then
         write (iout,9001)
 9001    format(/'0************ updated integral section',
     1           ' ************')
      end if
c
c     call jtime(it1)
c-----------------------------------c
c      make storage allocation      c
c-----------------------------------c
      npqrs=1
      npq = 0
      maxoo = 0
      maxo = 0
      maxb = 0
      ncm = 0
      nobt = 0
      naot=0
      do 100 l = 1, nsym
         maxoo = max(ncob(l),naob(l),maxoo)
         maxo = max(nob(l),maxo)
         maxb = max(nbf(l),maxb)
         npq = npq + nbf(l) * (nbf(l) + 1) / 2
         ncm = ncm + (ncob(l) + naob(l)) * nob(l)
         nobt = nobt + (naob(l) + ncob(l))
         naot=naot+naob(l)
 100  continue
c
cps      nij=(nobt*(nobt+1))/2
cps      nijs=nij*nij
cps      nijt= nij*nobt*nobt+nij
      nij=(naot*(naot+1))/2
      nijs=nij*nij
      nijt= nij*naot*naot+nij
c
      lg1=1
      lg2=lg1+nij
      lr = lg2 + nijs
      lg = lr + maxoo * maxo
c      ipqrs = lg + nint
      ipqrs = wpadti(lg + max(nint,nijt))
      lfpq = iadtwp(ipqrs + npqrs)
      ltpi = lfpq + npq
      lcm = ltpi + maxb
      last = wpadti(lcm + ncm)
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(last,cr,ncor,'mcgvca',0)
c
cpsdyn      if (last .le. ncor) go to 500
cpsdyn      write (iout,9901) last, ncor
cpsdyn 9901 format(//' ****** mcgvca  insufficient core (need,given0 ',2i10)
cpsdyn      stop
c
      call rzero(cr(1),iadtwp(last))
      ncor2=iadtwp(ncor)-lr
c
      call mcgvcb(nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,lok,mix,cr(lcm),
     $     b,nda1,lda1,nf35,nf36,nf37,ltrb,cr(lg),
     $     icr(ipqrs),cr(lr),cr(lfpq),cr(ltpi),
     $     cr(lg1),cr(lg2),ncor2,intape,buf,lbufso,
     $     rabcx,rabix,raibx,incor,sg)
c
c     call jtime(it2)
c     tim  = (it2 - it1) / 100.d0
c     write (iout,9099) tim
c9099 format('0   ****** tim  in this section ',10x,f10.2,' secs')
c
      return
      end
