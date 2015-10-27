*deck @(#)dxcatm.f	5.1 11/6/94
      subroutine dxcatm(nocont,nat,ntypes,nbtype,nnp,
     $                  start,nbf,nobf,ops,nderiv,
     $                  grad,d,dkay,
     $                  ndmat,d2e,nd2e,nd1e)
c***begin prologue     dxcatm.f
c***date written       940513   (yymmdd)
c***revision date      11/6/94 
c***keywords
c***author             martin, richard and russo,thomas (lanl)
c***source             @(#)dxcatm.f	5.1 11/6/94
c
c***purpose
c
c***description
c
c   module to form the gradient of the exchange-correlation energy with
c   respect to the atomic centers.
c***references
c
c***routines called    (none)
c
c***end prologue       dxcatm.f
c
      implicit none
c     --- input variables -----
      integer nat,ntypes,nbtype
      integer nnp,nbf,nderiv,ndmat,nd1e,nd2e
c     --- input arrays (unmodified) ---
      character*(*) ops
      integer nocont(nat,ntypes)
      integer noprim(nat,ntypes)
      integer start(nat,ntypes)
      integer nobf(ntypes)
      real*8 dkay(nbf,nbf,3),d(nnp)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 grad(3,nat)
      real*8 d2e(nd2e)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer iatom,jatom,itype,jtype,jtypmx
      integer nconti,ncontj
c
      common/io/inp,iout
c
 1000 format (5x,'the exchange-correlation contribution to the scf ',
     $        'gradients:')
 1010 format (5x,'the exchange-correlation contribution to the scf ',
     $           'force constants:')
c
c     --- zero the space to accumulate derivative
c         exchange-correlation integrals
      if (nderiv.eq.2) then
         call rzero(d2e,nd2e)
         call rzero(grad,3*nat)
      else if (nderiv.eq.1) then
         call rzero(grad,3*nat)
      end if
c
c     --- form the gradient contributions ---
      do 9 iatom=1,nat
         do 8 jatom=1,iatom
            do 7 itype=1,nbtype
               if (noprim(iatom,itype).le.0) go to 7
               if (iatom.ne.jatom) then
                  jtypmx=nbtype
               else
                  jtypmx=itype
               end if
               do 6 jtype=1,jtypmx
                  if (noprim(jatom,jtype).le.0) go to 6
c
                  nconti=nocont(iatom,itype)
                  ncontj=nocont(jatom,jtype)
c
c                 --- dot with the density, and sum into array -----
                  if (nderiv.ge.1.and.iatom.ne.jatom) then
                     call d1xc(nat,nbf,nnp,nbtype,iatom,jatom,
     $                         itype,jtype,nconti,ncontj,nobf,
     $                         start,d,dkay,grad)
                  end if
c
c                 --- second derivative integrals, and sum into array -----
                  if (nderiv.ge.2.and.iatom.ne.jatom) then
                     call lnkerr('second derivatives not yet')
                  end if
c
c
    6          continue
    7       continue
    8    continue
    9 continue
c
c
      return
      end
