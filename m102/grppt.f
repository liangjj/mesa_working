*deck @(#)grppt.f	5.1  11/6/94
      subroutine grppt(nocont,nobf,natoms,ntypes,ops)
c***begin prologue     grppt.f
c***date written       871116   (yymmdd)
c***revision date      11/6/94
c
c***keywords           ci gradients, mcscf gradients
c***author             saxe, paul (lanl)
c***source             @(#)grppt.f	5.1   11/6/94
c
c***purpose            to create pointer arrays for group ordering the
c  two-particle density matrix for the derivative programs.
c
c***description
c
c   the address for a density-matrix element g(i,j,k,l) is
c
c   address = gptij([bftgrp(i),bftgrp(j)])
c           + gptkl([bftgrp(k),bftgrp(l)])
c                *gklsiz([bftgrp(i),bftgrp(j)])
c           + bftcmp(i)
c              + grpsiz(bftgrp(i)) * (bftcmp(j)
c                 + grpsiz(bftgrp(j)) * (bftcmp(k)
c                    + grpsiz(bftgrp(k)) * bftcmp(l)))
c
c***references
c
c***routines called    (none)
c
c***end prologue       grppt.f
      implicit none
c     --- input variables -----
      integer inp, iout
      integer natoms,ntypes
c     --- input arrays (unmodified) ---
      integer nocont(natoms,ntypes)
      integer nobf(ntypes)
      character*(*) ops
c     --- input arrays (scratch) ---
      integer a
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer ngrp,nbf,grp,atom,n
      integer nbasis,bftgrp,gptkl,gptij,nnpgrp
      integer bftcmp,grpsiz,gklsiz
      integer top,ngot,mxcore, idum
      common/io/inp, iout
      pointer(p,a(1))
c
c     --- count up the number of momentum groups and basis functions
      ngrp=0
      nbf=0
      do 20 grp=1,ntypes
         do 10 atom=1,natoms
            if (nocont(atom,grp).gt.0) then
               n=nocont(atom,grp)
               ngrp=ngrp+1
               nbf=nbf+n*nobf(grp)
            end if
 10      continue
 20   continue
c
c     --- check on the number of basis functions ---
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbasis,0,' ')
c
      if (nbf.ne.nbasis) call lnkerr('mortal error with numbers of '//
     $     'basis functions')
c
c     --- allocate core, again ---
      nnpgrp=ngrp*(ngrp+1)/2
c
      bftgrp=1
      bftcmp=bftgrp+nbf
      grpsiz=bftcmp+nbf
      gptij=grpsiz+ngrp
      gptkl=gptij+nnpgrp
      gklsiz=gptkl+nnpgrp
      top=gklsiz+nnpgrp
      call iosys('read integer mxcore from rwf',1,mxcore,0,' ')
      if (top.gt.mxcore) then
         write (iout,1) top,mxcore
         call lnkerr('m102 core size problem')
      else
         call getmem(top,p,ngot,'grppt',0)
      end if
c
c      call getscm(top,a,ngot,'group pointers',0)
c
      call grppt1(nocont,nobf,natoms,ntypes,a(bftgrp),a(bftcmp),
     $            a(grpsiz),a(gptij),a(gptkl),a(gklsiz),ngrp,nnpgrp,nbf,
     $            ops)
      call getmem(-ngot,p,idum,'grppt',idum)
c
c
      return
    1    format (//,' $$$$$ m102 need ',i8,' words of core, have ',
     $           i8,//)
      end
