*deck @(#)grppt.f	1.1  11/20/92
      subroutine grppt(nocont,nobf,natoms,ntypes,a,ops)
c
c***begin prologue     grppt
c***date written       871116   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci gradients, mcscf gradients
c***author             saxe, paul (lanl)
c***source             @(#)grppt.f	1.1   11/20/92
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
c***end prologue       grppt
c
      implicit integer (a-z)
c
      character*(*) ops
      integer nocont(natoms,ntypes)
      integer nobf(ntypes)
      integer a(*)
c
c     ----- count up the number of momentum groups and basis functions
c
      ngrp=0
      nbf=0
      do 2 grp=1,ntypes
         do 1 atom=1,natoms
            if (nocont(atom,grp).gt.0) then
               n=nocont(atom,grp)
               ngrp=ngrp+1
               nbf=nbf+n*nobf(grp)
            end if
 1       continue
 2    continue
c
c     ----- check on the number of basis functions -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbasis,0,' ')
c
      if (nbf.ne.nbasis) call lnkerr('mortal error with numbers of '//
     $     'basis functions')
c
c     ----- allocate core, again -----
c
      nnpgrp=ngrp*(ngrp+1)/2
c
      bftgrp=1
      bftcmp=bftgrp+nbf
      grpsiz=bftcmp+nbf
      gptij=grpsiz+ngrp
      gptkl=gptij+nnpgrp
      gklsiz=gptkl+nnpgrp
      need=gklsiz+nnpgrp
c
      call getscm(need,a,ngot,'group pointers',0)
c
      call grppt1(nocont,nobf,natoms,ntypes,a(bftgrp),a(bftcmp),
     $     a(grpsiz),a(gptij),a(gptkl),a(gklsiz),ngrp,nnpgrp,nbf,
     $     ops)
c
c
      return
      end
