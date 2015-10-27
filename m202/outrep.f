*deck @(#)outrep.f	5.1  11/6/94
      subroutine outrep(natoms,iprmut)
c***begin prologue     outrep.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           symmetry, i/o
c***author             martin, richard (lanl)
c***source             @(#)outrep.f	5.1   11/6/94
c***purpose            writes symmetry information to the rwf.
c***description
c     call outrep(natoms,iprmut)
c
c     writes /repcom/, /repnam/, and iprmut to the rwf.
c***references
c***routines called    iosys(io)
c***end prologue       outrep.f
      implicit none
c     --- input variables -----
      integer natoms
c     --- input arrays (unmodified) ---
      integer iprmut(natoms,10)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer nsymop,nreps,lblrep
      integer ltbl,lennam,lsym,lenprm
      real*8 symops,chrtbl
      character*4 repnam
c
      common/repcom/nsymop,nreps,lblrep(32),chrtbl(10,16),symops(9,10)
      common/repnam/repnam(32)
c
c     --- this routine writes out /repcom/,/repnam/, and iprmut onto the rwf.
      call iosys('write integer "number of irreducible '//
     $           'representations" to rwf',1,nreps,0,' ')
      call iosys('write integer "number of symmetry operations" '//
     $           'to rwf',1,nsymop,0,' ')
      lennam=nreps*len(repnam(1))
      call iosys('write character repnam to rwf',lennam,0,0,repnam)
      ltbl=160
      call iosys('write real chrtbl to rwf',ltbl,chrtbl(1,1),0,' ')
      lsym=90
      call iosys('write real symops to rwf',lsym,symops(1,1),0,' ')
      lenprm=natoms*nsymop
      call iosys('write integer "atom permutation array" to rwf',
     $           lenprm,iprmut,0,' ')
c
c
      return
      end
