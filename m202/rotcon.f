*deck @(#)rotcon.f	5.1  11/6/94
      subroutine rotcon(iout,natoms,ian,c,ptgrp,bb,prtrot,atmass)
c***begin prologue     rotcon.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           rotational constants, isotopes
c***author             binkley, et al.(g82)
c***source             @(#)rotcon.f	5.1   11/6/94
c***purpose            computes and prints rotational constants for the
c                      principal isotopes.
c***description
c     call rotcon(iout,natoms,ian,c,ptgrp,bb,prtrot)
c***references
c***routines called    filmas(m202), mofi(m202), symnum(m202)
c***end prologue       rotcon.f
      implicit none
c     --- input variables -----
      integer iout,natoms
      logical prtrot
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      character*(*) bb(*)
      character*(*) ptgrp
      real*8 c(3,natoms)
      real*8 atmass(natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
      logical linear
      real*8 pmom(3)
      real*8 eigvec(9)
      real*8 zero,one,four,eight,agiga
      real*8 boltz,planck,tomet,tokg,pi
      real*8 pipi,con,rca,rcb,rcc,rtemp1,rtemp2,rtemp3
c
      data zero/0.0d+00/,one/1.0d+00/,four/4.0d+00/,eight/8.0d+00/
      data agiga/1.0d+09/
      save zero,one,four,eight,agiga
c
 1000 format(' rotational constants (ghz):',3f15.7)
 1010 format(' isotopes:')
 1020 format(5x,16a8)
c
c
c     --- get necessary constants.
      if(natoms.lt.2) return
      call iosys('read real boltzmann from rwf',1,boltz,0,' ')
      call iosys('read real planck from rwf',1,planck,0,' ')
      call iosys('read real m/bohr from rwf',1,tomet,0,' ')
      call iosys('read real kg/amu from rwf',1,tokg,0,' ')
      call iosys('read real pi from rwf',1,pi,0,' ')
      pipi = pi * pi
      con = planck / (boltz*eight*pipi)
      con = (con / tokg)  *  (planck / (tomet*tomet))
c
c     get masses of most common isotope and moments of inertia.
c
      call filmas(0,iout,ian,natoms,.true.,atmass,bb)
      call mofi(natoms,ian,c,atmass,pmom,eigvec)
c
c     compute rotational constants in ghz.  beware of linear molecules!
c
      rtemp1 = zero
      rtemp2 = zero
      rtemp3 = zero
      linear=ptgrp(2:2).eq.'*'
      if(.not.linear.and.pmom(1).ne.zero) rtemp1 = con / pmom(1)
      if(pmom(2).ne.zero) rtemp2 = con / pmom(2)
      if(pmom(3).ne.zero) rtemp3 = con / pmom(3)
      rca = boltz * rtemp1 / agiga / planck
      rcb = boltz * rtemp2 / agiga / planck
      rcc = boltz * rtemp3 / agiga / planck
      if(prtrot) then
         write(iout,1000) rca,rcb,rcc
         write(iout,1010)
         write(iout,1020) (bb(i),i=1,natoms)
      endif
c
c
      return
      end
