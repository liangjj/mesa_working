*deck @(#)bndprn.f	1.1 9/8/91
c***begin prologue     bndprn
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, print
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print bound orbitals
c***references         none
c
c***routines called    none
c***end prologue       bndprn
      subroutine bndprn(orbs,npnts,ncon,reg)
      implicit integer (a-z)
      common /io/ inp, iout
      character *80 title
      character *8 rowt, colt
      character *3 itoc
      real *8 orbs, rowv
      dimension orbs(npnts,ncon)
      title='bound orbitals region-'//itoc(reg)
      rowv=-99.d0
      colv=-99
      colt='con ao'
      rowt='point'
      call mprir(orbs,rowv,colv,npnts,ncon,npnts,ncon,title,rowt,colt,
     1            iout)
      return
      end
