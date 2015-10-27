*deck @(#)vprnt.f	1.1 9/8/91
c***begin prologue     vprnt
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, print potential
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print static potentials on grid
c***                   handles complex potentials
c***references         none
c
c***routines called    none
c***end prologue       vprnt
      subroutine vprnt(v,vc,npnts,ntri,reg)
      implicit integer (a-z)
      character *80 title
      character *3 itoc
      character *8 colt, rowt
      common /io/ inp, iout
      real *8 v, rowv
      complex *16 vc
      dimension v(npnts,ntri)
      title='potential matrices for region-'//itoc(reg)
      rowv=-99.d0
      colv=-99
      colt='state'
      rowt='point'
      call mprir(v,rowv,colv,npnts,ntri,npnts,ntri,title,rowt,colt,
     1           iout)
      return
      end
