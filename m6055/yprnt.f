*deck @(#)yprnt.f	1.1 9/8/91
c***begin prologue     yprnt
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, spherical harmonic print
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print spherical harmonics
c***references         none
c
c***routines called    none
c***end prologue       yprnt
      subroutine yprnt(ylm,npnts,lmax,mumax,reg)
      implicit integer (a-z)
      character *80 title
      character *3 itoc
      character *8 colt, rowt
      common /io/ inp, iout
      real *8 ylm, rowv
      dimension ylm(npnts,0:lmax,0:2*mumax)
      dimension lval(0:100)
      do 1 l=0,lmax
         lval(l)=l
    1 continue
      rowt='point'
      colt='l'
      rowv=-99.d0
      write (iout,10) reg
      mcnt=0
      do 20 m=0,mumax
         iup=2
         if (m.eq.0) then
             iup=1
         endif
         do 40 mm=1,iup
            title='mu-'//itoc(m)//' comp-'//itoc(mm)
            call mprir(ylm(1,m,mcnt),rowv,lval(m),npnts,lmax+1-m,npnts,
     1                 lmax+1,title,rowt,colt,iout)
            mcnt=mcnt+1
   40    continue
   20 continue
   10 format(/,5x,'ylm for region',1x,i4)
      return
      end
