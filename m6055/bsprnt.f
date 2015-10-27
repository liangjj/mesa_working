*deck @(#)bsprnt.f	1.1 9/8/91
c***begin prologue     bsprnt
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, bessel print
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print channel "bessel" functions
c***routines called    none
c***end prologue       bsprnt
      subroutine bsprnt(hp,hd,nchan,npnts,nlm,lch,maxlm,dimlm,
     1                  dimc,reg)
      implicit integer (a-z)
      character *80 title
      character *3 itoc
      character*8 rowt, colt
      common /io/ inp, iout
      complex *16 hp, hd
      real *8 rowv
      dimension hp(npnts,maxlm,nchan), hd(npnts,maxlm,nchan)
      dimension nlm(dimc), lch(dimlm,dimc)
      rowv=-99.d0
      rowt='point'
      colt='l'
      write (iout,10) reg
      do 20 ch1=1,nchan
         title='hp functions channel-'//itoc(ch1)
         call cmprir(hp(1,1,ch1),rowv,lch(1,ch1),npnts,nlm(ch1),npnts,
     1               maxlm,title,rowt,colt,iout)
         title='hd functions channel-'//itoc(ch1)
         call cmprir(hd(1,1,ch1),rowv,lch(1,ch1),npnts,nlm(ch1),npnts,
     1               maxlm,title,rowt,colt,iout)
   20 continue
      return
   10 format(/,5x,'free functions for region',1x,i4)
      end
