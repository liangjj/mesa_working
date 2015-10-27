*deck @(#)output.f	5.1  11/6/94
      subroutine output(vector,e,joutfg,idstor)
c
c  this routine writes out the results of the ci calculation.
c
      implicit real*8 (a-h,o-z)
      character*41 fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8
      character*24 dattim
      character*80 rtitle, ctitle, blabel
c
      common /headng/ rtitle, ctitle, blabel
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common/icntrl/ iciwrt, icipun, intwrt, ihamrd, ksym
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /c3/ ndet,nsef,idspcu,nroots,maxit,nw,maxesc,ianalz
      common /core/ ecore,thresh,ev,mask,nmask,limit1,limit2,intmxo
c
      dimension vector(nroots,*),e(*),joutfg(nbf),idstor(nw)
      dimension eev(7), emc(7)
c
c  field of max(29,nbf+15), then fields of 14
c
      mpr=max(29,nbf+15)
      icol=(132-mpr)/14
      msp=mpr-15
      mspl=(msp+1)/2
      mspr=msp-mspl
      write (fmt1,1400) mspl,mspr,icol
 1400 format('(//',i2,'x,10h    state ,',i2,'x,',i1,'i14)')
      write (fmt2,1500) mspl,mspr,icol
 1500 format('(/',i2,'x,15hvalence energy ,',i2,'x,',i1,'f14.8)')
      write (fmt3,1600) mspl,mspr,icol
 1600 format('(',i2,'x,15h core energy   ,',i2,'x,',i1,'f14.8)')
      write (fmt4,1700) mspl,mspr,icol
 1700 format('(',i2,'x,15h total energy  ,',i2,'x,',i1,'f14.8)')
      mspr=mspr-3
      write (fmt5,1800) mspl,mspr,icol
 1800 format('(',i2,'x,15h (e(i)-e(1)) ev,',i2,'x,',i1,'f14.5)')
      msp=mpr-28
      write (fmt6,1900) msp
 1900 format('(19h0 nsp configuration,',i2,'x,9hsym  ndgf/)')
      msp=mpr-nbf-12
      write (fmt7,2000) nbf,msp,icol
 2000 format('(i5,1x,',i3,'i1,i',i2,',i6,',i1,'f14.8)')
      write (fmt8,2100) mpr,icol
 2100 format('(i',i3,',',i1,'f14.8)')
c
      if(icipun.eq.0) go to 30
      open (icipun,file='punch',form='unformatted')
      write (icipun) rtitle, blabel, nroots, nsef
      do 10 i=1,nroots
      epc = e(i) + ecore
      write (icipun) i, epc
   10 write (icipun) (vector(i,j), j=1,nsef)
      close (icipun)
c
c     punch selected configurations if this is a selection run
c
c     open (iunts2,file='selcon')
c     write (iunts2) rtitle
c     write (iunts2) ntotfg
c     rewind iunts1
c     do 20 jjj=1,ntotfg
c     read (iunts1) jsym,nsefi,idstor
c     call unpck(idstor,joutfg,2,nbf)
c  20 write (iunts2) joutfg
c
   30 do 70 istart = 1,nroots,icol
      istop = min((istart - 1) + icol,nroots)
c
      rewind iunts1
c
      call btime(dattim,bsec)
      write (iw,1250) dattim, rtitle, blabel
 1250 format(1h1,23x,34hllnl/osu spin-orbit configuration-,
     1  19hinteraction program,5x,a24//
     2  5x,12hci output - ,a80//5x,22hintegral tape label - ,a80)
c
      do 40 i = istart,istop
      emc(i - (istart - 1)) = e(i) + ecore
   40 eev(i - (istart - 1)) = (e(i)-e(1))*ev
c
      write (iw,fmt1) (i, i=istart,istop)
      write (iw,fmt2) (e(i), i=istart,istop)
      write (iw,fmt3) (ecore, i=istart,istop)
      write (iw,fmt4) (emc(i - (istart - 1)), i=istart,istop)
      write (iw,fmt5) (eev(i - (istart - 1)), i=istart,istop)
      write (iw,fmt6)
c
      jstop = 0
c
      do 60  jjj = 1,ntotfg
c
      read (iunts1) jsym,nsefi,idstor
      call unpck(idstor,joutfg,2,nbf)
c
      lll = jstop + 1
      write (iw,fmt7) jjj, joutfg, jsym,
     1                lll, (vector(i,lll), i=istart,istop)
c
      jstart = jstop + 2
      jstop = jstop + nsefi
      do 50 lll = jstart,jstop
   50 write (iw,fmt8) lll, (vector(i,lll), i=istart,istop)
   60 continue
c
   70 continue
c
      return
c
      end
