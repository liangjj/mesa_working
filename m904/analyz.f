*deck @(#)analyz.f	5.1  11/6/94
      subroutine analyz(vector,e,row,h,dele,joutfg,idstor,istsef,
     2  indx,jfigs,jsyms)
c
c  calculate the energy contributions for each configuration for
c  each root and write out a summary of the ci calculation.
c
      implicit real*8 (a-h,o-z)
      character*41 fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8
      character*24 dattim
      character*80 rtitle, ctitle, blabel
c
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /headng/ rtitle, ctitle, blabel
      common/icntrl/ iciwrt, icipun, intwrt, ihamrd, ksym
      common /core/ ecore,thresh,ev,mask,nmask,limit1,limit2,intmxo
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /c3/ ndet,nsef,idspcu,nroots,maxit,nw,maxesc,ianalz
c
      dimension vector(nroots,*),e(*),row(*),h(*),dele(*),joutfg(nbf),
     1  idstor(nw),istsef(*),indx(*),jfigs(nw,*),jsyms(*)
c
*mdc*if cray
*      equivalence (rnext,next)
*mdc*else
      dimension nexta(2)
      equivalence (nexta(1),rnext), (nexta(2),next)
*mdc*endif
c
      parameter (deflt=1.0d-3)
c
c  field of max(29,nbf+15), then fields of 14
c
      mpr=max(29,nbf+15)
      msp=mpr-15
      mspl=(msp+1)/2
      mspr=msp-mspl
      write (fmt1,1400) mspl,mspr
 1400 format('(//',i2,'x,10h    state ,',i2,'x,i14)')
      write (fmt2,1500) mspl,mspr
 1500 format('(/',i2,'x,15hvalence energy ,',i2,'x,f14.8)')
      write (fmt3,1600) mspl,mspr
 1600 format('(',i2,'x,15h core energy   ,',i2,'x,f14.8)')
      write (fmt4,1700) mspl,mspr
 1700 format('(',i2,'x,15h total energy  ,',i2,'x,f14.8)')
      mspr=mspr-3
      write (fmt5,1800) mspl,mspr
 1800 format('(',i2,'x,15h (e(i)-e(1)) ev,',i2,'x,f14.5)')
      msp=mpr-28
      write (fmt6,1900) msp
 1900 format('(19h0 nsp configuration,',i2,'x,9hsym  ndgf/)')
      msp=mpr-nbf-12
      write (fmt7,2000) nbf,msp
 2000 format('(i5,1x,',i3,'i1,i',i2,',i6,3f14.8)')
      write (fmt8,2100) mpr
 2100 format('(i',i3,',2f14.8)')
c
c  read in the threshold for the energy selection
c
      thresh=deflt
      ipunch=0
      if (ianalz.ne.0) then
        read (5,'(d10.2,i5)') thresh,ipunch
        if((thresh.lt.1.0d-6).or.(thresh.gt.1.0d-1)) thresh = deflt
        write (iw,2300) thresh,ipunch
 2300   format (//5x,'configurations listed with energy ',
     1    'contributions greater than',f15.7//5x,'ipunch =',i10)
      endif
c
c  calculate the istsef array - contains the starting points for the
c  double-group function expansions
c
      rewind iunts1
      istsef(1) = 1
      do 110 ifig=1,ntotfg
      read (iunts1) jsym, nsefi
  110 istsef(ifig+1) = istsef(ifig) + nsefi
c
c  store the diagonal elements of the hamiltonian matrix in h
c
      rewind iunt2a
      read (iunt2a)
      next = 2
      do 120 isef=1,nsef
      call hin(row,next)
      h(isef) = row(next-1)
  120 rnext = row(next)
c
c  analyze each of the roots
c
      call btime(dattim,bsec)
      write (iw,2400) dattim, rtitle, blabel
 2400 format(1h1,29x,34hllnl/osu spin-orbit configuration-,
     1  19hinteraction program,5x,a24//
     2  5x,24hanalysis of ci output - ,a80//
     3  5x,22hintegral tape label - ,a80)
      do 290 iroot=1,nroots
c
c  write out the labels for the ci calculation
c
      emc = e(iroot) + ecore
      eev = (e(iroot)-e(1))*ev
      write (iw,fmt1) iroot
      write (iw,fmt2) e(iroot)
      write (iw,fmt3) ecore
      write (iw,fmt4) emc
      write (iw,fmt5) eev
      write (iw,fmt6)
c
c     return if root has not converged
c     if (abs(e(iroot)).lt.1.0d-5) return
c     if (ipunch.ne.0) write (ipunch,2500) iroot
c2500 format(5x,'root no.',i5,' dominant configurations')
c
c  calculate the energy contributions for each spatial configuration
c
      do 130 ifig=1,ntotfg
      dele(ifig) = 0.0d0
      istart = istsef(ifig)
      istop = istsef(ifig+1) - 1
      do 130 isef=istart,istop
      vect2 = vector(iroot,isef)**2
      if(vect2.eq.1.0d0) then
        dele(ifig) = e(iroot)
      else
        dele(ifig) = dele(ifig)-vect2*(h(isef)-e(iroot))/(1.0d0-vect2)
      endif
  130 continue
c
c  sort the energy contributions according to magnitude
c
      call sort(dele,indx,ntotfg)
c
c  write out all configurations with energy contributions greater
c  than thresh
c
      iabort = 0
      do 150 ifig=1,ntotfg
      if(abs(dele(ifig)).ge.thresh) go to 150
      nafigs = ifig - 1
      go to 160
  150 continue
      nafigs = ntotfg
  160 if(nafigs.gt.maxesc) then
        iabort = 1
        nafigs = maxesc
      endif
c
      rewind iunts1
      do 200 ifig=1,ntotfg
      read (iunts1) jsym, nsefi, idstor
      do 190 iifig=1,nafigs
      if(ifig.ne.indx(iifig)) go to 190
      do 180 j=1,nw
  180 jfigs(j,iifig)=idstor(j)
      jsyms(iifig)=jsym
      go to 200
  190 continue
  200 continue
c
      do 250 ifig=1,nafigs
      ndx = indx(ifig)
      call unpck(jfigs(1,ifig),joutfg,2,nbf)
c     write (iw,1400) ifig, ndx, joutfg
c     if (ipunch.gt.0) write (ipunch,'(40i2)') joutfg
      istart = istsef(ndx)
      istop = istsef(ndx+1) - 1
      do 240 isef=istart,istop
      vect2 = vector(iroot,isef)**2
      if(vect2.eq.1.0d0) then
        del = e(iroot)
        go to 220
      else
        del = vect2*(e(iroot)-h(isef))/(1.0d0-vect2)
      endif
      if(isef.eq.istart) go to 220
      write (iw,fmt8) isef, vector(iroot,isef), del
      go to 240
  220 write (iw,fmt7) ndx,joutfg,jsyms(ifig),isef,
     1  vector(iroot,isef),del,dele(ifig)
  240 continue
  250 continue
c
      if(iabort.eq.1) write (iw,2700) maxesc
 2700 format(///5x,43h*****  warning. the configuration list has ,
     1  21hbeen truncated  *****,i10)
c
      if(iciwrt.gt.0) write (iw,2800) (indx(ifig), dele(ifig),
     1  ifig=1,ntotfg)
 2800 format(///5x,'energy contributions for spatial configurations'
     2        //5x, 5(4x,21hnumber  delta(energy))
     3        //(5x, 5(i10,f15.8) ))
c
c  statistical analysis of the energy contributions
c
      write (iw,2900)
 2900 format(///5x,34hstatistical analysis of the energy
     2  //30x,9hnumber of,13x,10hsum of the/10x,13henergy cutoff,
     3  5x,14hconfigurations,5x,20henergy contributions/)
c
      ncut = 0
      sumcut = 0.0d0
      cutoff = 1.0d-8
c
      do 280 ifig=ntotfg,1,-1
      del = dele(ifig)
  260 if(abs(del).ge.cutoff) then
        write (iw,3000) cutoff, ncut, sumcut
        ncut = 0
        sumcut = 0.0d0
        cutoff = 10.0d0*cutoff
        go to 260
      endif
      ncut = ncut + 1
  280 sumcut = sumcut + del
c
  290 write (iw,3000) cutoff, ncut, sumcut
 3000 format(1p,d21.2,i16,0p,f25.8)
c
      return
c
      end
