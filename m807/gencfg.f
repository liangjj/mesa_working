*deck @(#)gencfg.f	5.1  11/6/94
      subroutine gencfg(inbcfg,msymbf,iorbs,norboc,nbcoc,joutfg,
     1  istore,nhash,nelar,norbar,indx,indxu,korboc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  read the data for the various configurations and control their
c  generation and storage; generate configurations by calling
c  config, and then check to see that they obey all restrictions by
c  calling chkcfg
c     n refers to configuration *n* in joutfg
c     k refers to orbital group *k*
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*mdc*if harris
*      integer*6 nhash
*mdc*endif
c
c  mesa
c
      common /io/inp,iout
c
      common /c1/ nbf,maxhsh,maxorb,ni,n,nbfp1,ibdcfg,icfgpn
      common /c2/ mxop,mnop,ntexct,maxcfg,kpar,mxopn
      common /c3/ nel,nbcfg,ncfgde,ncfgin,nbfp3,mxnset
      common /c4/ norb,neli,i,nelu,nexcit,nexctu,ifirst
      dimension inbcfg(nbfp3),msymbf(*),iorbs(*),norboc(*),nbcoc(*),
     1  joutfg(nbfp1,*),istore(*),nhash(*),nelar(*),norbar(*),
     2  indx(*),indxu(*)
      mxindx = maxcfg*nbfp1
      ni = 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     read in symmetry of basis functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      read (5,'(16i5)')  (msymbf(i),i=1,nbf)
      write (iout,117)   (ibf,msymbf(ibf),ibf=1,nbf)
  117 format('1',10x,'basis function number',10x,'symmetry type'/
     x     /(i21,i26))
      call       cidata(korboc)
      if( ncfgin .eq. 0 )   go to 153
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     the following particular configurations will start the list of
c     configurations
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write (iout,133)
  133 format(/11x,'specific configurations put on the list'//)
      do 154  j=1,nbf
  154 inbcfg(j) = 0
      ntexct = 1
      do 134  ii=1,ncfgin
      read (5,'(80i1)')       (joutfg(j,n),j=1,nbf)
      neln = 0
      do 152  ibf=1,nbf
  152 neln = neln + joutfg(ibf,n)
      if( neln .ne. nel )   go to 134
      call       chkcfg(inbcfg,msymbf,joutfg,nhash,0)
      if( ibdcfg .eq. 0 )   n = n + 1
  134 continue
      call       wrtcfg(joutfg)
  153 if( ncfgde . eq. 0 )   go to 124
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     the following configurations will be stored in the nhash table
c     but will not be listed in joutfg.  this prevents these
c     configurations from appearing in the final list.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write (iout,128)
  128 format(/11x,'the following configurations will not be al',
     x  'lowed to appear in the final list of configurations'//)
      do 129 l=1,ncfgde
      read (5,'(80i1)')       (inbcfg(j),j=1,nbf)
      if(newcfg(inbcfg,nhash).eq.0) write (iout,76) (inbcfg(j),j=1,nbf)
  129 continue
      write (iout,'(5(/))')
  124 if( nbcfg .eq. 0 )   return
c
      do 27  l=1,nbcfg
      lk = mxindx
c
c     read in the data for the next configuration
c
      read  (5,'(80i1)') inbcfg
      neln = 0
      do 119 ibf=1,nbf
  119 neln = neln + inbcfg(ibf)
      if( neln .ne. nel )   go to 27
      ntexct = inbcfg(nbfp1)
      nsets  = inbcfg(nbfp1+1)
      write (iout,'(40i2)') nsets
      write (iout,'(1x,80i1)') inbcfg
      if( nsets .gt. mxnset )   go to 39
      if( nsets ) 38,120,36
  120 do 37  j=1,nbf
   37 iorbs(j) = j
      iorbs(nbfp1) = 0
      go to 38
   36 read (5,'(26i3)')       (iorbs(j),j=1,nbfp1)
   38 write (iout,5)
      write (iout,76) (inbcfg(j),j=1,nbfp1)
      if( inbcfg(nbfp3) .le. 0 )   go to 46
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     read in a new set of orbital restrictions.  these will be
c     applied to all subsequent configurations generated until
c     another set of orbital restrictions is encountered.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      korboc = inbcfg(nbfp3)
      call       cidata(korboc)
      write (iout,'(///)')
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     count number of orbitals and electrons for each set
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   46 k = 1
      kk = 1
      nelar(k) = 0
      norbar(k) = 0
      do 6 j=1,nbfp1
      if( iorbs(j) )  160,150,9
  150 iorbs(j) = inbcfg(nbfp1)
      go to 170
  160 iorbs(j) = -iorbs(j)
  170 iorbs(j) = min( iorbs(j), ntexct )
      write (iout,2000 ) k,(iorbs(jj),jj=kk,j)
      iorbs(j) = -iorbs(j)
      kk = j + 1
      if( k .ge. nsets )   go to 10
      k = k + 1
      nelar(k) = 0
      norbar(k) = 0
      go to 6
    9 nelar(k) = nelar(k) + inbcfg( iorbs(j) )
      norbar(k) = norbar(k) + 1
    6 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     index arguments for the first call to config for this orbital
c     group.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   10 kk = 1
      k = 1
   19 i = 0
      ifirst = 0
      norb = norbar(k)
      norbkk = norb + kk
      neli = nelar(k)
      nelu = 0
      nexcit = -iorbs(norbkk)
      do 12  j=1,norb
      nbcoc(j) = inbcfg( iorbs(kk) )
   12 kk = kk + 1
      kk = kk + 1
      if( nsets .gt. k )   go to 11
      if( nsets .gt. 1 )   go to 48
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     generate and check configurations for basic configurations
c     which have only one orbital group.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mini=nbf
   15 do 8  j=1,nbf
    8 joutfg(j,n) = inbcfg(j)
   14 call       config(norboc,nbcoc)
      mini = min(mini,i)
      if( i .le. 0 )   go to 30
      do 13  j=1,norb
   13 joutfg( iorbs(j) ,n) = norboc(j)
      call       chkcfg(inbcfg,msymbf,joutfg,nhash,korboc)
      if( ibdcfg .eq. 1 )   go to 14
      n = n + 1
      if( n .gt. maxcfg)   go to 31
      go to 15
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     generate all configurations within the orbital group and store
c     in the end of joutfg
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   11 jj = 0
    7 call       config(norboc,nbcoc)
      if( i .le. 0 )   go to 16
      jj = jj + 1
      do 17  j=1,norb
      istore(lk) = norboc(j)
   17 lk = lk - 1
      go to 7
   16 write (iout,42) k,jj
      k = k + 1
      indx(k) = lk
      go to 19
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     generate the configurations for the last orbital group.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   48 kj = kk - norb - 1
      indx(1) = mxindx
      lk = lk/nbf
      mini=nbf
      do 49  j=1,nbf
   49 joutfg(j,n) = inbcfg(j)
   50 call       config(norboc,nbcoc)
      mini = min(mini,i)
      if( i .le. 0 )   go to 30
      do 51  j=1,norb
      joutfg(iorbs(kj),n) = norboc(j)
   51 kj = kj + 1
      kj = kj - norb
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     form all possible configurations from this specific
c     configuration for the last group, and the different member
c     configurations from the other orbital groups.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   47 k = 1
      jj = 1
   18 kk = indx(k)
   26 numorb = norbar(k)
   43 do 20  j=1,numorb
      joutfg( iorbs(jj) ,n) = istore(kk)
      jj = jj + 1
   20 kk = kk - 1
      jj = jj + 1
      indxu(k) = kk
      k = k + 1
      if( k .ge. nsets )   go to 21
      go to 18
   21 k = k - 1
      call       chkcfg(inbcfg,msymbf,joutfg,nhash,korboc)
      if( ibdcfg .eq. 1 )   go to 22
      n = n + 1
      if( n .gt. lk )   go to 31
      do 23  j=1,nbf
   23 joutfg(j,n) = joutfg(j,n-1)
   22 jj = jj - numorb - 1
   24 if( indxu(k) .gt. indx(k+1) )   go to 43
   25 k = k - 1
      if( k .lt. 1 )   go to 50
      jj = jj - norbar(k) - 1
      if( indxu(k) .le. indx(k+1) )   go to 25
      kk = indxu(k)
      go to 26
   30 call       wrtcfg(joutfg)
      write (iout,'(5(/))')
   27 continue
c
      return
   31 write (iout,32) n,mini
      stop
   39 write (iout,40) nsets,mxnset
      stop
    5 format(/11x,'basic configuration and excitation number'//)
   32 format(/11x,'joutfg is full.',i6,' configurations have been',
     x   ' generated.'//11x,'excitations are being made from ',
     x   'orbital',i4,' of the last group.')
   40 format(/11x,'nsets =',i2,10x,'max number of sets =',i2
     x      //11x,'execution terminated')
   42 format(/11x,'orbital group',i2,' :',i10,' configurations')
   76 format(11x,50i2/)
 2000 format(/11x,'orbital group',i2/,(20i3/))
      end
