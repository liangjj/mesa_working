*deck @(#)pm807.f	5.1  11/6/94
c
c*******************************************************************
c
c   this computer program contains work performed partially at
c   lawrence livermore national laboratory under
c   the auspices of the office of basic energy sciences,
c   division of chemical sciences, u.s. department of energy,
c   under contract w-7405-eng-48.
c
c   these programs may not be (re)distributed without the
c   written consent of the columbus program group.
c
c   since these programs are under development, correct results
c   are not guaranteed.
c
c*******************************************************************
c
       subroutine pm807
c
c pitzer's configuration generator program modified
c by m. braunstein to fit into mesa, nov. 91
c incude cases of c1, cs, c2 for even no. of electrons. june 92
c
c  configuration generator for ci calculations
c
c  version log:
c  04-nov-89 mdc revisions (rmp)
c  18-aug-89 mdc version (rmp)
c  1984-1985 written, starting with mqm23 (rmp,nww)
c
c  cmdc info:
c  keyword   description
c  -------   -----------
c  vax       vax code.
c  cray      cray code. includes "d" exponents in constants.
c  ibm       ibm code.
c  sun       sun code.
c  stellar   stellar code.
c  ibmmvs    mvs specific code.
c  crayctss  ctss specific code.
c
      real*8 time1, time2
*mdc*if harris
*      integer*6 nhash
*mdc*endif
      character*24 dattim
      character*80 title
      common /c1/ nbf,maxhsh,maxorb,ni,n,nbfp1,ibdcfg,icfgpn
      common /c2/ mxop,mnop,ntexct,maxcfg,kpar,mxopn
      common /c3/ nel,nbcfg,ncfgde,ncfgin,nbfp3,mxnset
c
c mesa
c
      common /io/ inp,iout
c
c  iasz=(nbf+1)*maxcfg+5*nbf+4
      parameter (iasz=272402)
      dimension ia(iasz)
      dimension nhash(24576),nelar(6),norbar(6),indx(6),indxu(6)
c
*mdc*if crayctss
*      call link('unit99=tty//')
*mdc*endif
c
c  open input file
c
*mdc*if harris
*c  use preconnected unit 5
*mdc*else
      open (5,file='genin',status='old')
*mdc*endif
c
c  open the output listing file
c
*mdc*if harris ibmmvs
c  use preconnected unit 6
*mdc*else
c mesa
c
c      open (6,file='genout')
*mdc*endif
      write (iout,90)
   90 format('1',28x,'program "cgdbg" 4.0b2'/
     1  29x,'columbus program system'/
     2  29x,'configuration generation program'//
     3  29x,'version date: 04-nov-89'//
     &  ' maintained by:'/
     &  t5,'russell m. pitzer'/
     &  t5,'department of chemistry'/
     &  t5,'the ohio state university'/
     &  t5,'columbus, oh  43210'/
     &  t5,'bitnet: ts0775 at ohstmvsa')
c
      maxopn = 8
      maxnel = 72
      maxorb = 100
      mxnset = 6
      maxhsh = 24576
c
  200 call btime(dattim,time1)
      n = 1
      read  (5,'(a80)') title
      write (iout,'(//10x,a80,5x,a24//)') title, dattim
      read  (5,'(16i5)') nbf,nel,ksym,nbcfg,ncfgin,ncfgde,icfgpn,
     x                   mxop,msymtp,korboc
      if( nel .lt. 1  .or.  nel .gt. maxnel )  then
        write (iout,400)  nel,maxnel
  400   format(////11x,'nel is in error'
     x           //11x,'nel,maxnel =',2i6)
        stop
      endif
      if( mxop .gt. maxopn .or.  mxop .eq. 0 ) mxop = maxopn
      mxop = mxop - mod(abs(nel-mxop),2)
c
c  set parity (kpar=0 for none or g, =1 for u) and mnop
c
      if( ksym .le. msymtp .and. mod(nel,2) .eq. 1 ) go to 410
      if( msymtp .ne. 1 )  go to 130
      if( msymtp .ne. 1 )  go to 130
c
c  c1
      kpar = 0
      if( ksym - 1 )  410,110,410
  110 mnop = 0
      go to 330
c
  130 if( msymtp .ne. 2 )  go to 210
c
c  cs or c2
      kpar = 0
      if( ksym - 1 )  410,140,150
  140 mnop = 0
      go to 330
  150 if( ksym - 3 )  160,410,410
  160 mnop = 2
      go to 330
c
  210 if( msymtp .ne. 4 )  go to 260
c
c  c2v or d2
      kpar = 0
      if( ksym - 1 )  410,220,230
  220 mnop = 0
      go to 330
  230 if( ksym - 5 )  240,250,410
  240 mnop = 2
      go to 330
  250 mnop = 1
      go to 330
c
  260 if( msymtp .ne. 8 )  go to 430
c
c  d2h
      if( ksym - 1 )  410,270,280
  270 mnop = 0
      go to 310
  280 if( ksym .gt. 4 )  go to 290
      mnop = 2
      go to 310
  290 if( ksym .gt. 8 )  go to 300
      mnop = 2
      go to 320
  300 mnop = 1
      if( ksym - 10 )  310,320,410
  310 kpar = 0
      go to 330
  320 kpar = 1
c
  330 mxopn = mnop
      nbfp1 = nbf + 1
      nbfp3 = nbf + 3
      do 340  icl = 1,maxhsh
  340 nhash(icl) = 0
c
c  allocate space to nbf-dependent arrays
c
      jnbcfg = 1
      nsymbf = jnbcfg+nbfp3
      korbs  = nsymbf+nbf
      iorboc = korbs +nbfp1
      ibcoc  = iorboc+nbf
      koutfg = ibcoc +nbf
      jstore = koutfg
c
      maxcfg = (iasz+1-koutfg)/nbfp1
      if( maxcfg .le. 0 )  go to 450
      if( nbcfg .gt. maxcfg  .or.  nbcfg .lt. 0 )  go to 470
      write (iout,350) nbf,nbcfg,nel,ksym,ncfgin,ncfgde,mnop,mxop,
     x              icfgpn,korboc,msymtp,maxcfg
  350 format(/
     x  11x,'nbf.................',i5,5x,'nbasic configs......',i5//
     x  11x,'nel.................',i5,5x,'sym. of state.......',i5//
     x  11x,'nconfigs inserted...',i5,5x,'nconfigs deleted....',i5//
     x  11x,'min open shells.....',i5,5x,'max open shells.....',i5//
     x  11x,'icfgpn..............',i5,5x,'korboc..............',i5//
     x  11x,'no. sym. types......',i5,5x,'maxcfg............',i7//)
      call gencfg(ia(jnbcfg),ia(nsymbf),ia(korbs),ia(iorboc),
     1  ia(ibcoc),ia(koutfg),ia(jstore),nhash,nelar,norbar,indx,
     2  indxu,korboc)
      n = n - 1
c
      if( icfgpn .le. 0 ) go to 370
      open (icfgpn,file='configs')
      write (icfgpn,'(a80)') title
      write (icfgpn,'(8i5)') nbf,nel,ksym,n,mxopn,msymtp
c
  370 call       final(ia(koutfg),msymtp)
      call btime(dattim,time2)
      ijsec=(time2-time1)
      write (iout,380) ijsec, dattim
  380 format(/i16,' seconds cpu time',5x,a24)
      read  (5,'(i5)') itest
      if(itest.ne.0) go to 200
c
c mesa
c      stop 'end of cgdbg'
c
c exit gracefully
c
      call chainx(0)
c
      return
c
c the following is error handling, in future use lnkerr
c
  410 write (iout,420)  ksym
  420 format(////11x,'ksym is in error'
     x         //11x,'ksym =',i6)
      stop
  430 write (iout,440)  msymtp
  440 format(////11x,'msymtp is in error.  only 1, 2, 4, or 8'
     x               ' may be used.'//11x,'msymtp =',i6)
      stop
  450 write (iout,460)  maxcfg
  460 format(////11x,'more space needed.  increase iasz.'
     x         //11x,'maxcfg =',i6)
      stop
  470 write (iout,480)  nbcfg,maxcfg
  480 format(////11x,'nbcfg is in error, or more space needed ',
     x  '(increase iasz)'//11x,'nbcfg,maxcfg =',2i6)
      stop
      end
