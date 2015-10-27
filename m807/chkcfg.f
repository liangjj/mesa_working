*deck @(#)chkcfg.f	5.1  11/6/94
      subroutine chkcfg(inbcfg,msymbf,joutfg,nhash,korboc)
*mdc*if harris
*      integer*6 nhash
*mdc*endif
      common /c1/ nbf,maxhsh,maxorb,ni,n,nbfp1,ibdcfg,icfgpn
      common /c2/ mxop,mnop,ntexct,maxcfg,kpar,mxopn
      dimension inbcfg(*),msymbf(*),joutfg(nbfp1,*),nhash(*),mtab(8,8)
      data mtab /1,2,3,4,5,6,7,8,
     1           2,1,4,3,6,5,8,7,
     2           3,4,1,2,7,8,5,6,
     3           4,3,2,1,8,7,6,5,
     4           5,6,7,8,1,2,3,4,
     5           6,5,8,7,2,1,4,3,
     6           7,8,5,6,3,4,1,2,
     7           8,7,6,5,4,3,2,1/
      ibdcfg = 0
      nopenu = 0
      nexctu = 0
      jsym = 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     check for number of opens, total number of excitations
c     allowed, and parity (if d2h).  save symmetry type.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 10 ibf = 1,nbf
      m = inbcfg(ibf) - joutfg(ibf,n)
      if( m .gt. 0 )  nexctu = nexctu + m
      if( joutfg(ibf,n) .ne. 1 )  go to 10
      nopenu = nopenu + 1
      jsym = mtab(jsym,msymbf(ibf))
   10 continue
      if( nexctu .le. ntexct  .and.  nopenu .le. mxop  .and.
     1    nopenu .ge. mnop    .and.  jsym/5 .eq. kpar )  go to 20
      ibdcfg = 1
      return
   20 joutfg(nbfp1,n) = jsym
      if( korboc .eq. 0 )  go to 30
      call       chorre(joutfg,korboc)
      if( ibdcfg .eq. 1 )   return
   30 ibdcfg = newcfg(joutfg(1,n),nhash)
      if( ibdcfg .eq. 0 )   mxopn = max(mxopn,nopenu)
      return
      end
