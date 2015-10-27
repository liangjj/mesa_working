*deck @(#)cidata.f	5.1  11/6/94
      subroutine cidata(korboc)
c
c mesa
c 
      common /io/inp,iout
      common /c1/ nbf,maxhsh,maxorb,ni,n,nbfp1,ibdcfg,icfgpn
      common /datarr/ kset(10),lmset(2,10),mmset(100)
      if( korboc .eq. 0 )    return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     read in orbital restrictions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if( korboc .lt. 0  .or.  korboc .gt. 10 )   go to 423
      write (iout,121)
      kkk = 1
      do 200 ii = 1,korboc
      read (5,'(16i5)')   kset(ii),lmset(1,ii),lmset(2,ii)
      lll = kkk + kset(ii) - 1
      if( lll .gt. maxorb )     go to 225
      read (5,'(16i5)')    (mmset(jj),jj=kkk,lll)
      write (iout,217)   kset(ii),lmset(1,ii),lmset(2,ii),
     +   (mmset(jj),jj=kkk,lll)
      kkk = kkk + kset(ii)
  200 continue
      return
  225 ii = ii - 1
      write (iout,241)   ii,korboc
      korboc = ii
      return
  423 write (iout,443) korboc
      stop
  121 format(////11x,'these orbital restrictions will be applied',
     x   ' to all following configurations'/11x,
     x   'until another set of restrictions is encountered')
  217 format(/11x,'no. of orbitals in set =',i5,5x,
     1   'limits of occupation of set =',i5,' to',i5//
     2   11x,'orbitals in set =',30i3/(28x,30i3))
  241 format(/11x,'number of test sets reduced to',i6,' from',i5)
  443 format(///11x,'korboc out of range,  korboc =',i5)
      end
