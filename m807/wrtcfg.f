*deck @(#)wrtcfg.f	5.1  11/6/94
      subroutine wrtcfg(joutfg)
c
c mesa
c
      common /io/inp,iout
c
      common /c1/ nbf,maxhsh,maxorb,ni,n,nbfp1,ibdcfg,icfgpn
      dimension joutfg(nbfp1,*)
      if( ni .ne. n )   then
        write (iout,'(/)')
        do 4 ii = ni,n-1
    4   write (iout,'(i10,2x,80i1)') ii,(joutfg(ibf,ii),ibf=1,nbfp1)
        ni = n
      endif
      return
      end
