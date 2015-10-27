*deck @(#)final.f	5.1  11/6/94
      subroutine final(joutfg,msymtp)
c
c  count the number of double-group-adapted functions and determinants
c  which will be generated from the list of spatial configurations
c
c  mesa
c
      common /io/ inp,iout
c
      common /c1/ nbf,maxhsh,maxorb,ni,n,nbfp1,ibdcfg,icfgpn
      dimension joutfg(nbfp1,*)
c
      if( icfgpn .le. 0 ) go to 30
c
c  write configuration list file
c
      do 10 ii = 1,n
   10 write (icfgpn,'(80i1)') (joutfg(ibf,ii), ibf=1,nbfp1)
c
c  count double-group functions and determinants
c
   30 ndet = 0
      ndgf = 0
      do 60  kspc = 1,n
      if( joutfg(1,kspc) .eq. 3 )  go to 60
c
c     (count the singly occupied orbitals)
      ksng = 0
      do 40  ibf = 1,nbf
   40 if( joutfg(ibf,kspc) .eq. 1 )  ksng = ksng + 1
c
      if( ksng .gt. 0 ) go to 50
c
c     (closed-shell case)
      ndet = ndet + 1
      ndgf = ndgf + 1
      go to 60
c
c     (open-shell cases)
   50 if( msymtp .eq. 1 ) then
        ndet = ndet + 2**ksng
        ndfg = ndfg + 2**ksng
      else if( msymtp .eq. 2 ) then
        ndet = ndet + 2**(ksng-1)
        ndfg = ndfg + 2**(ksng-1)
      else
        ndet = ndet + 2**(ksng-1)
        ndgf = ndgf + 4**((ksng-1)/2)
      endif
c
   60 continue
c
      write (iout,90) ndgf,ndet
   90 format(/i16,' double-group-adapted functions will be ',
     1   'generated from this list'//
     2   i16,' determinants will be generated from this list')
      return
      end
