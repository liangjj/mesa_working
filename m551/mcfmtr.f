*deck @(#)mcfmtr.f	1.1  11/30/90
      subroutine mcfmtr(nbf,nob,ncob,naob,cv,fpq,fij,tpi)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcfmtr.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cc
cc
      dimension cv(nbf,nob),fpq(2)
      dimension fij(nob,naob)
      dimension tpi(2)
c
      common /io/ inp,iout
c
c-----------------------------------------------------------------------
c
c --- description     this routine transforms fock matrix between
c                     basis functions to between active and all
c                     orbitals.
c
c --- input
c
c     nbf             no of basis functions.
c     nob             no of orbitals.
c     ncob            no of core orbitals.
c     naob            no of active orbitals.
c     cv(nbf,nob)     orbitals.
c     fpq(nbf*(nbf+1)/2) f-matrix between basis functions.
c
c --- working storage
c
c     tpi(nbf)
c
c --- output
c
c     fij(nob,naob)   transformed fock matrix.
c
c-----------------------------------------------------------------------
      if (naob .le. 0) return
c     do 1000 i = 1, nob
c1000 write (iout,1001) (cv(j,i),j=1,nbf)
c1001 format(' *mcfmtr cv '//4(1x,f18.6))
c
      lpq = 0
c
      do 100 ip = 1, nbf
         lpq = lpq + ip
         fpq(lpq) = fpq(lpq) * 0.5d0
 100  continue
c
      do 300 i = 1, naob
         ii = ncob + i
         lpq = 1
cc
         do 320 ip = 1, nbf
 320     tpi(ip) = 0.d0
cc
cc
         do 400 ip = 1, nbf
ccc
            do 450 iq = 1, ip
               tpi(ip) = tpi(ip) + fpq(lpq) * cv(iq,ii)
               tpi(iq) = tpi(iq) + fpq(lpq) * cv(ip,ii)
               lpq = lpq + 1
 450        continue
ccc
 400     continue
cc
         do 500 j = 1, nob
            t = 0.d0
ccc
            do 550 ip = 1, nbf
               t = t + tpi(ip) * cv(ip,j)
 550        continue
ccc
            fij(j,i) = t
 500     continue
cc
 300  continue
c
      lpq = 0
c
      do 700 ip = 1, nbf
         lpq = lpq + ip
         fpq(lpq) = fpq(lpq) + fpq(lpq)
 700  continue
c
c     do 800 i = 1, naob
c 800 write (iout,9000) (fij(j,i),j=1,nob)
c9000 format(' *mcfmtr fij '//4(1x,f16.8))
      return
      end
