      subroutine bforthog(hpvbt,lmtop,nbfmax,nchnl,nstate,hpvbtx,
     1hb, nchan,nlm,nbas,nmotot,iprint,nbch)
c
c construct free-bound matrix of (h-e) for free functions
c orthogonalized to bound basis
c
      implicit real*8 (a-h,o-z)
      real*8 hb(nbfmax,nbfmax)
      complex*16 hpvbtx(lmtop,nbfmax,nchnl)
      complex*16 hpvbt(lmtop,nbfmax,nchnl)
      integer nlm(nchnl),nscat(nchnl),nbas(nchnl)
      integer nbch(nbfmax,nchnl)
c
c do the orthogonalization
c
      do 400 ic=1,nchan
      nlmic=nlm(ic)
      nbkc=nbas(ic)
      do 402 kbc=1,nbkc
         jw=nbch(kbc,ic)
      do 402 ilm=1,nlmic
      do 402 jsc=1,nmotot
402   hpvbtx(ilm,jsc,ic) = hpvbtx(ilm,jsc,ic)
     1 -hpvbt(ilm,jw,ic)*hb(jw,jsc)
400   continue
c
      if(iprint.ne.0) then
      do 200 ic=1,nchan
      nlmic=nlm(ic)
      write(6,107)ic
107   format(//,' bound-(orth)free matrix for channel:',i4)
      nsjc=nmotot
      do 200 ilm=1,nlmic
200   write(6,101) ilm,(hpvbtx(ilm,j,ic),j=1,nsjc)
101   format(1x,i3,6e12.5,/,(4x,6e12.5))
      endif
      return
      end
