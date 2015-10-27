      program solver
c note if npvec .ne. nbtot, q must be zero
      parameter (nbig=1200, nsmall=100,nchnl=20)
      complex hdenom(nbig,nbig),smat(nsmall,nsmall),work(nbig)
      complex tmat(nsmall,nsmall),cdotu,cov(nbig,nsmall)
      complex htop(nbig,nsmall),htopp(nbig,nsmall),det,ai
      real kchan(nchnl),hpp(nbig,nbig),hqq(nbig,nbig),
     1    hpq(nbig,nsmall),hfix(nbig),xsecmat(nchnl,nchnl),
     $    eimag(100)
      integer nlm(nchnl),kpvt(nbig)
      equivalence (hqq,hpp),(hdenom,htop)
      data istatic/6hstatic/,ioft/8hoffshell/
      ai=(0.,1.)
      open(5,file='insolve')
      open(6,file='outsolve')
      open(7,file='bmesa',form='unformatted')
c tnr 7/98 read istat from 7
      read(7)istat
      open(77,file='btmat',form='unformatted')
      open(8,file='bsolve',form='unformatted')
      read(5,*)iprint,iopen,nec
      read(5,*)(eimag(i),i=1,nec)
c tnr 7/98 don't read istat from 5
c      read(5,55)istat
 55   format(a8)
      if(istat.eq.istatic)then
      open(99,file='bstatic',form='unformatted')
      endif
      read(7)nchan,nbtot,nfree,(nlm(i),i=1,nchan)
      read(8)npvec,npdim
      write(6,*)' npdim ',npdim
      if(npvec.ne.nbtot)then
      write(6,*)' npvec and nbtot are unequal!', nbtot, npvec
      endif
      write(6,112)nchan,nbtot,nfree,(nlm(i),i=1,nchan)
 112  format(' number of channels=',i3/
     ^ ' dimension of P-space is ',i5/
     ^ ' number of free (lm) channels is =',i4/
     ^ ' number of lm terms per channel:'/(20i3))
      ntot=nbtot+nfree
      if(nchan.gt.nchnl.or.nfree.gt.nsmall.or.
     ^ ntot.gt.nbig)then
      write(6,*)' check dimensions'
      stop
      endif
      nall=ntot+npdim
      iw=nbtot+1
      read(7)nener
      do 1 ie=1,nener
      do 100 i=1,nbig
      do 101 j=1,nbig
 101     hdenom(i,j)=0.
      do 102 j=1,nsmall
 102  htopp(i,j)=0.
 100     continue
      read(7)(kchan(ic),ic=1,nchan)
      write(6,111) (kchan(i),i=1,nchan)
111   format(///,' ********************************',
     # ' new energy ********************************',
     # //,' channel momenta:',6e12.5)
      if(nbtot.eq.npvec)then
      read(7)((hdenom(i,j),i=1,ntot),j=iw,ntot)
      do 114 i=1,ntot
      do 114 j=iw,ntot
 114  hdenom(j,i)=hdenom(i,j)   
      else
      read(7)((hdenom(i,j),i=1,ntot),j=1,ntot)
      endif
      call rdbinsqr(hpp,nbig,npvec,8)
      if(istat.eq.istatic)then
      call rdbinsqr(hpp,nbig,nbtot,99)
      endif
      if(nbtot.eq.npvec.or.istat.eq.istatic)then
      do 2 i=1,nbtot
      do 2 j=1,nbtot
 2       hdenom(i,j)=hpp(i,j)
      endif
      if(npdim.gt.0)then
      call rdbinsqr(hqq,nbig,npdim,8)
      do 3 i=1,npdim
      do 3 j=1,npdim
      ii=i+ntot
      jj=j+ntot
 3    hdenom(ii,jj)=hqq(i,j)

      do 4 i=1,npvec
      read(8)(hpq(j,i),j=1,npdim)
      do 5 j=1,npdim
      hdenom(i,j+ntot)=hpq(j,i)
 5    hdenom(j+ntot,i)=hpq(j,i)
 4    continue
      endif

      read(7)((htopp(i,j),i=1,ntot),j=1,nfree)
c tnr 7/98
      if(istat.eq.ioft)then
      read(7)nread
      write(6,*)' no. occupied plus scattering terms', nread
      read(7)((cov(i,j),i=1,nread),j=1,nfree)
      do 6 k=1,nbtot
      do 6 i=1,npdim
      do 6 j=1,nfree
 6       htopp(i+ntot,j)=htopp(i+ntot,j)+hpq(i,k)*cov(k,j)
      endif
      if(iprint.ne.0) then
      write(6,109)
109   format(//' numerator matrix of (h-e)')
      do 61 i=1,nfree
61    write(6,108) i,(htopp(j,i),j=1,20)
      endif
c
c save numerator and demoninator matrices
c
      do 50 i=1,nbig
         do 51 j=1,nsmall
 51         cov(i,j)=htopp(i,j)
         do 50 j=1,nbig
 50         hqq(i,j)=hdenom(i,j)
c
c inner loop on complex energies
c
      do 1000 iec=1,nec
         write(6,*)' Imaginary part of Z=',eimag(iec)
      do 54 i=1,nbig
         do 54 j=1,nbig
 54         hdenom(i,j)=hqq(i,j)
c add imaginary term to energy in hqq-E
      do 33 i=1,npdim
 33    hdenom(ntot+i,ntot+i)=hqq(ntot+i,ntot+i)-ai*eimag(iec)
      if(iprint.ne.0) then
      write(6,107)
107   format(//' denominator matrix of (h-e)')
      do 60 i=1,20
60    write(6,108) i,(hdenom(j,i),j=1,20)
108   format(1x,i3,6e12.5,/,(4x,6e12.5))
      write(6,*)' last part'
      do 660 i=1,20
660    write(6,108) i,(hdenom(j+ntot,i+ntot),j=1,20)
      endif
      call cgefs(hdenom,nbig,nall,htopp(1,1),1,ind,work,kpvt)
      do 666 i=2,nfree
 666  call cgefs(hdenom,nbig,nall,htopp(1,i),2,ind,work,kpvt)
      if(iprint.ne.0) then
      write(6,116)
 116  format(//' solution matrix )')
      do 662 i=1,nfree
 662    write(6,108) i,(htopp(j,i),j=1,20)
      endif
      do 103 i=1,nbig
      do 103 j=1,nsmall
 103  htop(i,j)=0.
      if(iec.eq.1)then
         read(7)((htop(i,j),i=1,ntot),j=1,nfree)
         do 52 i=1,nbig
            do 52 j=1,nsmall
 52            hpq(i,j)=htop(i,j)
         if(iprint.ne.0) then
            write(6,110)
 110        format(//' conjugate numerator matrix of (h-e)')
            do 62 i=1,nfree
 62            write(6,108) i,(htop(j,i),j=1,20)
         endif
      else
         do 53 i=1,nbig
            do 53 j=1,nsmall
 53            htop(i,j)=hpq(i,j)
      endif
      do 500 ilm=1,nfree
      do 500 jlm=1,nfree
      tmat(ilm,jlm)=cdotu(nall,htopp(1,jlm),1,htop(1,ilm),1)
500   continue
      if(iec.eq.1)then
         read(7)((smat(i,j),i=1,nfree),j=1,nfree)
         if(iprint.ne.0) then
            write(6,113)
 113        format(//' first Born matrix')
            do 63 i=1,nfree
 63            write(6,108) i,(smat(j,i),j=1,nfree)
         endif
      endif
      do 501 ilm=1,nfree
      do 501 jlm=1,nfree
      tmat(ilm,jlm)=-2.*(-tmat(ilm,jlm)+smat(ilm,jlm))
 501  continue
      write(6,115)
 115  format(//' T-Matrix')
      do 64 i=1,nfree
64    write(6,108) i,(tmat(j,i),j=1,nfree)
      istart = 0
      do 879 ic=1,nchan
      jstart = 0
      ni = nlm(ic)
      do 878 jc=1,nchan
      nj=nlm(jc)
      ist=istart+1
      ifin=istart+ni
      jst=jstart+1
      jfin=jstart+nj
      write(77)ic,jc,ni,nj,kchan(ic),kchan(jc)
c      write(88,177)ic,jc,ni,nj,kchan(ic),kchan(jc)
 177  format(4i5,2f20.10)
      write(77)((tmat(ii,jj),ii=ist,ifin),jj=jst,jfin)
c      write(88,277)((tmat(ii,jj),ii=ist,ifin),jj=jst,jfin)
 277  format(4e20.10)
      summod = 0.
      do 877 ilm=1,ni
      do 877 jlm=1,nj
      isub = istart + ilm
      jsub = jstart + jlm
877   summod = summod + cabs(tmat(isub,jsub))**2
      xsecmat(ic,jc) = 4.0*3.141592654*summod/kchan(ic)**2
      jstart = jstart + nlm(jc)
878   continue
      istart = istart + nlm(ic)
879   continue
       write(6,874)
874   format(//,' total cross sections: row index = initial chnl,',
     # ' column index = final chnl')
      write(6,871) (i,i=1,nchan)
871   format(4x,6(6x,i2,4x))
      do 873 i=1,nchan
      write(6,872) i, (xsecmat(i,j),j=1,nchan)
873   continue
872   format(1x,i2,1x,6e12.5/(4x,6e12.5))
 1000 continue
 1    continue
      call exit
      end
      subroutine rdbinsqr(amat,nrow,ncol,iunit)
      dimension amat(nrow,ncol)
      do 1 i=1,ncol
  1   read(iunit) (amat(j,i),j=1,ncol)
      return
      end
c*    subroutine matinv ( a, n, b, m, det, nmax )
c**   <matinv> -- argonne nat. lab. program anl f453s -- complex matrix
c**             inversion with accompanying solution of linear equations
c     inverts matrix a(n,n) (with dimension statement a(nmax,n)
c     inverse of a is returned in a
c     solves linear equations of the form ax=b with b dimensioned b(nmax
c     x is returned in b
c
       subroutine matinv ( a, n, b, m, det, nmax )
      complex      a, b, det, swap, t
      real  amax, temp
      real  cabs
      real  pivot
      dimension a(nmax,n), b(nmax,m)
      common / f402 / pivot(4096), index(4096)
c     -----------------------------------------------------------------
c     initialize determinant and pivot element array
c     -----------------------------------------------------------------
      det = 1e0
      do 20 i = 1, n
      pivot(i) = 0.
   20 continue
c     perform successive pivot operations ( grand loop )
c     -----------------------------------------------------------------
      do 550 i = 1, n
c     -----------------------------------------------------------------
c     search for pivot element and extend determinant partial product
c     -----------------------------------------------------------------
      amax = 0e0
      do 105 j = 1, n
      if ( pivot(j) .ne. 0. ) go to 105
      do 100 k = 1, n
      if ( pivot(k) .ne. 0. ) go to 100
      temp = cabs ( a(j,k) )
      if ( temp .lt. amax ) go to 100
      ir = j
      ic = k
      amax = temp
  100 continue
  105 continue
      index(i) = 4096*ir + ic
      j = ir
      t = a(j,ic)
      det = det*t
c     -----------------------------------------------------------------
c     return if matrix is singular (zero pivot) after column interchange
c     -----------------------------------------------------------------
      if ( cabs(det) .eq. 0e0 ) go to 600
c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
      pivot(ic) = amax
c     -----------------------------------------------------------------
c     interchange rows to put pivot element on diagonal
c     -----------------------------------------------------------------
      if ( ir .eq. ic ) go to 260
      det = -det
      do 200 k = 1, n
      swap = a(j,k)
      a(j,k) = a(ic,k)
      a(ic,k) = swap
  200 continue
      if ( m .le. 0  ) go to 260
      do 250 k = 1, m
      swap = b(j,k)
      b(j,k) = b(ic,k)
      b(ic,k) = swap
  250 continue
c     -----------------------------------------------------------------
c     divide pivot row by pivot element
c     -----------------------------------------------------------------
  260 do 350 k = 1, n
      if ( k .eq. ic ) a(ic,k) = 1.
      a(ic,k) = a(ic,k)/t
  350 continue
      if ( m .le. 0 ) go to 380
      do 370 k = 1, m
      b(ic,k) = b(ic,k)/t
  370 continue
c     -----------------------------------------------------------------
c     reduce non-pivot rows
c     -----------------------------------------------------------------
  380  do 550 j = 1, n
       if ( j .eq. ic ) go to 550
       t = a(j,ic)
       a(j,ic) = 0e0
       do 450 k = 1, n
       a(j,k) = a(j,k) - a(ic,k)*t
  450 continue
      if ( m .le. 0 ) go to 550
      do 500 k = 1, m
      b(j,k) = b(j,k) - b(ic,k)*t
  500 continue
  550 continue
      i = n
c     -----------------------------------------------------------------
c     interchange colums after all pivot operations have been performed
c     -----------------------------------------------------------------
  600 do 710 il = 1, i
      ii = i + 1 - il
      k = index(ii)/4096
      ic = index(ii) - 4096*k
      if ( k .eq. ic ) go to 710
      do 705 j = 1, n
      swap = a(j,k)
      a(j,k) = a(j,ic)
      a(j,ic) = swap
  705 continue
  710 continue
      return
      end
