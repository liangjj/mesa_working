*deck %W%  %G%
      subroutine mccgb(nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,lok,mix,g1,g2,
     $     nda1,lda1,nf35,nf36,
     $     cg,cgx,ipqrs,r,fpq,tpi,bufix,lbufso,
     $     rabcx,rabix,raibx,incor,sg,tg)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
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
      real*8 bufix(lbufso),sg(*),tg(*)
      real*8 rabcx(*),rabix(*),raibx(*)
      dimension ntab(21),mtab(21),ltab(21),ktab(21)
      dimension nbf(2),nob(2),ncob(2),naob(2),cv(2)
      dimension ipqrs(2),ihd(24),nfob(2)
      dimension locsym(2),len(2),lok(2),mix(2),g1(2),g2(2)
      dimension cg(2)
      dimension cgx(2),r(2),tpi(2),fpq(2)
c
      common /io/ inp,iout
c
c
c-----------------------------------------------------------------------
c
c --- description     this program oversees the construction of
c                     c*g.
c
c --- input
c
c     nsym            no of symmetries.
c     nbf(nsym)       number of basis functions.
c     nob(nsym)       number of orbitals.
c     ncob(nsym)      number of core orbitals.
c     naob(nsym)      number of active orbitals.
c     cv(--)          orbitals.
c     locsym(nsym)    pointers into the arrays len and loc for the
c                     update vectors.
c     len(--)         no of nonzero componenets of update vectors.
c     lok(--)         location of update vectors in arrays mix and cg.
c     mix(--)         indices of nonzero components of update vectors.
c     g1(--)          one-electron density matrix elements.
c     g2(--)          two-electron density matrix elements.
c     nda1            direct access dataset containing the core
c                     fock-matrix.
c     lda1            address of the f-matrix in nda1.
c     nf35            fortran number of dataset containing the
c                     two electron integrals needed in the construction
c                     of the updated integrals.
c     nf36            fortran number of dataset containing the
c                     integrals needed to make core contributions
c                     to the updated integrals.
c
c --- output
c
c     cg(--)          c*g, packed according to len and mix
c
c --- working storage
c
c     r(x)            x = max(naob(isym)) * max(nob(isym)).
c     cgx(x)          x = number of occupied orbitals times number of
c                         orbitals, summed over all symmetries.
c     ipqrs(npqrs)
c     fpq(x)          x = sum(nbf(isym)*nbf(isym)+nbf(isym))/2
c     tpi(x)          max(nbf(isym)).
c     ihead(14,npqrs)
c
c-----------------------------------------------------------------------
c
c-------------------------------c
c     set up pointer arrays     c
c-------------------------------c
      nobt = 0
      ncobt = 0
      nvc = 0
      npq = 0
      nim = 0
      nnn = 0
      ntab(1) = 0
      mtab(1) = 0
      ltab(1) = 0
      ktab(1) = 0
c
c       write(iout,*)' mccga  nbf nob nfo nao nco '
c        write(iout,*) nbf(1),nob(1),nfob(1),naob(1),ncob(1)
c
      do 100 l = 1, nsym
         nobt = nobt + nob(l)
         ncobt = ncobt + ncob(l)
         nvc = nvc + nbf(l) * nfob(l)
         ntab(l) = nvc
         nvc = nvc + nbf(l) * nob(l)
         npq = npq + nbf(l) * (nbf(l) + 1) / 2
         mtab(l+1) = npq
         nnn = nnn + naob(l) * (naob(l) + 1) / 2
         ktab(l+1) = nnn
         nooc = ncob(l) + naob(l)
         nim = nim + nob(l) * nooc
         ltab(l+1) = nim
 100  continue
      ntab(nsym+1)=nvc
c-------------------------------------------c
c      zero out slots for expanded c*g      c
c-------------------------------------------c
      nim = ltab(nsym+1)
      do 200 i = 1, nim
 200     cgx(i) = 0.d0
c-----------------------------------------------------c
c                                                     c
c     calculate one electron contributions to c*g     c
c                                                     c
c-----------------------------------------------------c
      npq = mtab(nsym+1)
c
c
      call iosys('read real mc_core_fock from mcscr',npq,fpq,0,' ')
c
c
c      write (iout,1900) (fpq(i),i=1,npq)
 1900 format(' *mccgb fpq '//4(1x,f16.8))
c
      do 1000 l = 1, nsym
         if (naob(l) .le. 0) go to 1000
         call mcfmtr(nbf(l),nob(l),ncob(l),naob(l),cv(ntab(l)+1),
     $               fpq(mtab(l)+1),r,tpi)
         lla = ltab(l) + ncob(l) * nob(l) + 1
         call mccg1e(naob(l),nob(l),cgx(lla),r,g1(ktab(l)+1))
 1000 continue
c
c     write (iout,1357) (cgx(ii),ii=1,nim)
c1357 format(' *mccgb cgx '//4(1x,f16.8))
c---------------------------------------c
c     add core contribtuions to c*g     c
c---------------------------------------c
      if (ncobt .le. 0) go to 2200
c
      call iosys('read integer "header abix" from rwf',14,ihd,0,' ')
      call iosys('rewind abix on rwf',0,0,0,' ')
      ml=ihd(1)
      m1=ihd(2)
      m2=ihd(3)
      m3=ihd(4)
      m4=ihd(5)
      n3=ihd(6)
      nm3=ihd(7)
      n4=ihd(8)
      nm4=ihd(9)
      n1=ihd(10)
      n2=ihd(11)
      n12=ihd(12)
      n34=ihd(13)
      ipqr=ihd(14)
c
c2000 read (nf36,end=2200) ml, m1, m2, m3, m4, n3, nm3, n4, nm4,
c    1      n1, n2, n12, n34, ipqr
c     write (iout,2001) ml,m1,m2,m3,m4,n3,nm3,n4,nm4,n1,n2,n12,n34,ipqr
c2001 format(' ml   m1   m2   m3   m4     ',1x,5i6/
c    1       ' n3   nm3  n4   nm4  n1   n2',1x,6i6/
c    2       ' n12  n34  ipqr             ',1x,3i6)
      lla = ltab(m1+1) + 1
      call mccgc(n1,n3,nm3,cgx(lla),nf36,r,g1(ktab(m1+1)+1),'abix',
     $     bufix,lbufso,rabix,incor)
c
c
c--------------------------------------------c
c     form density matrix block pointers     c
c--------------------------------------------c
 2200 continue
c
c-------------------------------------------------------c
c                                                       c
c     calculate two-electron contributions to c*g       c
c                                                       c
c-------------------------------------------------------c
c
 2500 continue
 3000 continue
      call iosys('read integer "header abcx" from rwf',14,ihd,0,' ')
      if(ihd(1).eq.-100)go to 4000
      ml=ihd(1)
      m1=ihd(2)
      m2=ihd(3)
      m3=ihd(4)
      m4=ihd(5)
      n1=ihd(6)
      nm1=ihd(7)
      n2=ihd(8)
      nm2=ihd(9)
      n3=ihd(10)
      n4=ihd(11)
      n12=ihd(12)
      n34=ihd(13)
      l1234=ihd(14)
c
c     write (iout,2501) nf35
c2501 format(' *mccgb nf35 ',1x,i6)
c
c3000  read (nf35,end=4000) ml, m1, m2, m3, m4, n1, nm1, n2, nm2,
c    1      n3, n4, n12, n34, l1234
c     write (iout,3001) ml,m1,m2,m3,m4,n1,nm1,n2,nm2,n3,n4,n12,n34,l1234
 3001 format(' ml   m1   m2   m3   m4     ',1x,5i6,/
     $     ' n1   nm1  n2   nm2  n3   n4',1x,6i6,/
     $     ' n12  n34  l1234            ',1x,3i6)
      if (n12*n34 .le. 0) go to 3000
c
c------------------
c      (aa aa) case
c------------------
      lla = ltab(m1+1) + ncob(m1+1) * nob(m1+1) + 1
c
      n1234=n12*(n12+1)/2
      do 19001 i=1,n1234
         tg(i)=g2(i)
19001 continue
c
      call mccg1(n1,nm1,cgx(lla),nf35,r,tg,bufix,lbufso,sg,
     $     rabcx,incor)
c52       go to 3000
c
 4000 continue
c
c51      call srew(nf36)
c52      call srew(nf35)
c
c--------------------------c
c     pack c*g into cg     c
c--------------------------c
c
      llc = 1
      do 600 l = 1, nsym
         nocc = ncob(l) + naob(l)
         if (nocc .eq. 0) go to 600
         lla = ltab(l) + 1
         llb = locsym(l) + 1
c
         call mccgpk(nob(l),nocc,cgx(lla),lok(llb),len(llb),mix,
     $        cg(llc),ll)
         llc = llc + ll
 600  continue
c
      return
      end
