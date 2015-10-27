*deck @(#)mcgvcb.f	5.1  11/6/94
      subroutine mcgvcb(nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,lok,mix,cm,b,
     $     nda1,lda1,nf35,nf36,nf37,ltrb,
     $     g,ipqrs,r,fpq,tpi,
     $     g1,g2,ncore,intape,buf,lbufso,
     $     rabcx,rabix,raibx,incor,sg)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcgvcb.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c                     this program oversees the construction of
c                     the updated one- and two-electron integrals,
c                     and stores them in a sequential dataset in
c                     alchemy transformed integral format.
c
c --- input
c
c     nsym            no of symmetries.
c     nbf(nsym)       number of basis functions.
c     nob(nsym)       number of orbitals.
c     nfob(nsym)      number of frozen orbitals
c     ncob(nsym)      number of core orbitals.
c     naob(nsym)      number of active orbitals.
c     cv(--)          orbitals.
c     npqrs           no of integral blocks.
c     locsym(nsym)    pointers into the arrays len and loc for the
c                     update vectors.
c     len(--)         no of nonzero componenets of update vectors.
c     lok(--)         location of update vectors in arrays mix and cm.
c     mix(--)         indices of nonzero components of update vectors.
c     b(--)           nonzero components of update vectors.
c     nda1            direct access dataset containing the core
c                     fock-matrix.
c     lda1            address of the f-matrix in nda1.
c     nf35            fortran number of dataset containing the
c                     two electron integrals needed in the construction
c                     of the updated integrals.
c     nf36            fortran number of dataset containing the
c                     integrals needed to make core contributions
c                     to the updated integrals.
c     nf37            fortran number of dataset containing the updated
c                     integrals.
c     ltrb            record length for nf37.
c
c --- working storage
c
c     r(x)            x = max(naob(isym)) * max(nob(isym)).
c     g(x)            x = number of active integrals.
c              ****** this array must be zeroed by the calling program.
c     cm(x)           x=sum(ncob(isym)+naob(isym))*nob(isym)
c     ipqrs(npqrs)
c     fpq(x)          x = sum(nbf(isym)*nbf(isym)+nbf(isym))/2
c     tpi(x)          max(nbf(isym)).
c     ihead(14,npqrs)
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
      dimension ntab(21),mtab(21),ltab(21),ktab(21)
      dimension nbf(2),nob(2),ncob(2),naob(2),cv(2)
      real*8 sg(*)
      real*8 buf(lbufso),rabcx(*),rabix(*),raibx(*)
      dimension ipqrs(2),ihd(33),nfob(2)
      dimension locsym(2),len(2),lok(2),mix(2),b(2),g1(2),g2(2)
      dimension g(2),r(2),tpi(2),fpq(2),cm(2)
      logical debug
c
      parameter (debug=.false.)
      common /io/ inp,iout
c
c
      if (debug) then
         write (iout,9876) (nbf(i),i=1,nsym)
         write (iout,9875) (nob(i),i=1,nsym)
         write (iout,9874) (ncob(i),i=1,nsym)
         write (iout,9873) (naob(i),i=1,nsym)
 9876    format('0*mcgvcb  nbf  ',8(1x,i4))
 9875    format('0*mcgvcb  nob  ',8(1x,i4))
 9874    format('0*mcgvcb  ncob ',8(1x,i4))
 9873    format('0*mcgvcb  naob ',8(1x,i4))
      end if
c
c-------------------------------c
c     set up pointer arrays     c
c-------------------------------c
      nobt = 0
      naobt = 0
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
      do 100 l = 1, nsym
         nobt = nobt + nob(l)
         naobt = naobt + naob(l)
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
      ntab(nsym+1) = nvc
c------------------------------c
c     write header on nf37     c
c------------------------------c
c not needed in cray version
c-------------------------------------c
c     expanded the update vectors     c
c-------------------------------------c
      lll = 0
      do 300 l = 1, nsym
         nocc = ncob(l) + naob(l)
         if (nocc .eq. 0) go to 300
         lla = ltab(l) + 1
         llb = locsym(l) + 1
cc
         call mcgupk(nob(l),nocc,ncob(l),b(lll+1),lok(llb),len(llb),mix,
     $        cm(lla),ll)
         lll = lll + ll
cc
 300  continue
c
      ntoter=nob(1)*nocc
      if (debug) then
         write(iout,*)'  b  vec '
         write(iout,70701)(b(i),i=1,ntoter)
70701    format(5(1x,e15.7))
c
         write(iout,*)' update vector '
         call vecout(cm,nob(1),nocc)
      end if
c
c-------------------------------------------------------c
c                                                       c
c     calculate updated one-electron integrals in g     c
c                                                       c
c-------------------------------------------------------c
      npq = mtab(nsym+1)
c
c -----   read fock matrix -----
c
c
      call iosys('read real mc_core_fock from mcscr',npq,fpq,0,' ')
c
      if (debug) then
         write (iout,1900) (fpq(i),i=1,npq)
 1900    format(' *mcgvcb fpq '//4(1x,f16.8))
      end if
c
      do 1000 l = 1, nsym
         if (naob(l) .le. 0) go to 1000
         call mcfmtr(nbf(l),nob(l),ncob(l),naob(l),cv(ntab(l)+1),
     $        fpq(mtab(l)+1),r,tpi)
         lla = locsym(l) + ncob(l) + 1
         llb = ltab(l) + ncob(l) * nob(l) + 1
         call mcg1e(naob(l),nob(l),lok(lla),len(lla),mix,cm(llb),
     $        r,g(ktab(l)+1))
 1000 continue
c----------------------------------------------------------------------c
c     add core contribtuions to the updated one-electron integrals     c
c----------------------------------------------------------------------c
      if (ncobt .le. 0) go to 2500
c
      call iosys('read integer "header abix" from rwf',14,ihd,0,' ')
      call iosys('rewind abix on rwf',0,0,0,' ')
      if(ihd(1).eq.-100)go to 2500
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
      if (debug) then
         write (iout,2001) ml,m1,m2,m3,m4,n3,nm3,n4,nm4,n1,n2,
     $                     n12,n34,ipqr
 2001    format(' ml   m1   m2   m3   m4     ',1x,5i6/
     $          ' n3   nm3  n4   nm4  n1   n1',1x,6i6/
     $          ' n12  n34  ipqr             ',1x,3i6)
      end if
c
      lla = locsym(m3+1) + 1
      llb = ltab(m3+1) + 1
      call mcgcor(n1,n3,nm3,lok(lla),len(lla),mix,cm(llb),nf36,r,
     $     g(ktab(m1+1)+1),'abix',buf,lbufso,rabix,incor)
c
c------------------------------------------------------c
c     write updated one-electron integrals on nf37(ibm)c
c     store in g1 on cray version                      c
c------------------------------------------------------c
c
 2500 continue
c
      nij=ktab(nsym+1)
c
      if (debug) then
          write(iout,*)' mcgvcb: g1 ',nij
          write(iout,70701)(g(i),i=1,nij)
      end if
c
      do 2550 i=1,nij
         g1(i) = g(i)
 2550 continue
c
c------------------------------------------c
c      zero out the portion of g used      c
c      save 1-e ints in g1 on cray         c
c------------------------------------------c
c
      nij = ktab(nsym+1)
      do 2600 i = 1, nij
 2600 g(i) = 0.d0
c--------------------------------------------c
c     form density matrix block pointers     c
c--------------------------------------------c
 2200 continue
c-------------------------------------------------------c
c                                                       c
c     calculate updated two-electron integrals in g     c
c                                                       c
c-------------------------------------------------------c
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
      if (debug) then
         write (iout,3001) ml,m1,m2,m3,m4,n1,nm1,n2,nm2,n3,n4,
     $                     n12,n34,l1234
 3001    format(' mcgvcb: ml   m1   m2   m3   m4     ',1x,5i6/
     $          ' n1   nm1  n2   nm2  n3   n4',1x,6i6/
     $          ' n12  n34  l1234            ',1x,3i6)
      end if
c
      if (n12*n34 .le. 0) go to 3000
c
c------------------
c      (aa aa) case
c------------------
      lla = locsym(m1+1) + ncob(m1+1) + 1
      llb = ltab(m1+1) + ncob(m1+1) * nob(m1+1) + 1
cc       lbhl=ipqrs(l1234)  ......g(lbhl)
      call mcg1(n1,nm1,lok(lla),len(lla),mix,cm(llb),nf35,r,g,
     $     buf,lbufso,rabcx,incor,sg)
c52       go to 3000
c
c------------------------------------------------------c
c     write updated two-electron integrals on nf37(ibm)c
c     square up in g2 on cray                          c
c------------------------------------------------------c
c
 4000 continue
c
      nijnij=nij*nij
      nijij=(nij*(nij+1))/2
c
      if (debug) then
         write(iout,*)' mcgvcb:  gg '
         write(iout,70701)(g(i),i=1,nijij)
      end if
c
      call trtosq(g2,g,nij,nijij)
c
      if (debug) then
         write(iout,*)' g2 '
         write(iout,70701)(g2(i),i=1,nijnij)
      end if
c
c ----- sort integrals into drt order -----
c
c..bhl..unicos
c      call sorti(g,g,ncore,g1,g2,intape)
c..bhl..unicos
      ncbhl=min(ncore,4*nij*(nij+1))
c..bhl..unicos
      call sorti(g,g,ncbhl,g1,g2,intape)
c
c
 4500 continue
c
c      call srew(nf36)
c52      call srew(nf35)
c
c
      return
      end
