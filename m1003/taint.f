*deck @(#)taint.f	5.1  11/6/94
      subroutine taint(nbf,nob,nco,nao,cm,fpq,buf,lbufso,
     1 rabcx,rabix,sg,tg,tpi,r,tac,taa,g1,g2,ndf,incor)
c
cc
cc
      implicit real*8(a-h,o-z)
      real*8 sg(*),tg(*),tac(*),taa(*)
      real*8 buf(lbufso),rabcx(*),rabix(*)
      dimension ihd(33)
      real*8 g1(*),g2(*)
      dimension r(2),tpi(2),fpq(2),cm(2)
      common /io/ inp,iout
c-----------------------------------------------------------------------
c
c --- description     this program oversees the construction of
c                     the updated one- and two-electron integrals,
c                     and stores them in a sequential dataset in
c                     alchemy transformed integral format.
c
c --- input
c
c     nbf(nsym)       number of basis functions.
c     nob(nsym)       number of orbitals.
c     ncob(nsym)      number of core orbitals.
c     naob(nsym)      number of active orbitals.
c     cv(--)          orbitals.
c     npqrs           no of integral blocks.
c     locsym(nsym)    pointers into the arrays len and loc for the
c                     update vectors.
c
c --- working storage
c
c     r(x)            x = max(naob(isym)) * max(nob(isym)).
c     g(x)            x = number of active integrals.
c              ****** this array must be zeroed by the calling program.
c     cm(x)           x=sum(ncob(isym)+naob(isym))*nob(isym)
c     fpq(x)          x = sum(nbf(isym)*nbf(isym)+nbf(isym))/2
c     tpi(x)          max(nbf(isym)).
c     ihead(14,npqrs)
c-----------------------------------------------------------------------
c
c..bhl.bug      npq=nob*nob
c
      nnob=nob*(nob+1)/2
c23456
c..bhl.bug  call iosys('read real mcscf_ao_core_fock from rwf',npq,fpq,0,' ')
c23456
      call iosys('read real mc_core_fock from mcscr',nnob,fpq,0,' ')
c
c       write (iout,1900) (fpq(i),i=1,nnob)
c1900 format(' *mcgvcb fpq '//4(1x,f16.8))
c
c      write(iout,*)'  fpq  '
c      call matout(fpq,nob,nob,nob,nob,iout)
c

      call mcfmtr(nbf,nob,nco,nao,cm,fpq,r,tpi)
      call tag1e(nao,nob,taa,r,g1,ndf,tg)
c
c
 1000  continue
c----------------------------------------------------------------------c
c     add core contribtuions to the updated one-electron integrals     c
c----------------------------------------------------------------------c
      if (nco .le. 0) go to 2500
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
c..bhl
c      write (iout,2001) ml,m1,m2,m3,m4,n3,nm3,n4,nm4,n1,n2,n12,n34,ipqr
c2001  format(' ml   m1   m2   m3   m4     ',1x,5i6/
c     1       ' n3   nm3  n4   nm4  n1   n1',1x,6i6/
c     2       ' n12  n34  ipqr             ',1x,3i6)
c..bhl
c
      call tagcor(n1,n3,nm3,tac,r,
     1      g1,'abix',buf,lbufso,rabix,incor,ndf)
c
c
c------------------------------------------------------c
c     store in g1 on cray version                      c
c------------------------------------------------------c
c
 2500 continue
c
      nij=nao*(nao+1)/2
c
c      write(iout,*)' mcgvcb: g1 ',nij
c      write(iout,70701)(g(i),i=1,nij)
c      do 2550 i=1,nij
c      g1(i) = g(i)
c2550  continue
c
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
c     write (iout,3001) ml,m1,m2,m3,m4,n1,nm1,n2,nm2,n3,n4,n12,n34,l1234
 3001 format(' mcgvcb: ml   m1   m2   m3   m4     ',1x,5i6/
     1       ' n1   nm1  n2   nm2  n3   n4',1x,6i6/
     2       ' n12  n34  l1234            ',1x,3i6)
       if (n12*n34 .le. 0) go to 3000
c
c------------------
c      (aa aa) case
c------------------
c
      call tag1(n1,nm1,taa,r,g2,
     #           buf,lbufso,rabcx,incor,sg,tg,ndf)
c
 4000 continue
c
 4500 continue
c
      return
      end
