*deck @(#)canon.f	5.1  11/6/94
      subroutine canon(int,is,js,ks,ls,ia,ja,ka,la,ijklpt,n,values,
     $     labels,lenbin,nso,cutoff)
c
c***begin prologue     canon
c***date written       871028   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           "canonical" indexing
c***author             saxe, paul (lanl)
c***source             @(#)canon.f	5.1   11/6/94
c
c***purpose            to put integrals and associated indices in bins
c                      using the symmetry blocked "canonical" order.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       canon
c
      implicit integer (a-z)
c
      integer is
      integer js
      integer ks
      integer ls
      integer ia
      integer ja
      integer ka
      integer la
      integer ijklpt(*)
      integer n
      integer labels(lenbin)
      integer nso(0:*)
      real*8 int
      real*8 values(lenbin)
      real*8 cutoff
c
      ioff(i,j)=i*(i-1)/2+j
      ioff2(i,j,k,l)=ioff(ioff(i,j),ioff(k,l))
c
c     ----- check for negligible integrals -----
c
      if (abs(int).lt.cutoff) return
c
c     ----- form the symmetries into canonical order -----
c
      isym=is
      jsym=js
      ksym=ks
      lsym=ls
      ibf=ia
      jbf=ja
      kbf=ka
      lbf=la
c
      if (isym.lt.jsym) then
         junk=isym
         isym=jsym
         jsym=junk
         junk=ibf
         ibf=jbf
         jbf=junk
      else if (isym.eq.jsym.and.ibf.lt.jbf) then
         junk=ibf
         ibf=jbf
         jbf=junk
      end if
      if (ksym.lt.lsym) then
         junk=ksym
         ksym=lsym
         lsym=junk
         junk=kbf
         kbf=lbf
         lbf=junk
      else if (ksym.eq.lsym.and.kbf.lt.lbf) then
         junk=kbf
         kbf=lbf
         lbf=junk
      end if
c
      ijsym=ioff(isym+1,jsym+1)
      klsym=ioff(ksym+1,lsym+1)
      if (ijsym.lt.klsym) then
         junk=isym
         isym=ksym
         ksym=junk
         junk=ibf
         ibf=kbf
         kbf=junk
         junk=jsym
         jsym=lsym
         lsym=junk
         junk=jbf
         jbf=lbf
         lbf=junk
         junk=ijsym
         ijsym=klsym
         klsym=junk
      else if (ijsym.eq.klsym.and.ioff(ibf,jbf).lt.ioff(kbf,lbf))
     $        then
         junk=ibf
         ibf=kbf
         kbf=junk
         junk=jbf
         jbf=lbf
         lbf=junk
      end if
c
c     ----- check for non-canonical symmetry labels -----
c
      if (jsym.gt.isym.or.lsym.gt.ksym.or.
     $     (isym.eq.ksym.and.lsym.gt.jsym)) then
         call lnkerr('non-canonical symmetry types')
      end if
c
c     ----- and non-canonical orbital labels -----
c
      if ((isym.eq.jsym.and.jbf.gt.ibf).or.
     #     (ksym.eq.lsym.and.lbf.gt.kbf).or.
     #     (isym.eq.ksym.and.jsym.eq.lsym.and.ibf.eq.kbf
     #     .and.lbf.gt.jbf)) then
         call lnkerr('non-canonical orbital types')
      end if
c
c     ----- some integrals must be placed two places -----
c
      if (isym.eq.ksym.and.jsym.eq.lsym) then
c
c        ----- ij;kl --> ij;kl & kl;ij
c
         ijkl=ijklpt(ioff2(isym+1,jsym+1,ksym+1,lsym+1))
         values(n+1)=int
         values(n+2)=int
c
         if (isym.eq.jsym) then
c
c           ----- [ii;ii] -----
c
            labels(n+1)=ijkl+
     #           (ioff(kbf,lbf)-1)*ioff(nso(isym),nso(jsym))+
     #           ioff(ibf,jbf)
            labels(n+2)=ijkl+
     #           (ioff(ibf,jbf)-1)*ioff(nso(ksym),nso(lsym))+
     #           ioff(kbf,lbf)
         else
c
c           ----- [ij;ij] -----
c
            labels(n+1)=ijkl+ibf+nso(isym)*(jbf-1+nso(jsym)*(kbf-1+
     #           nso(ksym)*(lbf-1)))
            labels(n+2)=ijkl+kbf+nso(ksym)*(lbf-1+nso(lsym)*(ibf-1+
     #           nso(isym)*(jbf-1)))
         end if

         n=n+2
      else
c
c        ----- ij;kl --> ij;kl
c
         ijkl=ijklpt(ioff2(isym,jsym,ksym,lsym))
         values(n+1)=int
         if (isym.eq.jsym.and.ksym.eq.lsym) then
c
c           ----- [ii;jj] -----
c
            labels(n+1)=ijkl+ioff(ibf,jbf)+ioff(nso(isym),nso(isym))*
     #           (ioff(kbf,lbf)-1)
         else
c
c           ----- [ij;kl] -----
c
            labels(n+1)=ijkl+ibf+nso(isym)*(jbf-1+nso(jsym)*(kbf-1+
     #           nso(kbf)*(lbf-1)))
         end if
         n=n+1
      end if
c
c
      return
      end
