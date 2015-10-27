*deck @@(#)toguga.f	5.1  11/6/94
      subroutine toguga(h,n1int,in,lenbuf,val,lab,bin,lenbin,rsort,
     #                  isort,lnsort,isym,jsym,ijsym,nsym,nnpsym,
     #                  ijklpt,nnqsym,nso,symoff,kadd,ladd,ijgrp,
     #                  ijadd,norbs,numij,bftorb,nbf,nmax,ngroup,
     #                  ops,ijpt,intin,intout,calc)
c
      implicit integer (a-z)
c
      character*(*) calc
      character*16 intin,intout
      real*8 h(n1int),in(lenbuf),val(lenbin),rsort(*)
      integer lab(lenbin),bin(lenbin),isym(nnpsym),jsym(nnpsym)
      integer ijsym(nnpsym),ijklpt(nnqsym),nso(nsym),symoff(nsym)
      integer kadd(norbs,nsym),ladd(norbs,nsym),ijgrp(numij)
      integer ijadd(numij),bftorb(nbf),ijpt(numij)
      character*(*) ops
      integer isort(lnsort)
c
      common /io/ inp,iout
c
      ioff(i,j)=i*(i-1)/2+j
c
c     ----- set up the symmetry offset array -----
c
      symoff(1)=1
      do 1 i=2,nsym
         symoff(i)=symoff(i-1)+nso(i-1)
    1 continue
c
c     ----- combine 'ijadd' and 'ijgrp' into a single offset -----
c
      do 2 i=1,numij
         ijadd(i)=ijadd(i)+(ijgrp(i)-1)*nmax
    2 continue
c
c     ----- fetch the one-electron integrals -----
c
      if (calc.eq.'ci') then
         call iosys('read real "mo one-electron integrals" from '
     $        //intin,n1int,h,0,' ')
      else if (calc.eq.'mcscf') then
         call iosys('read real "h mcscf" from '//intin,n1int,h,0,' ')
      end if
c
c
c     ----- initalize the sorting routines-----
c
      call sorter('start',isort,rsort,lnsort,nmax*ngroup,0,0,0,0,
     #             'guga integrals',intout,.false.)
c
      nbin=0
      pt1=1
      do 200 ijs=1,nnpsym
         is=isym(ijs)
         ni=nso(is)
         js=jsym(ijs)
         nj=nso(js)
         if (is.eq.js) then
            nij=ioff(ni,nj)
         else
            nij=ni*nj
         end if
c
c        ----- deal with one-electron integrals -----
c
         if (is.eq.js) then
            call srone(ni,is,bftorb(symoff(is)),h(pt1),ioff(ni,ni),
     #                 ijadd,numij,kadd(1,1),ladd(1,is),norbs,val,lab,
     #                 bin,lenbin,nbin,isort,rsort,lnsort,ops)
            pt1=pt1+ioff(ni,ni)
         end if
c
         do 100 kls=1,ijs
c
c           ----- check for totally symmetric integrals -----
c
            if (ijsym(ijs).ne.ijsym(kls)) go to 100
            ks=isym(kls)
            nk=nso(ks)
            ls=jsym(kls)
            nl=nso(ls)
            if (ks.eq.ls) then
               nkl=ioff(nk,nl)
            else
               nkl=nk*nl
            end if
            if (nij*nkl.eq.0) go to 100
            maxblk=max(ni,nj,nk,nl)
            ijkl=ioff(ioff(is,js),ioff(ks,ls))
            pt=ijklpt(ijkl)
c
c
c            write (0,10) ijkl,is,js,ks,ls,pt,ni,nj,nk,nl
c   10       format (' symmetry block ',i3,5x,4i1,i8,5x,4i3)
c
c
c           ----- divide up input buffer appropriately -----
c
            nin=lenbuf/nij
c
            if (nin.le.0) then
               call lnkerr('allocation error for input buffer')
            end if
c
c           ----- grab the two-electron integrals -----
c
            if (is.eq.js.and.ks.eq.ls.and.is.eq.ls) then
               call sriiii(ni,in,nin,val,lab,bin,lenbin,nbin,
     #                     nij,pt,ops,is,bftorb(symoff(is)),
     #                     ijadd,kadd,ladd,numij,norbs,nsym,
     #                     isort,rsort,lnsort,intin,calc)
            else if (is.eq.js.and.ks.eq.ls) then
               call sriijj(ni,nk,in,nin,val,lab,bin,lenbin,nbin,
     #                     nij,nkl,pt,ops,is,ks,bftorb(symoff(is)),
     #                     bftorb(symoff(ks)),
     #                     ijadd,kadd,ladd,numij,norbs,nsym,
     #                     isort,rsort,lnsort,intin,calc)
            else if (is.eq.ks.and.js.eq.ls) then
               call srijij(ni,nj,in,nin,val,lab,bin,lenbin,nbin,
     #                     nij,pt,ops,is,js,bftorb(symoff(is)),
     #                     bftorb(symoff(js)),
     #                     ijadd,kadd,ladd,numij,norbs,nsym,
     #                     isort,rsort,lnsort,intin,calc)
            else
               call srijkl(ni,nj,nk,nl,in,nin,val,lab,bin,lenbin,nbin,
     #                     nij,nkl,pt,ops,is,js,ks,ls,
     #                     bftorb(symoff(is)),bftorb(symoff(js)),
     #                     bftorb(symoff(ks)),bftorb(symoff(ls)),
     #                     ijadd,kadd,ladd,numij,norbs,nsym,
     #                     isort,rsort,lnsort,intin,calc)
            end if
  100    continue
  200 continue
c
c     ----- flush the last bin and finish the sort -----
c
  933 format (10x,i3,i10,e20.6)
cend
      call sorter('with bin',isort,rsort,0,nbin,lab,bin,val,0,0,0,0)
      call sorter('end',isort,rsort,0,0,0,0,0,0,0,0,0)
c
c      call iosys('read real "guga integrals" from '//intout,-1,rsort,
c     $     0,' ')
c      write (iout,435) (rsort(i),i=1,nmax)
c 435  format (1x,5f12.6)
c
c
      return
      end
