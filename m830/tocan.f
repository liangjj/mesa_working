*deck  @(#)tocan.f	5.1 11/6/94
      subroutine tocan(h,n1int,int,val,lab,bin,lenbin,rsort,
     #     isort,lnsort,isym,jsym,ijsym,nsym,nnpsym,
     #     ijklpt,nnqsym,nso,symoff,kadd,ladd,ijgrp,
     #     ijadd,norbs,numij,bftorb,nbf,nmax,ngroup,
     #     ops,ijpt,inunit,ounit,imngrp,imxgrp,jmngrp,jmxgrp,n2int,
     $     tinfil,toutfi,tout1,orbsym,orbtbf)
c
      implicit integer (a-z)
c
      character*16 inunit,ounit
      character*(*) tinfil,toutfi,tout1
      character*64 infile,outfil,out1
      real*8 h(n1int),int(nmax),val(lenbin),rsort(*)
      integer isort(lnsort)
      integer orbsym(norbs)
      integer orbtbf(norbs)
      integer imngrp(ngroup)
      integer imxgrp(ngroup)
      integer jmngrp(ngroup)
      integer jmxgrp(ngroup)
      integer lab(lenbin),bin(lenbin),isym(nnpsym),jsym(nnpsym)
      integer ijsym(nnpsym),ijklpt(nnqsym),nso(nsym),symoff(nsym)
      integer kadd(norbs,nsym),ladd(norbs,nsym),ijgrp(numij)
      integer ijadd(numij),bftorb(nbf),ijpt(numij)
      character*(*) ops
      real*8 cutoff
      logical logkey
      logical debug
      real*8 f2,f4,f8
c
      parameter (debug=.false.)
c 
      common /io/ inp,iout
c
      ioff(i,j)=i*(i-1)/2+j
c
c     ----- set up the weighting factors: to normalize the density
c           but not integrals (for testing purposes)
c
         f2=1.0d+00/2.0d+00
         f4=1.0d+00/4.0d+00
         f8=1.0d+00/8.0d+00
c..bhl
c      if (logkey(ops,'guga-to-mo=scaling',.true.,' ')) then
c         f2=1.0d+00/2.0d+00
c         f4=1.0d+00/4.0d+00
c         f8=1.0d+00/8.0d+00
c      else
c         f2=1.0d+00
c         f4=1.0d+00
c         f8=1.0d+00
c      end if
c..bhl
c
c     ----- 'cutoff' is the threshold to neglect integrals -----
c
      cutoff=1.0d-12
c
c     ----- set up the symmetry offset array -----
c
      symoff(1)=1
      do 50 i=2,nsym
         symoff(i)=symoff(i-1)+nso(i-1)
 50   continue
c
      call rzero(h,n1int)
c
c     ----- initalize the sorting routines -----
c
      infile=tinfil
      outfil=toutfi
      out1=tout1
c
      n=0
      nints=0
      call sorter('start',isort,rsort,lnsort,n2int,0,0,0,0,
     #             outfil,ounit,.false.)
      call iosys('rewind "'//infile//'" on '//inunit,0,0,0,' ')
c
      do 1000 group=1,ngroup
c
         call iosys('read real "'//infile//'" from '//inunit//
     $        ' without rewinding',nmax,int,0,' ')
c
         if(debug) then 
            write (iout,9001) (int(iq),iq=1,nmax)
 9001       format (5x,'m830: guga density matrices:',/,
     $             (5x,5f12.6))
         end if
c
         do 400 i=imxgrp(group),imngrp(group),-1
            if (i.eq.imxgrp(group)) then
               jmax=jmxgrp(group)
            else
               jmax=i
            end if
            if (i.eq.imngrp(group)) then
               jmin=jmngrp(group)
            else
               jmin=1
            end if
c
            is=orbsym(i)
            ia=orbtbf(i)
            do 300 j=jmax,jmin,-1
               js=orbsym(j)
               ja=orbtbf(j)
               ij=ijadd(ioff(i,j))
               ijs=xor(is,js)
               if (i.eq.j) then
                  ijcas=3
               else
                  ijcas=1
               end if
               do 200 k=j,1,-1
                  ks=orbsym(k)
                  ka=orbtbf(k)
                  ijk=ij+kadd(k,ijs+1)
                  ijks=xor(ijs,ks)
                  if (j.eq.k) then
                     ijkcas=ijcas+1
                  else
                     ijkcas=ijcas
                  end if
                  do 100 l=k,1,-1
c
c                    ----- watch out for ijjj, which is stored as iiij
c
                     if (l.ne.i.and.l.eq.j) go to 100
c
                     ls=orbsym(l)
                     if (ijks.ne.ls) go to 100
c
                     if (n+6.gt.lenbin) then
c
c                       ----- dump this bin to the sort routines -----
c
                        call sorter('with bin',isort,rsort,0,n,
     $                       lab,bin,val,0,0,0,0)
                        nints=nints+n
                        n=0
                     end if
c
                     la=orbtbf(l)
                     ijkl=ijk+ladd(l,ijks+1)
                     if (k.eq.l) then
                        case=ijkcas+4
                     else
                        case=ijkcas
                     end if
c
                     go to (1,2,3,4,5,99,7,8), case
 99                  continue
                        call lnkerr('bad case')
 1                   continue
c
c                       ----- [ij;kl] -----
c
                        call canon(int(ijkl+1)*f8,
     $                       is,ks,js,ls,ia,ka,ja,la,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        call canon(int(ijkl+2)*f8,
     $                       is,js,ks,ls,ia,ja,ka,la,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        call canon(int(ijkl+3)*f8,
     $                       is,ls,js,ks,ia,la,ja,ka,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        go to 10
 2                   continue
c
c                       ----- [ij;jl] -----
c
                        call canon(int(ijkl+1)*f8,
     $                       is,js,js,ls,ia,ja,ja,la,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        call canon(int(ijkl+2)*f4,
     $                       is,ls,js,js,ia,la,ja,ja,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        go to 10
 3                   continue
c
c                       ----- [ii;kl] -----
c
                        call canon(int(ijkl+1)*f8,
     $                       is,ks,is,ls,ia,ka,ia,la,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        call canon(int(ijkl+2)*f4,
     $                       is,is,ks,ls,ia,ia,ka,la,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        go to 10
 4                   continue
c
c                       ----- [ii;il] -----
c
                        call canon(int(ijkl+1)*f4,
     $                       is,is,is,ls,ia,ia,ia,la,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        call canon(int(ijkl+2)*f4,
     $                       is,ls,ls,ls,ia,la,la,la,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        h(ioff(max(ia,la),min(ia,la)))=int(ijkl+3)*f2
                        go to 10
 5                   continue
c
c                       ----- [ij;ll] -----
c
                        call canon(int(ijkl+1)*f8,
     $                       is,ls,js,ls,ia,la,ja,la,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        call canon(int(ijkl+2)*f4,
     $                       is,js,ls,ls,ia,ja,la,la,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        go to 10
 7                   continue
c
c                       ----- [il;il] -----
c
                        call canon(int(ijkl+1)*f4,
     $                       is,ls,is,ls,ia,la,ia,la,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        call canon(int(ijkl+2)*f2,
     $                       is,is,ls,ls,ia,ia,la,la,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        go to 10
 8                   continue
c
c                       ----- [ii;ii] -----
c
                        call canon(int(ijkl+1),is,is,is,is,ia,ia,ia,ia,
     $                       ijklpt,n,val,lab,lenbin,nso,cutoff)
                        h(ioff(ia,ia))=int(ijkl+2)
 10                  continue
 100              continue
 200           continue
 300        continue
 400     continue
 1000 continue
c
c     ----- flush the last bin and finish the sort -----
c
      call sorter('with bin',isort,rsort,0,n,lab,bin,val,0,0,0,0)
      call sorter('end',isort,rsort,0,0,0,0,0,0,0,0,0)
c
      call iosys('write real "'//out1//'" to '//ounit,n1int,h,0,' ')
c
c
      return
      end
