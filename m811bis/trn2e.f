*deck @(#)trn2e.f	1.2  7/30/91
      subroutine trn2e(values,nnp,ntriang,c,nbf,norbs,t1,t2,numij,
     #                 asort,lnsort,ngroup,
     #                 nmax,nsym,ijgrp,ijadd,kadd,ladd,orbsym,h,
     #                 levfrm,val,lab,
     #                 bin,lenbin,toguga,moval,tocan,ops)
c
c***begin prologue     trn2e
c***date written       yymmdd   (yymmdd)
c***revision date      880188   (yymmdd)
c
c   13 january 1988    bhl at brl
c        storing mo integrals in canonical order for ci optimization
c
c   14 october 1987    pws at lanl
c        modifying to be able to write out integrals before sorting
c        to guga order.
c
c    3 december 1986  pws at lanl
c        changing iosys opens to character.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)trn2e.f	1.2   7/30/91
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       trn2e
c
c
c***purpose: no-symmetry four index transformation of two-electron
c            integrals.
c
c paul saxe                 21 august 1984                  lanl
c
      implicit integer (a-z)
c
      character*(*) ops
      logical toguga,tocan
      real*8 moval(numij,ntriang)
      real*8 values(nnp,ntriang),c(nbf,norbs),t1(norbs,norbs)
      real*8 t2(norbs,nbf),asort(*),h(numij),sqrt2
      real*8 val(lenbin)
      integer lab(lenbin),bin(lenbin)
      integer ijgrp(numij),ijadd(numij),kadd(norbs,nsym)
      integer ladd(norbs,nsym),orbsym(norbs)
c
      common /io/ inp,iout
c
      sqrt2=sqrt(2.0d+00)
c
c
c     note that the one-electron integrals go onto 'guga integrals'
c     if toguga.
      if(tocan) then
         call iosys('write real "mo one-electron integrals" to tints',
     $        numij,h,0,' ')
      end if
c
c
      junk=max(33000,nnp**2/10)
      junk=min(junk,4000000)
      call iosys('open scratch as scratch on ssd',junk,0,0,' ')
c
c     ----- initialise the sorting routines -----
c
      call sorter('start',asort,asort,lnsort,nnp*numij,1024,0,0,0,
     #            'half','scratch',.true.)
c
c     ----- set up labels -----
c
      do 50 i=1,numij
         lab(i)=(i-1)*nnp
   50 continue
c
c     ----- loop through kl triangles of integrals, transforming -----
c
      call iosys('rewind "sorted ao integrals" on ints',0,0,0,' ')
      maxkl=0
      kl=0
      do 4 k=1,nbf
         do 3 l=1,k
            kl=kl+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
                minkl=maxkl+1
                maxkl=min(nnp,maxkl+ntriang)
                lnread=(maxkl-minkl+1)*nnp
                call iosys('read real "sorted ao integrals" from ints '
     #                    //'without rewinding',lnread,values,0,' ')
            end if
c
c     ----- perform first half-transformation -----
c
            call trtosq(t1,values(1,kl-minkl+1),nbf,nnp)
            call ebtc(t2,c,t1,norbs,nbf,nbf)
            call ebc (t1,t2,c,norbs,nbf,norbs)
c
c     ----- pass the half-transformed intgrals to the sort -----
c
c            klij=kl
c            do 2 i=1,norbs
c               do 1 j=1,i
c                  call sort(t1(i,j),klij)
c                  klij=klij+nnp
c    1          continue
c    2       continue
            do 21 i=1,numij
               lab(i)=lab(i)+1
   21       continue
            ij=0
            do 22 i=1,norbs
               do 23 j=1,i
                  ij=ij+1
                  val(ij)=t1(i,j)
   23          continue
   22       continue
            if(ij.ne.numij) then
                call lnkerr('ij.ne.numij in trn2e')
            endif
c           call sorter('with bin',asort,asort,0,0,0,0,0,0,0,0,
c    #                   val,lab,bin,numij,.true.)
            call sorter('with bin',asort,asort,0,numij,lab,bin,
     #                   val,0,0,0,.true.)
    3    continue
    4 continue
c
c     ----- finish sorting the half-transformed integrals to kl;ij
c            triangles on tape2
c
      call sorter('end',asort,asort,0,0,0,0,0,0,0,0,.true.)
c
c     ----- initialize the sort to guga order -----
c
      if (toguga) then
          junk=max(33000,nnp**2/10)
          junk=min(junk,4000000)
          call sorter('start',asort,asort,lnsort,ngroup*nmax,1024,
     #                 0,0,0,'guga integrals','gints',.true.)
      end if
      if (tocan) then
          call iosys('create real "mo two-electron integrals" on '//
     $               'tints',numij**2,0,0,' ')
      end if
c
      pt=0
c
c     ----- loop through ij triangles of integrals, transforming -----
c
      call iosys('rewind half on scratch',0,0,0,' ')
      maxij=0
      ij=0
      do 31 i=1,norbs
         ia=i*(i-1)/2
         do 30 j=1,i
            ja=j*(j-1)/2
            ij=ij+1
c
c           ----- check that this triangle of integrals is in core ----
c
            if (ij.gt.maxij) then
                minij=maxij+1
                maxij=min(numij,maxij+ntriang)
                lnread=(maxij-minij+1)*nnp
                call iosys('read real half from scratch '//
     #                     'without rewinding',lnread,values,0,' ')
            end if
c
c     ----- perform second half-transformation -----
c
            call trtosq(t1,values(1,ij-minij+1),nbf,nnp)
            call ebtc(t2,c,t1,norbs,nbf,nbf)
            call ebc (t1,t2,c,norbs,nbf,norbs)
c
c           ----- stuff back into arrays if not sorting -----
c
            if (tocan) then
                call sqtotr(moval(1,ij-minij+1),t1,norbs,numij)
            end if
c
            if(toguga) then
c
c              ----- pass the transformed integrals to the sort -----
c
               if (i.ne.j) then
c
c                ----- [ik;jl] & [ij;jl] & [ik;il] & [il;jl] & [il;il]
c
                  do 6 k=j,i
                     ija=ia+k
                     ijs=xor(orbsym(i),orbsym(k))
                     ijk=(ijgrp(ija)-1)*nmax+ijadd(ija)+kadd(j,ijs+1)+1
                     ijks=xor(ijs,orbsym(j))
                     if (k.eq.j) then
                         maxl=j-1
                     else
                         maxl=j
                     end if
                     if (pt+maxl.gt.lenbin) then
                         call sorter('with bin',asort,asort,0,pt,lab,
     $                                bin,val,0,0,0,.true.)
                         pt=0
                     end if
                     do 5 l=1,maxl
                        if (orbsym(l).eq.ijks) then
                            val(pt+1)=t1(k,l)
                            lab(pt+1)=ijk+ladd(l,ijks+1)
                            pt=pt+1
                        end if
 5                   continue
 6                continue
c
c                 ----- [ij;kl] & [ij;kk] -----
c
                  pij=(ijgrp(ij)-1)*nmax+ijadd(ij)
                  ijs=xor(orbsym(i),orbsym(j))
                  do 8 k=1,j-1
                     ijk=pij+kadd(k,ijs+1)+2
                     ijks=xor(ijs,orbsym(k))
                     if (pt+k.gt.lenbin) then
                         call sorter('with bin',asort,asort,0,pt,lab,
     $                                bin,val,0,0,0,.true.)
                         pt=0
                     end if
                     do 7 l=1,k
                        if (orbsym(l).eq.ijks) then
                            val(pt+1)=t1(k,l)
                            lab(pt+1)=ijk+ladd(l,ijks+1)
                            pt=pt+1
                        end if
 7                   continue
 8                continue
c
c                 ----- [il;jk] & [il;jj] -----
c
                  do 10 k=j+1,i-1
                     ija=ia+k
                     ijl=(ijgrp(ija)-1)*nmax+ijadd(ija)+
     #                    ladd(j,orbsym(j)+1)+3
                     ijs=xor(orbsym(i),orbsym(k))
                     ijls=xor(ijs,orbsym(j))
                     if (pt+k-j-1.gt.lenbin) then
                         call sorter('with bin',asort,asort,0,pt,lab,
     $                                bin,val,0,0,0,.true.)
                         pt=0
                     end if
                     do 9 l=j+1,k-1
                        if (orbsym(l).eq.ijls) then
                            val(pt+1)=t1(k,l)
                            lab(pt+1)=ijl+kadd(l,ijs+1)
                            pt=pt+1
                        end if
 9                   continue
                     if (orbsym(k).eq.ijls) then
                         if (pt+1.gt.lenbin) then
                             call sorter('with bin',asort,asort,0,pt,
     $                                    lab,bin,val,0,0,0,.true.)
                             pt=0
                        end if
                        val(pt+1)=t1(k,k)
                        lab(pt+1)=ijl+kadd(k,ijs+1)-1
                        pt=pt+1
                     end if
 10               continue
c
c                 ----- [il;ll] -----
c
                  if (orbsym(i).eq.orbsym(j)) then
                      if (pt+1.gt.lenbin) then
                          call sorter('with bin',asort,asort,0,pt,lab,
     $                                 bin,val,0,0,0,.true.)
                          pt=0
                     end if
                     val(pt+1)=t1(j,j)
                     lab(pt+1)=(ijgrp((i+1)*i/2)-1)*nmax+
     $                    ijadd((i+1)*i/2)+
     #                    kadd(i,1)+ladd(j,orbsym(i)+1)+2
                     pt=pt+1
                  end if
               else
c
c              ----- [ii;kl] & [ii;ll] -----
c
                  pij=(ijgrp(ij)-1)*nmax+ijadd(ij)
                  do 12 k=1,i-1
                     ijk=pij+kadd(k,1)+2
                     ijks=orbsym(k)
                     if (pt+k.gt.lenbin) then
                         call sorter('with bin',asort,asort,0,pt,lab,
     $                                bin,val,0,0,0,.true.)
                         pt=0
                     end if
                     do 11 l=1,k
                        if (orbsym(l).eq.ijks) then
                            val(pt+1)=t1(k,l)
                            lab(pt+1)=ijk+ladd(l,ijks+1)
                            pt=pt+1
                        end if
 11                  continue
 12               continue
c
c                 ----- [ii;il] & [ii;ii] -----
c
                  ijk=pij+kadd(i,1)+1
                  ijks=orbsym(i)
                  if (pt+2*i.gt.lenbin) then
                      call sorter('with bin',asort,asort,0,pt,lab,bin,
     #                             val,0,0,0,.true.)
                      pt=0
                  end if
                  do 13 l=1,i
                     if (orbsym(l).eq.ijks) then
                         ijkl=ijk+ladd(l,ijks+1)
                         val(pt+1)=t1(i,l)
                         lab(pt+1)=ijkl
                         if (l.lt.i) then
                             val(pt+2)=h(ia+l)
                             lab(pt+2)=ijkl+2
                        else
                             val(pt+2)=h(ia+l)
                             lab(pt+2)=ijkl+1
                        end if
                        pt=pt+2
                     end if
 13               continue
               end if
            end if
c
c           ----- write out transformed integrals if not sorting
c
            if (tocan.and.ij.eq.maxij) then
                lnwrit=(maxij-minij+1)*numij
                call iosys('write real "mo two-electron integrals"'//
     $                     ' to tints without rewinding',lnwrit,
     $                       moval,0,' ')
            end if
c
 30      continue
 31   continue
c
      if (toguga) then
c        ----- flush last bin -----
          call sorter('with bin',asort,asort,0,pt,lab,bin,
     #                 val,0,0,0,.true.)
c        ----- complete sort to guga order -----
          call sorter('end',asort,asort,0,0,0,0,0,0,0,0,.true.)
      end if
c
c
c     ----- destroy the scratch files -----
c
      call iosys('destroy scratch',0,0,0,' ')
c
      return
      end
