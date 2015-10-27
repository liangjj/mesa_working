*deck @(#)trn2e.f	5.1  11/6/94
      subroutine trn2e(values,nnp,ntriang,c,nbf,norbs,t1,t2,numij,
     #                 asort,lnsort,ngroup,
     #                 nmax,nsym,ijgrp,ijadd,kadd,ladd,orbsym,h,
     #                 levfrm,val,lab,
     #                 bin,lenbin,toguga,moval,ops)
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
c***source             @(#)trn2e.f	5.1   11/6/94
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
      logical toguga
      real*8 moval(numij,ntriang)
      real*8 values(nnp,ntriang),c(nbf,norbs),t1(norbs,norbs)
      real*8 t2(norbs,nbf),asort(*),h(numij),sqrt2
      real*8 val(lenbin)
      integer lab(lenbin),bin(lenbin)
      integer ijgrp(numij),ijadd(numij),kadd(norbs,nsym)
      integer ladd(norbs,nsym),orbsym(norbs)
c
      common /io/     itape5,itape6
c
      sqrt2=sqrt(2.0d+00)
c
c
      call iosys('read real "mo one-electron integrals" from scr',
     $            numij,h,0,' ')
c
c
c     ----- set up labels -----
c
      do 50 i=1,numij
         lab(i)=(i-1)*nnp
   50 continue
c
c     ----- loop through kl triangles of integrals, transforming -----
c
      call iosys('rewind "mo two-electron integrals" on scr',
     $            0,0,0,' ')
c
c     ----- initialize the sort to guga order -----
c
      junk=max(33000,nnp**2/10)
      junk=min(junk,4000000)
      call sorter('start',asort,asort,lnsort,ngroup*nmax,1024,
     #             0,0,0,'guga integrals','gints',.false.)
c
c
      pt=0
c
c     ----- loop through ij triangles of integrals, transforming -----
c
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
                call iosys('read real "mo two-electron integrals"'//
     #                     ' from scr without rewinding',lnread,
     $                       values,0,0)
            end if
c
            call trtosq(t1,values(1,ij-minij+1),norbs,numij)
c
c
c              ----- pass the transformed integrals to the sort -----
c
            if (i.ne.j) then
c
c                ----- [ik;jl] & [ij;jl] & [ik;il] & [il;jl] & [il;il]
                do 6 k=j,i
                   ija=ia+k
                   ijs=xor(orbsym(i),orbsym(k))
                   ijk=(ijgrp(ija)-1)*nmax+ijadd(ija)+kadd(j,ijs+1)+1
                   ijks=xor(ijs,orbsym(j))
c
                   if (k.eq.j) then
                       maxl=j-1
                   else
                       maxl=j
                   end if
c
                   if (pt+maxl.gt.lenbin) then
                       call sorter('with bin',asort,asort,0,pt,lab,
     $                              bin,val,0,0,0,.false.)
                       pt=0
                   end if
c
                   do 5 l=1,maxl
                      if (orbsym(l).eq.ijks) then
                          val(pt+1)=t1(k,l)
                          lab(pt+1)=ijk+ladd(l,ijks+1)
                          pt=pt+1
                      end if
 5                 continue
 6              continue
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
     $                              bin,val,0,0,0,.false.)
                       pt=0
                   end if
                   do 7 l=1,k
                      if (orbsym(l).eq.ijks) then
                          val(pt+1)=t1(k,l)
                          lab(pt+1)=ijk+ladd(l,ijks+1)
                          pt=pt+1
                      end if
 7                 continue
 8              continue
c
c                 ----- [il;jk] & [il;jj] -----
c
                do 10 k=j+1,i-1
                   ija=ia+k
                   ijl=(ijgrp(ija)-1)*nmax+ijadd(ija)+
     #                  ladd(j,orbsym(j)+1)+3
                   ijs=xor(orbsym(i),orbsym(k))
                   ijls=xor(ijs,orbsym(j))
                   if (pt+k-j-1.gt.lenbin) then
                       call sorter('with bin',asort,asort,0,pt,lab,
     $                              bin,val,0,0,0,.false.)
                       pt=0
                   end if
                   do 9 l=j+1,k-1
                      if (orbsym(l).eq.ijls) then
                          val(pt+1)=t1(k,l)
                          lab(pt+1)=ijl+kadd(l,ijs+1)
                          pt=pt+1
                      end if
 9                 continue
                   if (orbsym(k).eq.ijls) then
                       if (pt+1.gt.lenbin) then
                           call sorter('with bin',asort,asort,0,pt,
     $                                  lab,bin,val,0,0,0,.false.)
                           pt=0
                       end if
                        val(pt+1)=t1(k,k)
                        lab(pt+1)=ijl+kadd(k,ijs+1)-1
                        pt=pt+1
                   end if
 10             continue
c
c                 ----- [il;ll] -----
c
                if (orbsym(i).eq.orbsym(j)) then
                    if (pt+1.gt.lenbin) then
                        call sorter('with bin',asort,asort,0,pt,lab,
     $                               bin,val,0,0,0,.false.)
                        pt=0
                    end if
                    val(pt+1)=t1(j,j)
                    lab(pt+1)=(ijgrp((i+1)*i/2)-1)*nmax+
     $                         ijadd((i+1)*i/2)+
     #                         kadd(i,1)+ladd(j,orbsym(i)+1)+2
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
     $                              bin,val,0,0,0,.false.)
                       pt=0
                   end if
                   do 11 l=1,k
                      if (orbsym(l).eq.ijks) then
                          val(pt+1)=t1(k,l)
                          lab(pt+1)=ijk+ladd(l,ijks+1)
                          pt=pt+1
                      end if
 11                continue
 12             continue
c
c                 ----- [ii;il] & [ii;ii] -----
c
                ijk=pij+kadd(i,1)+1
                ijks=orbsym(i)
                if (pt+2*i.gt.lenbin) then
                    call sorter('with bin',asort,asort,0,pt,lab,bin,
     #                           val,0,0,0,.false.)
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
 13             continue
            end if
c
c           ----- write out transformed integrals if not sorting
c
c
 30      continue
 31   continue
c
c
c        ----- flush last bin -----
c
      call sorter('with bin',asort,asort,0,pt,lab,bin,
     #             val,0,0,0,.false.)
c
c        ----- complete sort to guga order -----
c
      call sorter('end',asort,asort,0,0,0,0,0,0,0,0,.false.)
c
c
c
      return
      end
