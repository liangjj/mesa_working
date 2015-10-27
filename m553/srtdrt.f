*deck @(#)srtdrt.f	5.1  11/6/94
      subroutine srtdrt(values,nnp,ntriang,c,nbf,norbs,t1,t2,numij,
     $     ngroup,nmax,nsym,ijgrp,ijadd,kadd,ladd,orbsym,h,
     $     levfrm,ijww,klww,ijxx,klxx,nijvir,val,lab,
     $     maxipt,asort,lnsort,itape2,st2,itap91,st91,lenbin,bin)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)srtdrt.f	5.1   11/6/94
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
c
      implicit integer (a-z)
c
      real*8 values(nnp,ntriang),c(nbf,norbs),t1(norbs,norbs)
      real*8 t2(norbs,nbf),h(numij),sqrt2
      real*8 val(2),asort(2),bin(2)
      integer lab(2)
      integer ijgrp(numij),ijadd(numij),kadd(norbs,nsym)
      integer ladd(norbs,nsym),orbsym(norbs)
      integer ijww(numij),klww(nijvir),ijxx(numij),klxx(nijvir)
c
      common /io/ inp,iout
c
c
      sqrt2=sqrt(2.0d+00)
c
c     ----- write the triangles of integrals to mcscr for l912
c
      call iosys('write real "h mcscf" to mcscr',nnp,h,0,' ')
      call iosys('write real "g mcscf" to mcscr',nnp**2,values,0,' ')
c
c     ----- initialize the sort to guga order -----
c
c     write(iout,2238)
 2238 format(/'  initializing sort to guga order ')
cps      call srew(itap91)
cps      call getwad(itap91,st91)
c     write(iout,3333) ngroup,nmax,lnsort
 3333 format(//' in srtdrt- calling sorter, ngroup,nmax,lnsort',
     $     3i7)
cps      call sorter(asort,asort,lnsort,ngroup*nmax,1024,itape2,st2,
cps     #            itap91,st91,0)
      call sorter('start',asort,asort,wpadti(lnsort),ngroup*nmax,
     $            0,0,0,0,'mcscf_ints','ints',.false.)
c
      pt=0
c
c     ----- loop through ij triangles of integrals -----
c
cps          call setwa(itape2,st2)
cinc      maxij=0
c
      ij=0
      do 31 i=1,norbs
         ia=i*(i-1)/2
         do 30 j=1,i
            ja=j*(j-1)/2
            ij=ij+1
c
            call trtosq(t1,values(1,ij),nbf,nnp)
cmp
c
c     ----- pass the transformed integrals to the sort -----
c===============================================================c
            if (i.ne.j) then
c
c     ----- [ik;jl] & [ij;jl] & [ik;il] & [il;jl] & [il;il] -----
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
                     write(iout,2239) pt,maxl,lenbin
 2239                format(/' calling vsort pt,maxl,lenbin :',3i5)
cps                      call vsort(val,lab,bin,pt)
                     call sorter('with bin',asort,asort,0,pt,lab,bin,
     $                    val,0,0,0,.false.)
                     pt=0
                  end if
                  do 5 l=1,maxl
                     if (orbsym(l).eq.ijks) then
                        val(pt+1)=t1(k,l)
                        lab(pt+1)=ijk+ladd(l,ijks+1)
                        pt=pt+1
                     end if
    5             continue
    6          continue
c
c     ----- [ij;kl] & [ij;kk] -----
c
               pij=(ijgrp(ij)-1)*nmax+ijadd(ij)
               ijs=xor(orbsym(i),orbsym(j))
               do 8 k=1,j-1
                  ijk=pij+kadd(k,ijs+1)+2
                  ijks=xor(ijs,orbsym(k))
                  if (pt+k.gt.lenbin) then
c
                     write(iout,2237) pt,k,lenbin
 2237                format(/' calling vsort pt,k,lenbin :',3i5)
cps                      call vsort(val,lab,bin,pt)
                     call sorter('with bin',asort,asort,0,pt,lab,bin,
     $                    val,0,0,0,.false.)
                     pt=0
                  end if
                  do 7 l=1,k
                     if (orbsym(l).eq.ijks) then
                        val(pt+1)=t1(k,l)
                        lab(pt+1)=ijk+ladd(l,ijks+1)
                        pt=pt+1
                     end if
    7             continue
    8          continue
c
c     ----- [il;jk] & [il;jj] -----
c
               do 10 k=j+1,i-1
                  ija=ia+k
                  ijl=(ijgrp(ija)-1)*nmax+ijadd(ija)+
     $                 ladd(j,orbsym(j)+1)+3
                  ijs=xor(orbsym(i),orbsym(k))
                  ijls=xor(ijs,orbsym(j))
                  if (pt+k-j-1.gt.lenbin) then
                     write(iout,2236) pt,k,j,lenbin
 2236                format(/' calling vsort pt,k,j,lenbin :',4i5)
cps                      call vsort(val,lab,bin,pt)
                     call sorter('with bin',asort,asort,0,pt,lab,bin,
     $                    val,0,0,0,.false.)
                     pt=0
                  end if
                  do 9 l=j+1,k-1
                     if (orbsym(l).eq.ijls) then
                        val(pt+1)=t1(k,l)
                        lab(pt+1)=ijl+kadd(l,ijs+1)
                        pt=pt+1
                     end if
    9             continue
                  if (orbsym(k).eq.ijls) then
                     if (pt+1.gt.lenbin) then
                        write(iout,5239) pt,lenbin
 5239                   format(/' calling vsort pt,lenbin :',2i5)
cps                         call vsort(val,lab,bin,pt)
                        call sorter('with bin',asort,asort,0,pt,lab,
     $                       bin,val,0,0,0,.false.)
                        pt=0
                     end if
                     val(pt+1)=t1(k,k)
                     lab(pt+1)=ijl+kadd(k,ijs+1)-1
                     pt=pt+1
                  end if
 10            continue
c
c     ----- [il;ll] -----
c
               if (orbsym(i).eq.orbsym(j)) then
                  if (pt+1.gt.lenbin) then
                     write(iout,4239) pt,lenbin
 4239                format(/' calling vsort pt,lenbin :',2i5)
cps                      call vsort(val,lab,bin,pt)
                     call sorter('with bin',asort,asort,0,pt,lab,bin,
     $                    val,0,0,0,.false.)
                     pt=0
                  end if
                  val(pt+1)=t1(j,j)
                  lab(pt+1)=(ijgrp((i+1)*i/2)-1)*nmax+ijadd((i+1)*i/2)+
     $                 kadd(i,1)+ladd(j,orbsym(i)+1)+2
                  pt=pt+1
               end if
            else
c
c     ----- [ii;kl] & [ii;ll] -----
c
               pij=(ijgrp(ij)-1)*nmax+ijadd(ij)
               do 12 k=1,i-1
                  ijk=pij+kadd(k,1)+2
                  ijks=orbsym(k)
                  if (pt+k.gt.lenbin) then
                     write(iout,7239) pt,k,lenbin
 7239                format(/' calling vsort pt,k,lenbin :',3i5)
cps                      call vsort(val,lab,bin,pt)
                     call sorter('with bin',asort,asort,0,pt,lab,bin,
     $                    val,0,0,0,.false.)
                     pt=0
                  end if
                  do 11 l=1,k
                     if (orbsym(l).eq.ijks) then
                        val(pt+1)=t1(k,l)
                        lab(pt+1)=ijk+ladd(l,ijks+1)
                        pt=pt+1
                     end if
 11               continue
 12            continue
c
c     ----- [ii;il] & [ii;ii] -----
c
               ijk=pij+kadd(i,1)+1
               ijks=orbsym(i)
               if (pt+2*i.gt.lenbin) then
                  write(iout,8239) pt,i,lenbin
 8239             format(/' calling vsort pt,i,lenbin :',3i5)
cps                   call vsort(val,lab,bin,pt)
                  call sorter('with bin',asort,asort,0,pt,lab,bin,val,
     $                 0,0,0,.false.)
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
 13            continue
            end if
 30      continue
 31   continue
c
c     ----- flush last bin -----
c
c         write(iout,8876)
 8876 format(/' flushing last bin ')
cps          call vsort(val,lab,bin,pt)
      call sorter('with bin',asort,asort,0,pt,lab,bin,val,
     $     0,0,0,.false.)
c
c     ----- complete sort to guga order -----
c
cps          call endsrt
      call sorter('end',asort,asort,0,0,0,0,0,0,0,0,.false.)
c
      if(pt.lt.maxipt)go to 9999
c     write(iout,9911) pt,maxipt
 9911 format(/,' mctrans.. routine drtsrt:  pt gt maxipt .. stop ',/,
     $     '  pt  maxipt ',2i8)
      call lnkerr(' ')
 9999 continue
c
c     write(iout,8000) pt
 8000 format(/,' integrals in  drt order  ipt=',i8)
c
c      write(iout,8001)((val(ib),lab(ib)),ib=1,pt)
 8001 format(5(3x,f12.8,1x,i8))
      return
      end
