*deck @(#)trn2dm.f	5.1  11/6/94
      subroutine trn2dm(out,nnp,ntriang,c,nbf,t1,t2,asort,lnsort,
     #     nsym,val,lab,bin,lenbin,in,nactiv,nnpact,ncore,ci,mcscf)
c
c***begin prologue     trn2dm
c***date written       871117   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)trn2dm.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       trn2dm
c
c
c***purpose: no-symmetry four index transformation of the 2pdm
c
c paul saxe                 21 august 1984                  lanl
c
      implicit integer (a-z)
c
      character*3 answer
      character*32 fin
      character*32 fout
      logical ci
      logical mcscf
      real*8 in(nnpact,ntriang)
      real*8 out(nnp,ntriang)
      real*8 c(nbf,nbf)
      real*8 t1(nbf,nbf)
      real*8 t2(nbf,nbf)
      real*8 asort(*)
      real*8 val(lenbin)
      integer lab(lenbin)
      integer bin(lenbin)
c
      common /io/     itape5,itape6
c
      if (mcscf) then
         fin='mcscf mo 2pdm'
         fout='mcscf ao 2pdm'
      else
         fin='mo 2pdm'
         fout='ci ao 2pdm'
      end if
c
      junk=max(33000,nnp**2/10)
      junk=min(junk,4000000)
      call iosys('open scratch as scratch on ssd',junk,0,0,' ')
c
c     ----- initialise the sorting routines -----
c
      call sorter('start',asort,asort,lnsort,nnp*nnpact,1024,0,0,0,
     #            'half','scratch',.false.)
c
c     ----- set up labels -----
c
      do 50 i=1,nnp
         lab(i)=(i-1)*nnpact
   50 continue
c
c     ----- loop through kl triangles of the 2pdm, transforming -----
c
      call iosys('rewind "'//fin//'" on moden',0,0,0,' ')
      maxkl=0
      kl=0
      do 4 k=1,nactiv
         do 3 l=1,k
            kl=kl+1
c
c     ----- check that this triangle of the 2pdm is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnpact,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnpact
               call iosys('read real "'//fin//'" from moden '
     #                    //'without rewinding',lnread,in,0,' ')
            end if
c
c     ----- perform first half-transformation -----
c
            call trtosq(t1,in(1,kl-minkl+1),nactiv,nnpact)
            call ebc (t2,c(1,ncore+1),t1,nbf,nactiv,nactiv)
            call ebct(t1,t2,c(1,ncore+1),nbf,nactiv,nbf)
c
c     ----- pass the half-transformed 2pdm to the sort -----
c
            do 21 i=1,nnp
               lab(i)=lab(i)+1
   21       continue
            ij=0
            do 22 i=1,nbf
               do 23 j=1,i
                  ij=ij+1
                  val(ij)=t1(i,j)
   23          continue
   22       continue
c
            call sorter('with bin',asort,asort,0,nnp,lab,bin,
     #                   val,0,0,0,.false.)
    3    continue
    4 continue
c
c     ----- finish sorting the half-transformed 2pdm to kl;ij
c            triangles
c
      call sorter('end',asort,asort,0,0,0,0,0,0,0,0,.false.)
c
c     ----- create the output file -----
c
      call iosys('does "'//fout//'" exist on aoden',0,0,0,answer)
      if (answer.eq.'no') then
         call iosys('create real "'//fout//'" on aoden',nnp**2,0,0,' ')
      end if
      call iosys('rewind all on aoden',0,0,0,' ')
      pt=0
c
c     ----- loop through ij triangles of the 2pdm, transforming -----
c
      call iosys('rewind half on scratch',0,0,0,' ')
      maxij=0
      ij=0
      do 31 i=1,nbf
         ia=i*(i-1)/2
         do 30 j=1,i
            ja=j*(j-1)/2
            ij=ij+1
c
c           ----- check that this triangle of the 2pdm is in core ----
c
            if (ij.gt.maxij) then
               minij=maxij+1
               maxij=min(nnp,maxij+ntriang)
               lnread=(maxij-minij+1)*nnpact
               call iosys('read real half from scratch '//
     #                    'without rewinding',lnread,in,0,' ')
            end if
c
c     ----- perform second half-transformation -----
c
            call trtosq(t1,in(1,ij-minij+1),nactiv,nnpact)
            call ebc (t2,c(1,ncore+1),t1,nbf,nactiv,nactiv)
            call ebct(t1,t2,c(1,ncore+1),nbf,nactiv,nbf)
            call sqtotr(out(1,ij-minij+1),t1,nbf,nnp)
c
c           ----- write out the ao 2pdm elements -----
c
            if (ij.eq.maxij) then
               lnwrit=(maxij-minij+1)*nnp
               call iosys('write real "'//fout//'"'//
     $              ' to aoden without rewinding',lnwrit,out,0,' ')
            end if
c
 30      continue
 31   continue
c
c     ----- destroy the scratch files -----
c

      call iosys('destroy scratch',0,0,0,' ')
c
      return
      end
