*deck @(#)grppt1.f	5.1  11/6/94
      subroutine grppt1(nocont,nobf,natoms,ntypes,bftgrp,bftcmp,
     $                  grpsiz,gptij,gptkl,gklsiz,ngrp,nnpgrp,nbf,ops)
c***begin prologue     grppt1.f
c***date written       871116   (yymmdd)
c***revision date      11/6/94
c
c***keywords           ci gradients, mcscf gradients
c***author             saxe, paul (lanl)
c***source             @(#)grppt1.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       grppt1.f
      implicit none
c     --- input variables -----
      integer natoms,ntypes,ngrp,nnpgrp,nbf
c     --- input arrays (unmodified) ---
      integer nocont(natoms,ntypes)
      integer nobf(ntypes)
      character*(*) ops
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer bftgrp(nbf)
      integer bftcmp(nbf)
      integer grpsiz(ngrp)
      integer gptij(nnpgrp)
      integer gptkl(nnpgrp)
      integer gklsiz(nnpgrp)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer grp,bf,atom,type,i,j,ij
      integer n,cmpnt,k,l,kl
      logical logkey
c
      common /io/ inp,iout
c
c     --- work out the group and component of each basis function ---
      grp=0
      bf=0
      do 4 atom=1,natoms
         do 3 type=1,ntypes
            n=nocont(atom,type)
            if (n.gt.0) then
               grp=grp+1
               cmpnt=0
               grpsiz(grp)=n*nobf(type)
               do 2 i=1,n
                  do 1 j=1,nobf(type)
                     cmpnt=cmpnt+1
                     bf=bf+1
                     bftgrp(bf)=grp
                     bftcmp(bf)=cmpnt
 1                continue
 2             continue
            end if
 3       continue
 4    continue
c
      if (ngrp.ne.grp) call lnkerr('ngrp.ne.grp')
      if (nbf.ne.bf) call lnkerr('nbf.ne.bf')
c
c     --- set up sizes of kl pairs of groups, and pointers to them
      gptkl(1)=0
      gklsiz(1)=grpsiz(1)*grpsiz(1)
      kl=1
      do 6 k=2,ngrp
         do 5 l=1,k
            kl=kl+1
            gklsiz(kl)=grpsiz(k)*grpsiz(l)
            gptkl(kl)=gptkl(kl-1)+gklsiz(kl-1)
 5       continue
 6    continue
c
c     --- and ij pair pointers ---
      gptij(1)=0
      ij=1
      do 8 i=2,ngrp
         do 7 j=1,i
            ij=ij+1
            gptij(ij)=gptij(ij-1)+gptkl(ij)*gklsiz(ij-1)
 7       continue
 8    continue
c
c     --- the number of elements in the grouped 2pdm ---
      n=gptij(nnpgrp)+gptkl(nnpgrp)*gklsiz(nnpgrp)+gklsiz(nnpgrp)**2
c
c
      call iosys('write integer "number of grouped 2pdm elements" '//
     $           'to rwf',1,n,0,' ')
      call iosys('write integer "number of momentum groups" to rwf',
     $           1,ngrp,0,' ')
      call iosys('write integer "basis functions to groups" to rwf',
     $            nbf,bftgrp,0,' ')
      call iosys('write integer "basis functions to components" to rwf',
     $            nbf,bftcmp,0,' ')
      call iosys('write integer "momentum group size" to rwf',
     $            ngrp,grpsiz,0,' ')
      call iosys('write integer "group ij pointers" to rwf',
     $            nnpgrp,gptij,0,' ')
      call iosys('write integer "group kl pointers" to rwf',
     $            nnpgrp,gptkl,0,' ')
      call iosys('write integer "group pair blocksize" to rwf',
     $            nnpgrp,gklsiz,0,' ')
c
c
      if (logkey(ops,'print=gradient=pointer-arrays',.false.,' ')) then
         write (iout,9) (bf,bftgrp(bf),bftcmp(bf),bf=1,nbf)
 9       format (/,5x,'basis function   momentum group  component',/,
     $        5x,'--------------   --------------  ---------',/,
     $        (5x,i8,i17,i13))
c
         write (iout,10) (ij,gklsiz(ij),gptij(ij),gptkl(ij),ij=1,nnpgrp)
 10      format (/,5x,'group pair  size   ij pointer  kl pointer',/,
     $        5x,'----------  ----   ---------   ----------',/,
     $        (5x,i6,i9,i9,i14))
         write (iout,11) (i,grpsiz(i),i=1,ngrp)
 11      format (/,5x,'momentum group   number of components',/,
     $        5x,'--------------   --------------------',/,
     $        (5x,i7,i20))
         write (iout,12) n
 12      format (/,5x,'number of grouped 2pdm elements: ',i7)
      end if
c
c
      return
      end
