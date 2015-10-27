*deck @(#)expand.f	5.1  11/6/94
      subroutine expand(c,t,nbf,nsmall,set,eigval,root)
      implicit integer(a-z)
      real*8 c(nbf,nbf),t(*),eigval(*)
      character*128 namchk
      character*4 itoc
      common /io/ inp,iout
c
      call rzero(eigval,nbf)
c
c     prepare the checkpoint file for action.
c
      call iosys('read character "checkpoint filename" from rwf',
     $        0,0,0,namchk)
      call iosys('open chk as old',0,0,0,namchk)
c
c     retrieve the old basis parameters from the chk.
c
      call iosys('read integer "number of basis functions" from chk',
     $        1,nbasj,0,' ')
c
      nbnb=nbasj*nbasj
      nb=nbasj
c bhl 7/10/90
      if(set.eq.2) then
         write(iout,*)' reading    no vectors from chk'
         call iosys('read real "no vector '//itoc(root)//'" from chk',
     $           nbnb,t,0,' ')
         call iosys('read real "no occ '//itoc(root)//'" from chk',
     $           nb,eigval,0,' ')
      else if(set.eq.1) then
         write(iout,*)' reading mcscf vectors from chk'
         call iosys('read real "mcscf vector" from chk',nbnb,t,0,' ')
         call iosys('read real "mcscf orbital energies" from chk',
     $           nb,eigval,0,' ')
      else
         write(iout,*)' reading   scf vectors from chk'
         call iosys('read real "scf vector" from chk',nbnb,t,0,' ')
         call iosys('read real "orbital energies" from chk',
     $           nb,eigval,0,' ')
      end if
c
      nsmall=nbf-nbasj
      nbf2=nbf*nbf
      call rzero(c,nbf2)
c
      ix=1
      do 1 i=1,nbasj
         call scopy(nbasj,t(ix),1,c(1,i),1)
         ix=ix+nbasj
  1   continue
      do 2 i=nbasj+1,nbf
         c(i,i)=1.d0
  2   continue
c
      return
      end
