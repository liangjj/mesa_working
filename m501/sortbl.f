*deck @(#)sortbl.f	5.1  11/6/94
      subroutine sortbl(c,eig,dipole,ct,symtyp,sympt,nbf,ops)
      dimension c(nbf,nbf),ct(nbf,nbf),dipole(nbf,nbf),eig(nbf)
      dimension eval(2),evec(2,2),temp(3),tv(2,5),origin(3)
      integer sympt(nbf),symtyp(nbf)
      real*8 nucprt
      character*(*) ops
      logical logkey
c
      dimension scr(30,2),scrdip(2,2)
c
      common /io/inp,iout
c
      data thrshv/1.d-3/,thrshe/1.d-3/,thrshd/1.d-5/
      save thrshv,thrshe,thrshd
c
c         write(iout,*)'  sortbl ',thrshe
c
      call iosys('write real "nosym orbitals" to rwf',nbf*nbf,c,
     $    0,' ')
      call iosys('write real "nosym energies" to rwf',nbf,eig,
     $    0,' ')
c
         nnp=nbf*(nbf+1)/2
c
      call iosys('read integer "number of alpha electrons" from rwf',
     $     1,nae,0,' ')
      call iosys('read integer "number of beta electrons" from rwf',
     $     1,nbe,0,' ')
c
      ncore=nbe
      nactiv=nae-nbe
c
      ncore=intkey(ops,'linear=ncore',ncore,' ')
      nactiv=intkey(ops,'linear=nactive',nactiv,' ')
c
      nocc=ncore+nactiv
      nvirt=nbf-nocc
c
      if(logkey(ops,'linear=ints',.false.,' ')) then
      call iosys('read real e2xx from rwf',3,origin,0,' ')
      call iosys('read real e2xx from rwf without rewinding',
     #   1,nucprt,0,' ')
      call iosys('read real e2xx from rwf without rewinding',
     #   nnp,ct,0,' ')
      int=1
      else
      int=0
      end if
c
      call iosys('read integer sympt from rwf',nbf,sympt,
     #   0,' ')
cc
c         write(iout,*)' sympt '
c         write(iout,101)(sympt(i),i=1,nbf)
c  101    format(2x,5i5)
cc
         call trtosq(dipole,ct,nbf,nnp)
c.         write(iout,*)' ao dipole ints '
c.         call matout(dipole,nbf,nbf,nbf,nbf,iout)
c
ccc
c     transform the ao dipole integrals to the mo basis
ccc
c
c.         do 91 i=1,2
c.          do 92 j=1,nbf
c.           scr(j,i)=0.0
c.           do 93 k=1,nbf
c.           scr(j,i)=scr(j,i)+dipole(j,k)*c(k,2+i)
c.  93       continue
c.  92      continue
c.  91     continue
c
c.         write(iout,*)' scf(j,3) and scr(j,4) '
c.         call matout(scr,nbf,2,nbf,2,iout)
c
c.         do 81 i=1,2
c.          do 82 j=1,2
c.            scrdip(i,j)=0.0
c.              do 83 k=1,nbf
c.               scrdip(i,j)=scrdip(i,j)+c(k,i+2)*scr(k,j)
c.  83          continue
c.  82       continue
c.  81      continue
c
c.         write(iout,*)'  mo dipole 2by2 '
c.         call matout(scrdip,2,2,2,2,iout)
c
         call ebc(ct,dipole,c,nbf,nbf,nbf)
c.         write(iout,*)' half-transformed dipole ints '
c.         call matout(ct,nbf,nbf,nbf,nbf,iout)
         call ebtc(dipole,c,ct,nbf,nbf,nbf)
c.         write(iout,*)' mo dipole ints '
c.         call matout(dipole,nbf,nbf,nbf,nbf,iout)
ccc
c     rotate degenerate orbitals
ccc
         ndeg=2
         nbf2=2*nbf
         i=0
c
  1     continue
c
        i=i+1
c
        if(i+1.gt.nbf)go to 2
c
c.        write(iout,*)' eig(i) eig(i+1)  ',eig(i),eig(i+1),i
c
        if(abs(eig(i)-eig(i+1)).gt.thrshe) go to 1
c
         symtyp(i)=-1
         symtyp(i+1)=-1
c
c.         write(iout,*)' degenerate orbitals '
c
         if(int.eq.1) then
c
          temp(1)=dipole(i,i)
          temp(2)=dipole(i+1,i)
          temp(3)=dipole(i+1,i+1)
c
          if(abs(temp(2)).gt.thrshd) then
c
          call givens(ndeg,ndeg,ndeg,temp,tv,eval,evec)
c
c.         write(iout,*)' eval ',eval(1),eval(2)
          call ebc(ct,c(1,i),evec,nbf,ndeg,ndeg)
          call scopy(nbf2,ct,1,c(1,i),1)
c
          end if
         end if
c
         i=i+1
c
        go to 1
c
  2    continue
ccc
c   determine the symmetry of orbitals
ccc
       do 5 i=1,nbf
         j=isamax(nbf,c(1,i),1)
         symtyp(i)=sympt(j)
c.         write(iout,*)' orbital symmetry ',i,symtyp(i)
  5    continue
c
      if(logkey(ops,'linear=d2h',.false.,' ')) then
        id2h=1
      else
        id2h=0
      end if
      if(logkey(ops,'linear=localize',.false.,' ')) then
        localz=1
       if(logkey(ops,'linear=print',.false.,' ')) then
        prtd2h=1
       else
        prtd2h=0
       end if
      else
        localz=0
      end if
c
      thrd2h=fpkey(ops,'linear=thrd2h',thrshv,' ')
c
        write(iout,*)' d2h  ',id2h
c
        call srtv(c,ct,nbf,nbf,symtyp,eig,dipole,id2h,sympt,
     $ localz,prtd2h,thrd2h)
c
c..
c  the following code will gather core orbitals by symmetry
c  then active orbitals by symmetry etc.
c..
c.       if(ncore.ne.0) then
c.        call srtv(c,ct,nbf,ncore,symtyp,eig,dipole)
c.       end if
c.       if(nactiv.ne.0) then
c.        call srtv(c(1,ncore+1),ct(1,ncore+1),nbf,nactiv,symtyp(ncore+1),
c.     #    eig(ncore+1),dipole(ncore+1,1))
c.       end if
c.       if(nvirt.ne.0) then
c.        call srtv(c(1,nocc+1),ct(1,nocc+1),nbf,nvirt,symtyp(nocc+1),
c.     #    eig(nocc+1),dipole(nocc+1,1))
c..     end if
c
       call scopy(nbf,dipole,1,eig,1)
       call scopy(nbf*nbf,ct,1,c,1)
c
c.       call matout(c,nbf,nbf,nbf,nbf,iout)
c
       return
       end
