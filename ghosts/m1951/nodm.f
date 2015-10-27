*deck %W%  %G%
      subroutine nodm(nbf,d,nos,occ,t1,t2,root)
      implicit integer (a-z)
c
      real*8 d(nbf,nbf),nos(nbf,nbf),occ(nbf),t1(nbf,nbf),t2(nbf,nbf)
      character*4 itoc
c
c     ----- read in the natural orbitals and occupations (diagonal dm)
c
      call iosys('read real "no vector '//itoc(root)//'" from rwf',
     #            nbf**2,nos,0,' ')
      call iosys('read real "no occ '//itoc(root)//'" from rwf',
     #            nbf,occ,0,' ')
c
c     ----- put the diagonalized density matrix back into a matrix
c
      call rzero(d,nbf**2)
      do 1 i=1,nbf
         d(i,i)=occ(i)
    1 continue
c
c     ----- and transform to the ao basis -----
c
      call ebc(t1,nos,d,nbf,nbf,nbf)
      call ebct(d,t1,nos,nbf,nbf,nbf)
c
c
      return
      end
