*deck @(#)testh.f	5.1  11/6/94
      subroutine testh(hpp,hopt,c0,t,npvec,iu,mdim,iout)
      implicit real*8(a-h,o-z)
      real*8 hpp(npvec,npvec),hopt(npvec,npvec),c0(mdim),
     $          t(npvec)
c
c
      call iosys('read real "ci root 1" from rwf',
     $           mdim,c0,0,' ')
c
c
      do 1 i=1,npvec
         do 2 j=1,npvec
            hpp(i,j)=hpp(i,j)-hopt(i,j)
  2      continue
  1   continue
c
c
      write(iout,*)'  scattering hamiltonian '
      call matout(hpp,npvec,npvec,npvec,npvec,iout)
      call ebc(t,hpp,c0,npvec,npvec,1)
c
      write(iout,*)' p-space eigenvector '
      call matout(c0,npvec,1,npvec,1,iout)
      write(iout,*)' test vector '
      call matout(t,npvec,1,npvec,1,iout)
c
      return
      end
