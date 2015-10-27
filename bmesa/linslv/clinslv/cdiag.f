*deck cdiag.f
c***begin prologue     cdiag
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diagonalization
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for complex diagonalization.
c***                   
c***references         
c
c***routines called    
c***end prologue       cdiag
      subroutine cdiag(ham,ur,ul,eig,work,n)
      implicit integer (a-z)
      complex*16 ham, ur, ul, eig, cdotc, ovlp, ovlpc
      real*8 work
      character*80 title
      dimension ham(n,n), ur(n,n), ul(n,n), eig(n), work(*)
      common/io/inp, iout 
      call iosys('write real cham to lamdat',2*n*n,ham,0, ' ')
      call cgeev(ham,n,n,eig,ur,n,work,1,info)
      title='eigenvalues of zeroth order hamiltonian'
      call prntcm(title,eig,n,1,n,1,iout)
      call iosys('read real cham from lamdat',2*n*n,ham,0, ' ')
      do 10 i=1,n
         do 20 j=1,n
            ul(i,j) = conjg(ham(j,i))
 20      continue
 10   continue
      call cc2opy(ul,ham,n*n)
      call cgeev(ham,n,n,eig,ul,n,work,1,info)
      title='eigenvalues of zeroth order transposed hamiltonian'
      call prntcm(title,eig,n,1,n,1,iout)
      call iosys('read real cham from lamdat',2*n*n,ham,0, ' ')
      do 30 i=1,n
         ovlp = cdotc(n,ul(1,i),1,ur(1,i),1)
         ovlp=1.d0/sqrt(ovlp)
         ovlpc=conjg(ovlp)
         do 40 j=1,n
            ul(j,i)=ovlpc*ul(j,i)
            ur(j,i)=ovlp*ur(j,i)
 40      continue
 30   continue 
      return
      end       
