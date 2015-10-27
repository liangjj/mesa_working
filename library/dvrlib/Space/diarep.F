*deck diarep.f
c***begin prologue     diarep
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            diagonalization of DVR/FEM hamiltonian.
c***                   
c***references         
c
c***routines called    
c***end prologue       diarep
      subroutine diarep(ham,eigval,eigvec,work,n,prn)
      implicit integer (a-z)
      real*8 ham, eigval, eigvec, work
      character*80 title
      logical prn
      dimension ham(n,n), eigval(n), eigvec(n,n)
      dimension work(*), prn(*)
      common/io/inp, iout
      call copy(ham,eigvec,n*n)
      call dsyev('v','l',n,eigvec,n,eigval,work,5*n,info)
      if(prn(1)) then
         title='eigenvalues'
         call prntrm(title,eigval,n,1,n,1,iout)
      endif
      if(prn(2)) then
         title='eigenvectors'
         call prntrm(title,eigvec,n,n,n,n,iout)
      endif
      return
      end       



