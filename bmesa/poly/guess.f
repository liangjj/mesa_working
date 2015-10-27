*deck guess.f
      subroutine guess(hambuf,eigbuf,vec,eig,work,n,m,nroots,prnt)
      implicit integer (a-z)
      real*8 hambuf, eigbuf, vec, eig, work
      character*80 title
      logical prnt
      dimension hambuf(n,*), eigbuf(*), vec(n,*), eig(*), work(*)
      common/io/inp, iout
      do 10 i=1,m
         do 20 j=1,m
            vec(i,j)=hambuf(i,j)
 20      continue
 10   continue            
      call tred2(n,m,vec,eig,work,vec)
      call tql2(n,m,eig,work,vec,ierr)
      title='guess eigenvalues'
      call prntrm(title,eig,m,1,m,1,iout)
      if(prnt) then
         title='guess eigenvectors'
         call prntrm(title,vec,n,m,n,n,iout)
      endif 
      do 30 i=1,n         
         eigbuf(i)=hambuf(i,i)
         hambuf(i,i)=0.d0
 30   continue
      return
      end       

