*deck rdiag.f
c***begin prologue     rdiag
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diagonalization
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for real diagonalization.
c***                   
c***references         
c
c***routines called    
***end prologue       rdiag
      subroutine rdiag(mat,eig,vec,work,rep,iter,n,m,prnt)
      implicit integer (a-z)
      character*80 title
      character*3 itoc
      real*8 mat, eig, vec, work, rep
      logical prnt
      dimension mat(n,*), eig(*), work(*), vec(n,*)
      common/io/inp, iout
      call dsyev('v','l',m,mat,n,eig,work,5*m,info)
      if(info.ne.0) then
         call lnkerr('error from direct diagonalization routine')
      endif  
      do 10 i=1,m
         do 20 j=1,m
            vec(j,i)=mat(j,i)
 20      continue
 10   continue    
      if(prnt) then
         do 30 i=1,nend
            work(i)=eig(i)+rep
 30      continue   
         title='eigenvalues of small matrix iteration = '
     1          //itoc(iter)
         call prntfm(title,work,nend,1,nend,1,iout)
      endif
      return
      end       

