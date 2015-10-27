*deck gentrl.f
c***begin prologue     gentrl
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           trial vectors
c***author             schneider, barry (nsf)
c***source             
c***purpose            find trial vectors based on diagonal elements
c***                   
c***references         
c
c***routines called    
c***end prologue       gentrl
      subroutine gentrl(eig,vec,work,ind,n,ntrial)
      implicit integer (a-z)
      real*8 eig, vec, work, tmp
      character*80 title
      dimension work(n), ind(n)
      dimension vec(n,*), eig(n)
      common/io/inp, iout
      do 10 i=1,n
         ind(i)=i
 10   continue   
      call copy(eig,work,n)
      do 20 ii=2,n
         i=ii-1
         k=i
         tmp=work(i)
         i1=ind(i)
         do 30 j=ii,n
            if(work(j).lt.tmp) then
               k=j
               tmp=work(j)
            endif   
 30      continue
         if(k.ne.i) then
            ind(i)=ind(k)
            ind(k)=i1
            work(k) = work(i)
            work(i) = tmp
         endif
 20   continue   
      do 40 i=1,ntrial
         call rzero(vec(1,i),n)
         eig(i)=work(i)
         vec(ind(i),i)=1.d0
 40   continue   
c      title='trial eigenvalues'
c      call prntrm(title,eig,ntrial,1,ntrial,1,iout)
c      title='trial vectors'
c      call prntrm(title,vec,n,ntrial,n,ntrial,iout)
      return
      end       







