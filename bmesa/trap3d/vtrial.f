*deck vtrial.f
c***begin prologue     vtrial
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            guess vectors based on separable hamiltonian
c***                   for davidson routine.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vtrial
      subroutine vtrial(u01,u02,u03,eig01,eig02,eig03,vec,eig,
     1                  ind,nd,dim,n,nroots,prnt)
      implicit integer (a-z)
      real*8 u01, u02, u03, eig01, eig02, eig03, vec
      real*8 eig, tmp
      character*80 title
      logical prnt
      dimension nd(3)
      dimension u01(nd(1),nd(1)), u02(nd(2),nd(2)), u03(nd(3),nd(3))
      dimension eig01(nd(1)), eig02(nd(2)), eig03(nd(3))
      dimension vec(n,nroots), eig(nroots), ind(n,*)
      common/io/inp, iout
c
c       the ind array is destroyed by this routine and needs to be copied or
c       regenerated if it is used again.
c        
      if(dim.eq.1) then
         call copy(eig01,eig,nroots)
         call copy(u01,vec,n*nroots)
      elseif(dim.eq.2) then           
         count=0
         do 10 i=1,nd(1)
            do 20 j=1,nd(2)
               count=count+1
               eig(count) = eig01(i) + eig02(j)
 20         continue   
 10      continue
         do 30 ii=2,n
            i=ii-1
            k=i
            tmp=eig(i)
            i1=ind(i,1)
            j1=ind(i,2)
            do 40 j=ii,n
               if(eig(j).lt.tmp) then
                  k=j
                  tmp=eig(j)
               endif   
 40         continue
            if(k.ne.i) then
               ind(i,1)=ind(k,1)
               ind(i,2)=ind(k,2)
               ind(k,1)=i1
               ind(k,2)=j1
               eig(k) = eig(i)
               eig(i) = tmp
            endif
 30      continue                 
         do 50 i=1,nroots
            count=0
            do 60 j=1,nd(1)
               do 70 k=1,nd(2)
                  count=count+1
                  vec(count,i)=u01(j,ind(i,1))*u02(k,ind(i,2))
 70            continue
 60         continue
 50      continue
      else
         count=0
         do 80 i=1,nd(1)
            do 90 j=1, nd(2)
               do 100 k=1,nd(3)
                  count=count+1
                  eig(count) = eig01(i) + eig02(j) + eig03(k)
 100           continue
 90         continue   
 80      continue
         do 200 ii=2,n
            i=ii-1
            k=i
            tmp=eig(i)
            i1=ind(i,1)
            j1=ind(i,2)
            k1=ind(i,3)
            do 210 j=ii,n
               if(eig(j).lt.tmp) then
                  k=j
                  tmp=eig(j)
               endif   
 210        continue
            if(k.ne.i) then
               ind(i,1)=ind(k,1)
               ind(i,2)=ind(k,2)
               ind(i,3)=ind(k,3)
               ind(k,1)=i1
               ind(k,2)=j1
               ind(k,3)=k1
               eig(k) = eig(i)
               eig(i) = tmp
            endif
 200     continue   
         do 300 i=1,nroots
            count=0
            do 310 j=1,nd(1)
               do 320 k=1,nd(2)
                  do 330 l=1,nd(3)
                     count=count+1
                     vec(count,i) = u01(j,ind(i,1)) * 
     1                              u02(k,ind(i,2)) * 
     2                              u03(l,ind(i,3))
 330              continue     
 320           continue
 310        continue
 300     continue   
      endif
      if(prnt) then
         title='guess eigenvalues'
         call prntfm(title,eig,nroots,1,nroots,1,iout)
         title='guess eigenvectors'
         call prntrm(title,vec,n,nroots,n,nroots,iout)
      endif                     
      return
      end       






