*deck srteig.f
c***begin prologue     srteig
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            eigenvector for composite hamiltonian
c***                   
c***references         
c
c***routines called    
c***end prologue       srteig
      subroutine srteig(u01,u02,u03,eig01,eig02,eig03,eig,vec,t0,ind,
     1                  nd,dim,n,nroot,prnt)
      implicit integer (a-z)
      real*8 u01, u02, u03, eig01, eig02, eig03
      real*8 eig, vec, t0
      logical prnt
      character*80, title
      character*16 fptoc
      dimension nd(dim), ind(n,3)
      dimension u01(nd(1),nd(1)), u02(nd(2),nd(2)), u03(nd(3),nd(3))
      dimension eig01(*), eig02(*), eig03(*), vec(*), eig(*)
      common/io/inp, iout
c
      if(dim.eq.1) then
         call copy(eig01,eig,n)
         do 10 i=1,n
            vec(i) = u01(i,nroot)
 10      continue   
      elseif(dim.eq.2) then           
         count=0
         do 20 i=1,nd(1)
            do 30 j=1,nd(2)
               count=count+ 1
               eig(count) = eig01(i) + eig02(j)
 30         continue   
 20      continue
         do 40 ii=2,n
            i=ii-1
            k=i
            i1=ind(i,1)
            j1=ind(i,2)
            tmp=eig(i)
            do 50 j=ii,n
               if(eig(j).lt.tmp) then
                  k=j
                  tmp=eig(j)
               endif   
 50         continue
            if(k.ne.i) then
               ind(i,1) = ind(k,1)
               ind(i,2) = ind(k,2)
               ind(k,1) = i1
               ind(k,2) = j1            
               eig(k) = eig(i)
               eig(i) = tmp
            endif
 40      continue
         count=0            
         do 60 i=1,nd(1)
            do 70 j=1,nd(2)
               count=count+1
               vec(count) = u01(i,ind(nroot,1))*u02(j,ind(nroot,2))
 70         continue
 60      continue
      else
         count=0
         do 100 i=1,nd(1)
            do 110 j=1,nd(2)
               do 120 k=1,nd(3)
                  count=count+1
                  eig(count) = eig01(i) + eig02(j) + eig03(k)
 120           continue
 110        continue   
 100     continue
         do 130 ii=2,n
            i=ii-1
            k=i
            tmp=eig(i)
            i1=ind(i,1)
            j1=ind(i,2)
            k1=ind(i,3)
            do 140 j=ii,n
               if(eig(j).lt.tmp) then
                  k=j
                  tmp=eig(j)
               endif   
 140        continue
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
 130     continue   
         count=0
         do 200 j=1,nd(1)
            do 210 k=1,nd(2)
               do 220 l=1,nd(3)
                  count=count+1
                  vec(count) = u01(i,ind(nroot,1)) * 
     1                         u02(j,ind(nroot,2)) * 
     2                         u03(k,ind(nroot,3))
 220           continue     
 210        continue
 200     continue
      endif
      if(prnt) then
         title='eigenvector at t0 = '//fptoc(t0)//' energy = '
     1                             //fptoc(eig(nroot))
         call prntrm(title,vec,n,1,n,1,iout)
      endif                     
      return
      end       








