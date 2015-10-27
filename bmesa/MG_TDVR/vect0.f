*deck vect0.f
c***begin prologue     vect0
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            projection of zeroth order eigenvectors on to
c***                   dvr basis.
c***                   
c***references         
c
c***routines called    
c***end prologue       vect0
      subroutine vect0(u01,u02,u03,eig01,eig02,eig03,psi0,energy,
     1                 ind,nd,dim,n,nroot)
      implicit integer (a-z)
      real*8 u01, u02, u03, eig01, eig02, eig03
      real*8 eig, tmp
      real*8 psi0, energy
      character*16 fptoc
      character*80 title
      dimension nd(3)
      dimension u01(nd(1),nd(1)), u02(nd(2),nd(2)), u03(nd(3),nd(3))
      dimension eig01(nd(1)), eig02(nd(2)), eig03(nd(3))
      dimension psi0(n,2), ind(n,*)
      pointer(p,eig(1))
      common/io/inp, iout
      need=wptoin(n)
      call memory(need,p,ngot,'vect0',0)
c
c       the ind array is destroyed by this routine and needs to be copied or
c       regenerated if it is used again.
c        

      energy=eig01(nroot)
      if(dim.eq.1) then
         do 10 i=1,nd(1)
            psi0(i,1) = u01(i,nroot)
            psi0(i,2) = 0.d0
 10      continue   
      elseif(dim.eq.2) then
         call setnd2(ind,nd(2),nd(1),n)           
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
            tmp=eig(i)
            i1=ind(i,1)
            j1=ind(i,2)
            do 50 j=ii,n
               if(eig(j).lt.tmp) then
                  k=j
                  tmp=eig(j)
               endif   
 50         continue
            if(k.ne.i) then
               ind(i,1) = ind(k,1)
               ind(i,2) = ind(k,2)
               ind(k,1)=i1
               ind(k,2)=j1
               eig(k) = eig(i)
               eig(i) = tmp
            endif
 40      continue
         count=0
         do 60 i=1,nd(1)
            do 70 j=1,nd(2)
               count=count+1
               psi0(count,1)=u01(i,ind(nroot,1)) * 
     1                       u02(j,ind(nroot,2))
               psi0(count,2)= 0.d0
 70         continue
 60      continue
      else
         call setnd3(ind,nd(3),nd(2),nd(1),n)
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
         do 200 i=1,nd(1)
            do 210 j=1,nd(2)
               do 220 k=1,nd(3)
                  count=count+1
                  psi0(count,1) = u01(i,ind(nroot,1)) * 
     1                            u02(j,ind(nroot,2)) * 
     2                            u03(k,ind(nroot,3))
                  psi0(count,2) = 0.d0
 220           continue
 210        continue
 200     continue
         energy=eig(nroot)
      endif
      write(iout,1) nroot, energy
      call memory(-ngot,p,idum,'vect0',idum)
      return
 1    format(/,5x,'initial state = ',i3,/,5x,
     1            'energy        = ',e15.8)
      end       






