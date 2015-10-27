c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Cvect0}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck cvect0.f
***begin prologue     cvect0
***date written       960615   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           hamiltonian
***author             schneider, barry (nsf)
***source             3-dim
***purpose            projection of zeroth order eigenvectors on to
***                   dvr basis.
***references         
***routines called    
***end prologue       cvect0
      subroutine cvect0(u01,u02,u03,eig01,eig02,eig03,
     1                 psi0,eig,energy,ind,nd,dim,n,nroot)
      implicit integer (a-z)
      real*8 u01, u02, u03, eig01, eig02, eig03
      real*8 eig, tmp
      real*8 energy
      complex*16 psi0 
      character*16 fptoc
      character*80 title
      dimension nd(3)
      dimension u01(nd(1),nd(1)), u02(nd(2),nd(2)), u03(nd(3),nd(3))
      dimension eig01(nd(1)), eig02(nd(2)), eig03(nd(3))
      dimension psi0(n), ind(n,*), eig(n)
      common/io/inp, iout

c       The eig and ind arrays are destroyed by this routine.

      call czero(psi0,n)         
      energy=eig01(nroot)
      if(dim.eq.1) then
         do 10 i=1,n
            psi0(i) = u01(i,nroot) 
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
               psi0(count) = u01(i,ind(nroot,1)) * 
     1                       u02(j,ind(nroot,2))
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
                  psi0(count) = u01(i,ind(nroot,1)) * 
     1                          u02(j,ind(nroot,2)) * 
     2                          u03(k,ind(nroot,3))
 220           continue
 210        continue
 200     continue
         energy=eig(nroot)
      endif
      write(iout,1) nroot, energy
      return
 1    format(/,5x,'initial state = ',i3,/,5x,
     1            'energy        = ',e15.8)
      end       






