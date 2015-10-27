*deck eigord.f
c***begin prologue     eigord
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            order eigenvalues and corresponding eigenvectors.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       eigord
      subroutine eigord(eig,reig,imeig,ceig,vec,vecl,vecr,cvecl,cvecr,
     1                  n,dim,type)
      implicit integer (a-z)
      character*(*) type
      character*80 title
      real*8 eig, reig, imeig, vec, vecl, vecr, tempr, tempi
      complex*16 ceig, cvecl, cvecr, ctmp
      dimension eig(n), reig(n), imeig(n), vec(dim,n)
      dimension vecl(dim,*), vecr(dim,*), ceig(n), cvecl(dim,*)
      dimension cvecr(dim,*)
      common/io/inp, iout
      if(type.eq.'complex') then 
         do 10 ii=2,n
            i=ii-1
            k=i
            ctmp=ceig(i)
            tempr=real(ctmp)
            do 20 j=ii,n
               if(real(ceig(j)).lt.tempr) then
                  k=j
                  ctmp=ceig(j)
                  tempr=real(ctmp)
               endif   
 20         continue
            if(k.ne.i) then
               ceig(k) = ceig(i)
               ceig(i) = ctmp
               call cswap(n,cvecl(1,i),1,cvecl(1,k),1)
               call cswap(n,cvecr(1,i),1,cvecr(1,k),1)
            endif
 10      continue           
      elseif(type.eq.'real-symmetric') then
         do 30 ii=2,n
            i=ii-1
            k=i
            tempr=eig(i)
            do 40 j=ii,n
               if(eig(j).lt.tempr) then
                  k=j
                  tempr=eig(j)
               endif   
 40         continue
            if(k.ne.i) then
               eig(k) = eig(i)
               eig(i) = tempr
               call sswap (n,vec(1,i),1,vec(1,k),1)
            endif
 30      continue
       elseif(type.eq.'real-unsymmetric') then
         do 50 ii=2,n
            i=ii-1
            k=i
            tempr=reig(i)
            tempi=imeig(i)
            do 60 j=ii,n
               if(reig(j).lt.tempr) then
                  k=j
                  tempr=reig(j)
                  tempi=imeig(j)
               endif   
 60         continue
            if(k.ne.i) then
               reig(k) = reig(i)
               reig(i) = tempr
               imeig(k) = imeig(i)
               imeig(i) = tempi
               call sswap (n,vecl(1,i),1,vecl(1,k),1)
               call sswap (n,vecr(1,i),1,vecr(1,k),1)
            endif
 50      continue                
      endif  
      return
      end       

