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
      subroutine vtrial(ham01,cham01,ham02,cham02,ham03,cham03,
     1                  u01,cu01,u02,cu02,u03,cu03,
     2                  eig01,ceig01,eig02,ceig02,eig03,ceig03,
     3                  vec,cvec,eig,ceig,ind,n1,n2,n3,dim,
     2                  n,nroots,prnt,mattyp)
      implicit integer (a-z)
      logical prnt
      character*(*) mattyp
      real*8 ham01, ham02, ham03
      real*8 u01, u02, u03, eig01, eig02, eig03, vec
      real*8 eig, tmp, temp
      complex*16 cham01, cham02, cham03
      complex*16 cu01, cu02, cu03, ceig01, ceig02, ceig03, cvec
      complex*16 ceig, ctmp
      character*80 title
      dimension ham01(n1,n1), ham02(n2,n2), ham03(n3,n3)
      dimension u01(n1,n1), u02(n2,n2), u03(n3,n3)
      dimension eig01(n1), eig02(n2), eig03(n3)
      dimension vec(n,nroots), eig(nroots), ind(n,*)
      dimension cham01(n1,n1), cham02(n2,n2), cham03(n3,n3)
      dimension cu01(n1,n1), cu02(n2,n2), cu03(n3,n3)
      dimension ceig01(n1), ceig02(n2), ceig03(n3)
      dimension cvec(n,nroots), ceig(n)
      common/io/inp, iout
c
c       the ind array is destroyed by this routine and needs to be copied or
c       regenerated if it is used again.
c        
      if(dim.eq.1) then
         if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
            call cc2opy(ceig01,ceig,n)
            do 10 ii=2,n
               i=ii-1
               k=i
               ctmp=ceig(i)
               temp=real(ctmp)
               i1=ind(i,1)
               do 20 j=ii,n
                  if(real(ceig(j)).lt.temp) then
                     k=j
                     ctmp=ceig(j)
                     temp=real(ctmp)
                  endif   
 20            continue
               if(k.ne.i) then
                  ind(i,1)=ind(k,1)
                  ind(k,1)=i1
                  ceig(k) = ceig(i)
                  ceig(i) = ctmp
               endif
 10         continue           
            do 30 i=1,nroots
               do 40 j=1,n
                  cvec(j,i)=cu01(j,ind(i,1))
 40            continue
 30         continue
         else
            call copy(eig01,eig,nroots)
            call copy(u01,vec,n*nroots)
         endif   
      elseif(dim.eq.2) then           
         count=0
         if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
            do 100 i=1,n1
               do 200 j=1,n2
                  count=count+1
                  ceig(count) = ceig01(i) + ceig02(j)
 200           continue   
 100        continue
            do 300 ii=2,n
               i=ii-1
               k=i
               ctmp=ceig(i)
               temp=real(ctmp)
               i1=ind(i,1)
               j1=ind(i,2)
               do 400 j=ii,n
                  if(real(ceig(j)).lt.temp) then
                     k=j
                     ctmp=ceig(j)
                     temp=real(ctmp)
                  endif   
 400           continue
               if(k.ne.i) then
                  ind(i,1)=ind(k,1)
                  ind(i,2)=ind(k,2)
                  ind(k,1)=i1
                  ind(k,2)=j1
                  ceig(k) = ceig(i)
                  ceig(i) = ctmp
               endif
 300        continue                 
            do 500 i=1,nroots
               count=0
               do 600 j=1,n1
                  do 700 k=1,n2
                     count=count+1
                     cvec(count,i)=cu01(j,ind(i,1))*cu02(k,ind(i,2))
 700              continue
 600           continue
 500        continue 
         else          
            do 800 i=1,n1
               do 900 j=1,n2
                  count=count+1
                  eig(count) = eig01(i) + eig02(j)
 900           continue   
 800        continue
            do 1000 ii=2,n
               i=ii-1
               k=i
               tmp=eig(i)
               i1=ind(i,1)
               j1=ind(i,2)
               do 2000 j=ii,n
                  if(eig(j).lt.tmp) then
                     k=j
                     tmp=eig(j)
                  endif   
 2000          continue
               if(k.ne.i) then
                  ind(i,1)=ind(k,1)
                  ind(i,2)=ind(k,2)
                  ind(k,1)=i1
                  ind(k,2)=j1
                  eig(k) = eig(i)
                  eig(i) = tmp
               endif
 1000       continue                 
             do 3000 i=1,nroots
               count=0
               do 4000 j=1,n1
                  do 5000 k=1,n2
                     count=count+1
                     vec(count,i)=u01(j,ind(i,1))*u02(k,ind(i,2))
 5000             continue
 4000          continue
 3000       continue
         endif
      else
         count=0
         if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
            do 6000 i=1,n1
               do 7000 j=1,n2
                  do 8000 k=1,n3
                     count=count+1
                     ceig(count) = ceig01(i) + ceig02(j) + ceig03(k)
 8000             continue
 7000          continue   
 6000       continue
            do 9000 ii=2,n
               i=ii-1
               k=i
               ctmp=ceig(i)
               temp=real(ctmp)
               i1=ind(i,1)
               j1=ind(i,2)
               k1=ind(i,3)
               do 10000 j=ii,n
                  if(real(ceig(j)).lt.temp) then
                     k=j
                     ctmp=ceig(j)
                     temp=real(ctmp)
                  endif   
10000          continue
               if(k.ne.i) then
                  ind(i,1)=ind(k,1)
                  ind(i,2)=ind(k,2)
                  ind(i,3)=ind(k,3)
                  ind(k,1)=i1
                  ind(k,2)=j1
                  ind(k,3)=k1
                  ceig(k) = ceig(i)
                  ceig(i) = ctmp
               endif
 9000       continue   
            do 11000 i=1,nroots
               count=0
               do 12000 j=1,n1
                  do 13000 k=1,n2
                     do 14000 l=1,n3
                        count=count+1
                        cvec(count,i) = cu01(j,ind(i,1)) * 
     1                                  cu02(k,ind(i,2)) * 
     2                                  cu03(l,ind(i,3))
14000                continue     
13000             continue
12000          continue
11000       continue
         else
            do 15000 i=1,n1
               do 16000 j=1,n2
                  do 17000 k=1,n3
                     count=count+1
                     eig(count) = eig01(i) + eig02(j) + eig03(k)
17000             continue
16000          continue   
15000       continue
            do 18000 ii=2,n
               i=ii-1
               k=i
               tmp=eig(i)
               i1=ind(i,1)
               j1=ind(i,2)
               k1=ind(i,3)
               do 19000 j=ii,n
                  if(eig(j).lt.tmp) then
                     k=j
                     tmp=eig(j)
                  endif   
19000          continue
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
18000       continue   
            do 20000 i=1,nroots
               count=0
               do 21000 j=1,n1
                  do 22000 k=1,n2
                     do 23000 l=1,n3
                        count=count+1
                        vec(count,i) = u01(j,ind(i,1)) * 
     1                                 u02(k,ind(i,2)) * 
     2                                 u03(l,ind(i,3))
23000                continue     
22000             continue
21000          continue
20000       continue
         endif     
      endif
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
         if(prnt) then
            title='guess eigenvalues'
            call prntcm(title,ceig,nroots,1,nroots,1,iout)
            title='guess eigenvectors'
            call prntcm(title,cvec,n,nroots,n,nroots,iout)
         endif
      else
         if(prnt) then
            title='guess eigenvalues'
            call prntfm(title,eig,nroots,1,nroots,1,iout)
            title='guess eigenvectors'
            call prntfm(title,vec,n,nroots,n,nroots,iout)
         endif            
      endif                     
      return
      end       






