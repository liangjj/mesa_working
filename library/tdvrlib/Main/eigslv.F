*deck eigslv.f
c***begin prologue     eigslv
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
***end prologue       eigslv
      subroutine eigslv(matrix,eigr,eigi,vecl,vecr,n)
      implicit integer (a-z)
      character*80 title
      real*8 matrix, eigr, eigi, vecl, vecr, work, tempr, tempi
      dimension matrix(n,*), eigr(n), eigi(n), vecl(n,*), vecr(n,*)
      common/io/inp, iout
      pointer(pwork,work(1))
      lwork=5*4*n
      need=wpadti(lwork+lwork)
      call getmem(need,pwork,ngot,'work',0)
      call dgeev('v','v',n,matrix,n,eigr,eigi,vecl,n,vecr,n,work,lwork,
     1           info)
      if(info.ne.0) then
         call lnkerr('error from direct diagonalization routine')
      endif
      do 10 ii=2,n
         i=ii-1
         k=i
         tempr=eigr(i)
         tempi=eigi(i)
         do 20 j=ii,n
            if(eigr(j).lt.tempr) then
               k=j
               tempr=eigr(j)
               tempi=eigi(j)
            endif   
 20      continue
         if(k.ne.i) then
            eigr(k) = eigr(i)
            eigr(i) = tempr
            eigi(k) = eigi(i)
            eigi(i) = tempi
            call sswap (n,vecl(1,i),1,vecl(1,k),1)
            call sswap (n,vecr(1,i),1,vecr(1,k),1)
         endif
 10   continue                
      write(iout,1)
      do 30 i=1,n
         write(iout,2) eigr(i), eigi(i)
 30   continue   
      call getmem(-ngot,pwork,idum,'work',idum)
      return
 1    format(/,1x,'eigenvalues of full hamiltonian',
     1       /,10x,'   real part   ',10x,'imaginary part')
 2    format(10x,e15.8,10x,e15.8)
      end       









