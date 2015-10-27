*deck filham.f
c***begin prologue     filham
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            fill hamiltonian blocks in direct in-core
c***                   diagonalization.
c***references         
c
c***routines called    
c***end prologue       filham
      subroutine filham(ham,tr,tth,eigr,i,j,nr,n)
      implicit integer (a-z)
      real*8 ham, tr, tth, eigr
      dimension ham(n,*), tr(nr,nr), eigr(nr)
      common/io/inp, iout
      do 10 k=1,nr
         ham(k,k) = tth/(eigr(k)*eigr(k))
 10   continue
      if(i.eq.j) then
         do 20 k=1,nr
            do 30 l=1,k
               ham(k,l) = ham(k,l) + tr(k,l) 
               ham(l,k)=ham(k,l)
 30         continue
 20      continue
      endif   
      return
      end       
