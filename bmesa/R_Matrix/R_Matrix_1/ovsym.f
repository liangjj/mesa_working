*deck ovsym.f
c***begin prologue     ovsym
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            projection of perturbed symmetrized states
c***                   onto unperturbed states.
c***                   
c***references         
c
c***routines called    
c***end prologue       ovsym
      subroutine ovsym(ov,scr,vec,vec0,n1,n,prn)
      implicit integer (a-z)
      real*8 ov, scr, vec, vec0, fac
      character*80 title
      logical prn
      dimension ov(n1,n1,n), scr(n1,n1,n), vec(n,n)
      dimension vec0(n1,n1)
      common/io/inp, iout
      fac=1.d0/sqrt(2.d0)
c
c     fill the matrix.
c           <a(1)b(2) | Psi(q) >
c
c     since the basis states of Psi are symmetric, there are identical 
c     contributions  coming from switching the a and b orbitals.  Note also 
c     that the projection on to the product basis of h0 states will 
c     automatically select the symmetric component since Psi(q) is symmetric.
c
      call rzero(scr,n1*n1*n)
      do 10 k=1,n1
         do 20 i=1,k-1
            ind= k*(k-1)/2 + i
            do 30 t=1,n1
               do 40 q=1,n
                  scr(k,t,q) = scr(k,t,q) + vec(ind,q)*fac*vec0(i,t)
 40            continue
 30         continue   
 20      continue      
         ind=k*(k+1)/2
         do 50 t=1,n1
            do 60 q=1,n
               scr(k,t,q) = scr(k,t,q) + vec(ind,q)*vec0(k,t)
 60         continue
 50      continue   
         do 70 i=k+1,n1
            ind=i*(i-1)/2 +k   
            do 80 t=1,n1
               do 90 q=1,n
                  scr(k,t,q) = scr(k,t,q) + vec(ind,q)*fac*vec0(i,t)
 90            continue
 80         continue   
 70      continue   
 10   continue        
      nprd=n1*n 
      call ebtc(ov,vec0,scr,n1,n1,nprd)
      if(prn) then
         title='< psi0(i) psi0(j)| psi(q) >'
         do 100 i=1,n
            call prntfm(title,ov(1,1,i),n1,n1,n1,n1,iout)
 100     continue
      endif          
      return      
      end       






