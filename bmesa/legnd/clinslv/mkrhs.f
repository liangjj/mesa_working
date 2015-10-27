*deck mkrhs.f
c***begin prologue     mkrhs
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            right hand side for inhomogeneous, time-dependent hamiltonian.
c***                   
c***references         
c
c***routines called    
c***end prologue       mkrhs
      subroutine mkrhs(p,q,wt,rhs,n)
      implicit integer (a-z)
      logical dollar
      real*8 p, q, wt
      complex*16 rhs
      character*80 chrkey, cpass, card, title
      character*24 type
      dimension p(n,*), q(n), wt(n), rhs(n) 
      common/io/inp, iout
      call czero(rhs,n)
      if ( dollar('$rhs',card,cpass,inp) ) then      
         type=chrkey(card,'type-right-hand-side','one',' ')
      endif
      if(type.eq.'one') then
         do 10 i=1,n
            do 20 j=1,n
               rhs(i) = rhs(i) +wt(j)*p(j,i)
 20         continue
 10      continue   
      elseif(type.eq.'q') then
         do 30 i=1,n
            do 40 j=1,n
               rhs(i) = rhs(i) +wt(j)*q(j)*p(j,i)
 40         continue
 30      continue   
      elseif(type.eq.'exponential') then
         do 50 i=1,n
            do 60 j=1,n
               rhs(i) = rhs(i) +wt(j)*exp(-q(j))*p(j,i)
 60         continue
 50      continue   
      elseif(type.eq.'damped-sine') then
         do 70 i=1,n
            do 80 j=1,n
               rhs(i) = rhs(i) +wt(j)*exp(-q(j))*sin(q(j))*p(j,i)
 80         continue
 70      continue   
      else
         call lnkerr('error in rhs type')
      endif
      title='type right hand side = '//type
      call prntcm(title,rhs,n,1,n,1,iout)
      return
      end       


