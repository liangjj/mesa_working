*deck calrhs.f
c***begin prologue     calrhs
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            calculate right hand side.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       calrhs
      subroutine calrhs(p,crhs,rhs,q,wt,n,mattyp,type,prnt)
      implicit integer (a-z)
      real*8 p, rhs, q, wt
      complex*16 crhs
      character*80 title
      character*(*) type
      character*(*) mattyp
      logical prnt
      dimension p(n,n), crhs(n), rhs(n), q(n), wt(n)
      common/io/inp, iout
      write(iout,1) type
 1    format(/,5x,'type right hand side = 'a16)      
      if(mattyp.eq.'complex') then
         call czero(crhs,n)
         if(type.eq.'delta') then
            crhs(n) = p(n,n)
         elseif(type.eq.'one') then
            do 10 i=1,n
               crhs(i) = p(i,i)*wt(i)
 10         continue
         elseif(type.eq.'x') then
            do 20 i=1,n
               crhs(i) = p(i,i)*wt(i)*q(i)
 20         continue
         endif                                     
         title='right hand side'
         call prntcm(title,crhs,n,1,n,1,iout)
      else
         call rzero(rhs,n)
         if(type.eq.'delta') then
            rhs(n) = p(n,n)
         elseif(type.eq.'one') then
            do 30 i=1,n
               rhs(i) = p(i,i)*wt(i)
 30         continue
         elseif(type.eq.'x') then
            do 40 i=1,n
               rhs(i) = p(i,i)*wt(i)*q(i)
 40         continue
         endif              
         title='right hand side'
         call prntrm(title,rhs,n,1,n,1,iout)
      endif            
      return
      end       
