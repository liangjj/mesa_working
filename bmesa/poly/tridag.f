*deck tridag.f
c***begin prologue     tridag
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            tridiagonal approximation to error.
c***references         
c
c***routines called    
c***end prologue       tridag
      subroutine tridag(resid,d,ham,root,dum,n)
      implicit integer (a-z)
      real*8 resid, d, ham, root, dum
      dimension resid(n), d(n), ham(n,n), dum(n,3)
      common/io/inp, iout 
      do 10 i=1,n
         dum(i,2) = root-d(i)
 10   continue   
      do 20 i=2,n
         dum(i,1) = - ham(i,i-1)
 20   continue
      do 30 i=1,n-1
         dum(i,3) = - ham(i,i+1)
 30   continue
      call sgtsl(n,dum(1,1),dum(1,2),dum(1,3),resid,info)
      if(info.ne.0) then
         call lnkerr('error in tridiagonal solve')
      endif   
      return
      end       
