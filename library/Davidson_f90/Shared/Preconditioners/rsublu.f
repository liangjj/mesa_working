*deck rsublu.f
c***begin prologue     rsublu
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            lu factor a real hmiltonian and solve a set
c***                   complex linear equations.  
c***references         
c
c***routines called    
c***end prologue       rsublu
      subroutine rsublu(h,tmp,ilst,n,prn,trial)
      implicit integer (a-z)
      real*8 h, tmp
      character*80 title 
      logical prn, trial
      dimension h(n,n), tmp(n), ilst(n)
      common/io/inp, iout
c      call sgefa(h,n,n,ilst,info)      
      call sgetrf(n,n,h,n,ilst,info)
      if(info.ne.0) then
         call lnkerr('error in real linear solve')
      endif
      if(prn) then 
         title='LU decomposition'
         call prntcm(title,h,n,n,n,n,iout)
      endif
      if(trial) then
c         call sgesl(h,n,n,ilst,tmp,0)
         call sgetrs('n',n,1,h,n,ilst,tmp,n,info)
      endif
      return
      end       

