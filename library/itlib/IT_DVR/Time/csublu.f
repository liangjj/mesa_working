*deck csublu.f
c***begin prologue     csublu
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            lu factor a complex hmiltonian and solve a set
c***                   complex linear equations.  
c***references         
c
c***routines called    
c***end prologue       csublu
      subroutine csublu(h,tmp,ilst,n,prn,trial)
      implicit integer (a-z)
      complex*16 h, tmp
      character*80 title 
      logical prn, trial
      dimension h(n,n), tmp(n), ilst(n)
      common/io/inp, iout
      title='input matrix'
      call prntcm(title,h,n,n,n,n,iout)
      call cgetrf(n,n,h,n,ilst,info)      
      if(info.ne.0) then
         call lnkerr('error in complex linear solve')
      endif
      if(prn) then 
         title='LU decomposition'
         call prntcm(title,h,n,n,n,n,iout)
      endif
      if(trial) then
         call cgetrs('n',n,1,h,n,ilst,tmp,n,info)
      endif
      return
      end       

