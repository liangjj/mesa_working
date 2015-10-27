*deck dcslve.f
c***begin prologue     dcslve
c***date written       970531   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           davidson, trial, vector
c***author             schneider, barry (nsf)
c***source             
c***purpose            preconditioning using a diagonal approximation.
c***
c***description        solve the equation AX = Resid approximately
c                      using the diagonal of A.                     
c***references         
c
c***routines called    
c***end prologue       dcslve
      subroutine dcslve(energy,diag,resid,vec,n,m,iter,prnt)
      implicit integer (a-z)
      real*8 energy
      complex*16  diag, resid, vec, test
      real*8 zero, nrzero, one
      logical prnt
      character*4 itoc
      character*80 title  
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      dimension diag(n), resid(n,m), vec(n,m)
      common/io/inp, iout
      do 10 i=1,m
         do 20 j=1,n
            test = energy - diag(j)
            if(abs(test).ge.nrzero) then
               vec(j,i) = resid(j,i)/test
            else
              vec(j,i) = one
            endif
 20      continue
 10   continue   
      if(energy.eq.0.d0) then
         call vneg(vec,vec,2*m*n)
      endif
      if(prnt) then
         title='new trial vectors iteration = '//itoc(iter)
         call prntcm(title,vec,n,m,n,m,iout)
      endif
      return
      end       





