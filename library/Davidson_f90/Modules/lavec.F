*deck lavec.f
c***begin prologue     lavec
c***date written       970531   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           davidson, trial, vector
c***author             schneider, barry (nsf)
c***source             
c***purpose            set up new trial vector based on some zeroth
c***                   order model.  this is equivalent to a preconditioning
c***                   of the matrix.
c***references         
c
c***routines called    
c***end prologue       lavec
      subroutine lavec(energy,diag,resid,n,m,iter,prnt)
      implicit integer (a-z)
      real*8 energy, diag, resid
      real*8 test, zero, nrzero, one
      logical prnt
      character*4 itoc
      character*80 title  
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      dimension diag(n), resid(n,m)
      common/io/inp, iout
      do 10 i=1,m
         do 20 j=1,n
            test = energy - diag(j)
            if(abs(test).ge.nrzero) then
               resid(j,i) = resid(j,i)/test
            else
               resid(j,i) = one
            endif
 20      continue
 10   continue   
      if(prnt) then
         title='new trial vectors iteration = '//itoc(iter)
         call prntrm(title,resid,n,m,n,m,iout)
      endif
      return
      end       





