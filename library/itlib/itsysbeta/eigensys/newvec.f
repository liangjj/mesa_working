*deck newvec.f
c***begin prologue     newvec
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
c***end prologue       newvec
      subroutine newvec(eig,diag,resid,vec,n,m,iter,prnt)
      implicit integer (a-z)
      real*8 eig, diag, resid, vec
      real*8 t1, t2
      real*8 test, zero, nrzero, one
      logical prnt
      character*80 title
      character*5 itoc  
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      dimension eig(n), diag(n), resid(n,m), vec(n,m)
      common/io/inp, iout
      do 10 i=1,m
         do 20 j=1,n
            test=eig(i) - diag(j)
            if(abs(test).ge.nrzero) then
               vec(j,i) = resid(j,i)/test
            else
               vec(j,i) = one
            endif
 20      continue
 10   continue   
      if(prnt) then
         title='new trial vectors iteration = '//itoc(iter)
         call prntrm(title,vec,n,m,n,m,iout)
      endif
      return
      end       





