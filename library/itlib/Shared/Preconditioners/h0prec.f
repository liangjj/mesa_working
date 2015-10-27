*deck h0prec.f
c***begin prologue     h0prec
c***date written       970531   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           davidson, trial, vector
c***author             schneider, barry (nsf)
c***source             
c***purpose            preconditioning using the transformation which
c***                   diagonalizes an approximation to A.
c***
c***description        solve AX = Resid approximately using an
c                      approximation for A in its diagonal representation.
c***references         
c
c***routines called    
c***end prologue       h0prec
      subroutine h0prec(energy,diag,resid,vec,eig01,eig02,eig03,u01,u02,
     1                  u03,t1,t2,dim,n1,n2,n3,n,m,drctv,iter,prnt)
      implicit integer (a-z)
      real*8 energy, diag, resid, vec, eig01, eig02, eig03
      real*8 u01, u02, u03, t1, t2
      real*8 test, zero, nrzero, one
      logical drctv, prnt
      character*4 itoc
      character*80 title  
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      dimension diag(n), resid(n,m), vec(n,m)
      dimension eig01(n1), eig02(n2), eig03(n3)
      dimension u01(n1,n1), u02(n2,n2), u03(n3,n3)
      dimension t1(*), t2(*)
      common/io/inp, iout
c
c         transform residual from dvr to h0 representation
c 
      call dvr2h0(resid,vec,u01,u02,u03,t1,t2,n1,n2,n3,n,m,dim)

c
c         using the new and diagonal representation solve the residual 
c         equation approximately
c
      if(dim.eq.1) then
         do 10 i=1,m
            do 20 j=1,n1
               test = energy - eig01(j)
               if(abs(test).ge.nrzero) then
                   vec(j,i) = vec(j,i)/test
               else
                   vec(j,i) = one
               endif
 20         continue
 10      continue
      elseif(dim.eq.2) then   
         do 30 i=1,m
            count=0
            do 40 j=1,n1
               do 50 k=1,n2
                  count=count+1
                  test = energy - eig01(j) - eig02(k)
                  if(abs(test).ge.nrzero) then
                     vec(count,i) = vec(count,i)/test
                  else
                     vec(count,i) = one
                  endif
 50            continue
 40         continue   
 30      continue   
      elseif(dim.eq.3) then
         do 60 i=1,m
            count=0
            do 70 j=1,n1
               do 80 k=1,n2
                  do 90 l=1,n3                      
                     count=count+1
                     test = energy - eig01(j) - eig02(k) - eig03(l)
                     vec(count,i) = vec(count,i)/test
 90               continue
 80            continue   
 70         continue
 60      continue   
      endif
      if(energy.eq.0.d0) then
         call vneg(vec,vec,m*n)
      endif
      call h02dvr(resid,vec,u01,u02,u03,t1,t2,n1,n2,n3,n,m,dim)
      call copy(resid,vec,n*m)
      if(prnt) then
         title='new trial vectors iteration = '//itoc(iter)
         call prntrm(title,vec,n,m,n,m,iout)
      endif
      return
      end       





