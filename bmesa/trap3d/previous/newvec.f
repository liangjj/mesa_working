*deck newvec.f
c***begin prologue     newvec
c***date written       970531   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           davidson, trial, vector
c***author             schneider, barry (nsf)
c***source             
c***purpose            set up new trail vector based on some zeroth
c***                   order model.  this is equivalent to a preconditioning
c***                   of the matrix.
c***references         
c
c***routines called    
c***end prologue       newvec
      subroutine newvec(eig,diag,resid,vec,eig01,eig02,eig03,u01,u02,
     1                  u03,t1,t2,t3,dim,nd,n,m,drctv)
      implicit integer (a-z)
      real*8 eig, diag, resid, vec, eig01, eig02, eig03, u01, u02, u03
      real*8 t1, t2, t3
      real*8 test, zero, nrzero, one
      logical drctv  
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      dimension nd(3)
      dimension eig(n), diag(n), resid(n,m), vec(n,m)
      dimension eig01(nd(1)), eig02(nd(2)), eig03(nd(3))
      dimension u01(nd(1),nd(1)), u02(nd(2),nd(2)), u03(nd(3),nd(3))
      dimension t1(*), t2(*), t3(*)
      common/io/inp, iout
      if(.not.drctv) then
         do 10 i=1,m
            do 20 j=1,n
               test=eig(i) - diag(j)
c               if(abs(test).ge.nrzero) then
                  vec(j,i) = resid(j,i)/test
c               else
c                  vec(j,i) = one
c               endif
 20         continue
 10      continue   
      else
         if(dim.eq.1) then
c
c         transform residual from dvr to h0 representation
c 
            call dvr2h0(resid,vec,u01,u02,u03,t1,t2,t3,nd,n,m,dim)

c
c         using the new and diagonal representation solve the residual equation approximately
c
            do 30 i=1,m
               do 40 j=1,nd(1)
                  test=eig(i) - eig01(j)
c                             if(abs(test).ge.nrzero) then
                  vec(j,i) = vec(j,i)/test
c                             else
c                             vec(j,i) = one
c                             endif
 40            continue
 30         continue   
c
c           backtransform to dvr representation
c
            call h02dvr(resid,vec,u01,u02,u03,t1,t2,t3,nd,n,m,dim)
            call copy(resid,vec,n*m)
         elseif(dim.eq.2) then      
c
c        transform residual from  dvr to h0 representation
c
            call dvr2h0(resid,vec,u01,u02,u03,t1,t2,t3,nd,n,m,dim)
c
c         using the new and diagonal representation solve the residual equation approximately
c
            do 50 i=1,m
               count=0
               do 60 j=1,nd(1)
                  do 70 k=1,nd(2)
                     count=count+1
                     test=eig(i) - eig01(j) - eig02(k)
c                             if(abs(test).ge.nrzero) then
                     vec(count,i) = vec(count,i)/test
c                             else
c                             vec(count,i) = one
c                             endif
 70               continue
 60            continue   
 50         continue   
c
c          backtransform to dvr representation
c
            call h02dvr(resid,vec,u01,u02,u03,t1,t2,t3,nd,n,m,dim)
            call copy(resid,vec,n*m)
         elseif(dim.eq.3) then
            call lnkerr('not yet implimented')
         endif
      endif
      return
      end       
