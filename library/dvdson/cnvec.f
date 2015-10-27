*deck cnvec.f
c***begin prologue     cnvec
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
c***end prologue       cnvec
      subroutine cnvec(eig,diag,resid,vec,eig01,eig02,eig03,u0l1,u0l2,
     1                  u0l3,u0r1,u0r2,u0r3,t1,t2,dim,n1,n2,n3,n,m,
     2                  drctv,iter,prnt)
      implicit integer (a-z)
      complex*16 eig, diag, resid, vec, eig01, eig02, eig03
      complex*16 u0l1, u0l2, u0l3, u0r1, u0r2, u0r3, t1, t2, test
      real*8 zero, nrzero, one
      logical drctv, prnt
      character*4 itoc
      character*80 title  
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      dimension eig(n), diag(n), resid(n,m), vec(n,m)
      dimension eig01(n1), eig02(n2), eig03(n3)
      dimension u0l1(n1,n1), u0l2(n2,n2), u0l3(n3,n3)
      dimension u0r1(n1,n1), u0r2(n2,n2), u0r3(n3,n3)
      dimension t1(*), t2(*)
      common/io/inp, iout
c                 DEBUG INFORMATION      
c      title='u0l1'
c      call prntcm(title,u0l1,n1,n1,n1,n1,iout)
c      title='u0l2'
c      call prntcm(title,u0l2,n2,n2,n2,n2,iout)
c      title='u0r1'
c      call prntcm(title,u0r1,n1,n1,n1,n1,iout)
c      title='u0r2'
c      call prntcm(title,u0r2,n2,n2,n2,n2,iout) 
c      call cebtc(t1,u0l1,u0r1,n1,n1,n1)
c      title='unit matrix'
c      call prntcm(title,t1,n1,n1,n1,n1,iout)
c      call cebtc(t1,u0l2,u0r2,n2,n2,n2)
c      title='unit matrix'
c      call prntcm(title,t1,n2,n2,n2,n2,iout)          
      if(.not.drctv) then
         do 10 i=1,m
            do 20 j=1,n
               test=eig(i) - diag(j)
               if(abs(test).ge.nrzero) then
                  vec(j,i) = resid(j,i)/test
               else
                  vec(j,i) = one
               endif
 20         continue
 10      continue   
      else
         if(dim.eq.1) then
c
c         transform residual from dvr to h0 representation
c 
            call cdvr2h0(resid,vec,u0l1,u0l2,u0l3,t1,t2,
     1                   n1,n2,n3,n,m,dim)

c
c         using the new and diagonal representation solve the residual 
c         equation approximately
c
            do 30 i=1,m
               do 40 j=1,n1
                  test=eig(i) - eig01(j)
                  if(abs(test).ge.nrzero) then
                      vec(j,i) = vec(j,i)/test
                  else
                      vec(j,i) = one
                  endif
 40            continue
 30         continue   
c
c           back transform to dvr representation
c
            call ch02dvr(resid,vec,u0r1,u0r2,u0r3,t1,t2,
     1                   n1,n2,n3,n,m,dim)
            call cc2opy(resid,vec,n*m)
         elseif(dim.eq.2) then      
c
c        transform residual from  dvr to h0 representation
c
            call cdvr2h0(resid,vec,u0l1,u0l2,u0l3,t1,t2,
     1                   n1,n2,n3,n,m,dim)
c
c         using the new and diagonal representation solve the 
c         residual equation approximately
c
            do 50 i=1,m
               count=0
               do 60 j=1,n1
                  do 70 k=1,n2
                     count=count+1
                     test=eig(i) - eig01(j) - eig02(k)
                     if(abs(test).ge.nrzero) then
                        vec(count,i) = vec(count,i)/test
                     else
                         vec(count,i) = one
                     endif
 70               continue
 60            continue   
 50         continue
c
c          back transform to dvr representation
c
            call ch02dvr(resid,vec,u0r1,u0r2,u0r3,t1,t2,
     1                   n1,n2,n3,n,m,dim)
            call cc2opy(resid,vec,n*m)
         elseif(dim.eq.3) then
c
c        transform residual from  dvr to h0 representation
c
            call cdvr2h0(resid,vec,u0l1,u0l2,u0l3,t1,t2
     1                   n1,n2,n3,n,m,dim)
c
c         using the new and diagonal representation solve the 
c         residual equation approximately
c
            do 80 i=1,m
               count=0
               do 90 j=1,n1
                  do 100 k=1,n2
                     do 110 l=1,n3                      
                        count=count+1
                        test=eig(i) - eig01(j) - eig02(k) - eig03(l)
                        vec(count,i) = vec(count,i)/test
 110                 continue
 100              continue   
 90            continue
 80         continue   
c
c          back transform to dvr representation
c
            call ch02dvr(resid,vec,u0r1,u0r2,u0r3,t1,t2,
     1                   n1,n2,n3,n,m,dim)
            call cc2opy(resid,vec,n*m)
         endif
      endif
      if(prnt) then
         title='new trial vectors iteration = '//itoc(iter)
         call prntcm(title,vec,n,m,n,m,iout)
      endif
      return
      end       





