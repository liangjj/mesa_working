*deck seidel.f
c***begin prologue     seidel
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            perform a set of gauss-seidel relaxations.
c***references         
c
c***routines called    
c***end prologue       seidel
      subroutine seidel(resid,d,ham,root,vec,n,iter)
      implicit integer (a-z)
      real*8 resid, d, ham, root, vec, nrzero
      real*8 temp, test
      dimension resid(n), d(n), ham(n,n), vec(n)
      common/io/inp, iout 
      data nrzero / 1.0d-06 /
      do 10 i=1,n
         test = root-d(i)
         if(abs(test).ge.nrzero) then
            vec(i)=resid(i)/test
         else
            vec(i)=1.d0
         endif
 10   continue
      do 20 niter=1,iter 
         do 30 i=1,n
            temp=resid(i)
            do 40 j=1,n
               temp = temp +ham(i,j)*vec(j)
 40         continue
            test= root-d(i)
            if(abs(test).ge.nrzero) then
               vec(i)=temp/test
            else
               vec(i)=1.d0
            endif
 30         continue
 20   continue
      call copy(vec,resid,n)
      return
      end       
