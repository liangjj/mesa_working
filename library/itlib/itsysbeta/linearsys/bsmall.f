*deck bsmall.f
c***begin prologue     bsmall
c***date written       980420   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           small davidson matrix
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       bsmall
      subroutine bsmall(h,htmp,b,btmp,vec,hvec,rhs,energy,nrhs,
     1                  n,begin,end,maxvec,action)
      implicit integer (a-z)
      real*8 h, htmp, b, btmp, vec, hvec, rhs, energy, ddot
      character*(*) action
      dimension h(maxvec,*), htmp(maxvec,*)
      dimension b(maxvec,nrhs), btmp(maxvec,nrhs)
      dimension  vec(n,*), hvec(n,*), rhs(n,nrhs)
      common/io/inp, iout
      if(action.eq.'initialize') then
         do 10 i=1,end
            do 20 j=1,i
               h(i,j) = - ddot(n,vec(1,i),1,hvec(1,j),1)
               h(j,i) = h(i,j)
               htmp(i,j) = h(i,j)
               htmp(j,i) = h(j,i)
 20         continue
            h(i,i) = energy + h(i,i)
            htmp(i,i) = h(i,i)
            do 30 j=1,nrhs
               b(i,j) = ddot(n,vec(1,i),1,rhs(1,j),1)
               btmp(i,j) = b(i,j)
 30         continue   
 10      continue   
      elseif(action.eq.'fill') then
         do 40 i=1,end
            do 50 j=begin,end
               h(i,j) = - ddot(n,vec(1,i),1,hvec(1,j),1)
               h(j,i) = h(i,j)
 50         continue
            do 60 j=1,nrhs
                b(i,j) = ddot(n,vec(1,i),1,rhs(1,j),1)   
 60         continue    
 40      continue
         do 70 i=begin,end
            h(i,i) = energy + h(i,i)
 70      continue   
         do 80 i=1,end   
            do 90 j=1,i
               htmp(i,j) = h(i,j)
               htmp(j,i) = h(i,j)
 90         continue
            do 100 j=1,nrhs
                btmp(i,j) = b(i,j)
 100        continue   
 80      continue  
      endif
      return
      end

