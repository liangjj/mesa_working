*deck hsmall.f
c***begin prologue     hsmall
c***date written       980420   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           small davidson matrix
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       hsmall
      subroutine hsmall(b,btmp,vec,hvec,n,begin,end,maxvec,action)
      implicit integer (a-z)
      real*8 b, btmp, vec, hvec, sdot
      character*(*) action
      dimension b(maxvec,*), btmp(maxvec,*), vec(n,*), hvec(n,*)
      common/io/inp, iout
      if(action.eq.'initialize') then
         do 10 i=1,end
            do 20 j=1,i
               b(i,j) = sdot(n,vec(1,i),1,hvec(1,j),1)
               b(j,i) = b(i,j)
               btmp(i,j) = b(i,j)
               btmp(j,i) = b(j,i)
 20         continue
 10      continue   
      elseif(action.eq.'fill') then
         do 30 i=1,end
            do 40 j=begin,end
               b(i,j) = sdot(n,vec(1,i),1,hvec(1,j),1)
               b(j,i) = b(i,j)
 40         continue
 30      continue
         do 50 i=1,end   
            do 60 j=1,i
               btmp(i,j) = b(i,j)
               btmp(j,i) = b(i,j)
 60         continue
 50      continue     
      endif
      return
      end

