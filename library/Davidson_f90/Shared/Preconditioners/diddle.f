*deck diddle.f
      subroutine diddle(ham,n)
      implicit integer (a-z)
      real*8 ham, eps 
      dimension ham(n,n)
      common/io/inp, iout
      eps=.1d0
      do 10 i=1,n
         ham(i,i)= ham(i,i) + eps/i
 10   continue   
      return
      end       
