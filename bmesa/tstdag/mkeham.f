*deck mkeham.f
      subroutine mkeham(ham,n)
      implicit integer (a-z)
      real*8 ham
      character*80 title
      dimension ham(n,n)
      common/io/inp, iout 
      do 10 i=1,n
         do 20 j=1,n
            ham(i,j)=-.2d0
 20      continue   
 10   continue   
      do 30 i=1,n
         ham(i,i)=-(n-i+1)
 30   continue   
      return
      end       

