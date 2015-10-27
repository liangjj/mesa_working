*deck @(#)icopy.f	1.1 9/6/91 
      subroutine icopy(ivec,jvec,n)
      implicit integer (a-z)
      integer ivec(n), jvec(n)
c
c
      do 100 i=1,n
         jvec(i)=ivec(i)
  100 continue
c
c
      return
      end
