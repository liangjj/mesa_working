*deck @(#)fixlag.f	3.1  11/20/92
      subroutine fixlag(nbf,a,b)
      real*8 a(nbf,nbf),b(nbf,nbf)
c
      call scopy(nbf*nbf,a,1,b,1)
c
      do 1 i=1,nbf
        do 2 j=1,nbf
        a(j,i)=a(j,i)+b(i,j)
  2     continue
  1   continue
c
      do 3 i=1,nbf
        do 4 j=1,nbf
        a(j,i)=a(j,i)*.5
  4     continue
  3   continue
c
      return
      end
