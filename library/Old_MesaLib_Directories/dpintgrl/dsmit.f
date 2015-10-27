*deck  @(#)dsmit.f	5.1 11/6/94
      subroutine dsmit(c,s,n)
c
c     returns an n*n triangular matrix c such that c(transpose)sc=i,
c     where i is a n*n identity matrix
c
      implicit double precision (a-h,o-z)
c
      dimension c(10,10),s(10,10),v(10),y(10)
c
      data a0, a1 /0.0d+00, 1.0d+00/
      save a0,a1
c
      common/io/inp,iout
c
      do 10 i=1,n
         do 10 j=1,i
            c(i,j)=a0
   10 continue
      do 100 j=1,n
         kmax=j-1
         fac=s(j,j)
         if(kmax.ne.0) then
            do 20 k=1,kmax
               v(k)=a0
               y(k)=s(k,j)
   20       continue
            do 50 k=1,kmax
               dot=a0
               do 30 i=1,k
                  dot=c(i,k)*y(i)+dot
   30          continue
               do 40 i=1,k
                  v(i)=v(i)-dot*c(i,k)
   40          continue
               fac=fac-dot*dot
   50       continue
         endif
         fac=a1/sqrt(fac)
         c(j,j)=fac
         if(kmax.ne.0) then
            do 70 k=1,kmax
               c(k,j)=fac*v(k)
   70       continue
         endif
  100 continue
c
c
      return
      end
