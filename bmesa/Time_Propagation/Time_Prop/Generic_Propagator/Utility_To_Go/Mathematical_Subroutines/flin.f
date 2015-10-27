*deck @(#)flin.f	5.1  11/6/94
      subroutine flin(a,idim,in,im,det)
c***begin prologue     flin
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           linear, simultaneous, equations
c***author             saxe, paul (lanl)
c***source             @(#)flin.f	5.1   11/6/94
c***purpose            solves b*x=c for x.
c***description
c                      call flin(a,idim,in,im,det)
c                      solves the linear equation b(in,in)*x(in,im)=c(in,im)
c                      for x.
c                      the input storage is a little strange. b and c are
c                      packed one after the other in the array a(in,(in+im)).  c
c                      occupied by c.  i believe idim=in.  the determinant
c                      of b is returned in det.
c
c***references
c***routines called    sdot(math), abs
c***end prologue       flin
      implicit real*8 (a-h,o-z)
c
c
c     linear simultaneous equation
c
c     a(in*in) * x(in*im) = b(in*im)
c
c     a & b should be stored on a(in*(in+im))
c     solution x will be stored on b part in dimension a.
c
      dimension a(idim,1)
      real*8 zero,one
      parameter (zero=0.0d+00,one=1.0d+00)
c
      n=in
      nr=im
      jmax=n+nr
      sign=one
c m is the stage of elimination
      do 49 m=1,n
         temp=zero
         do 41 i=m,n
            if(m.gt.1)a(i,m)=a(i,m)-sdot(m-1,a(i,1),idim,a(1,m),1)
            if(abs(a(i,m)).le.temp)goto 41
            temp=abs(a(i,m))
            imax=i
 41      continue
         if(temp.le.0.0)goto 999
         if(imax.eq.m)goto 45
         sign=-sign
         do 44 j=1,jmax
            stor=a(m,j)
            a(m,j)=a(imax,j)
            a(imax,j)=stor
 44      continue
 45      continue
         jj=m+1
         if(jj.gt.jmax)goto 49
         if(m.gt.1)goto 47
         do 46 j=jj,jmax
            a(1,j)=a(1,j)/a(1,1)
 46      continue
         d=a(1,1)
         goto 49
 47      continue
         do 48 j=jj,jmax
            a(m,j)=(a(m,j)-sdot(m-1,a(m,1),idim,a(1,j),1))/a(m,m)
 48      continue
         d=d*a(m,m)
 49   continue
      if(nr.eq.0)return
      do 59 i=1,nr
         npi=n+i
         do 58 k=2,n
            j=n+1-k
            a(j,npi)=a(j,npi)-sdot(k-1,a(j,j+1),idim,a(j+1,npi),1)
 58      continue
 59   continue
      det=d*sign
      return
c on zero pivot, set det=0.and return to calling program nov 1972
 999  det=zero
c
c
      return
      end
