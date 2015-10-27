*deck @(#)geneig.f	5.1  11/6/94
      subroutine geneig(k,p,q,y,tol)
c.01/13/77
c  solution of the generalized eigenvalue problem
c              p * y = q * y * r  ,
c  where p and q are given symmetric k*k matrices and q is positive
c  definite, y is the resulting matrix of eigenvectors, such that
c              y' * q * y = 1
c  (y' is the transpose of y, and 1 is a k*k unit matrix), and r is
c  a diagonal matrix of eigenvalues.
c
c  the matrices p and q are stored as linear arrays (length=k*(k+1)/2)
c  by rows of the lower triangle, while y is stored as a linear
c  array (length=k**2) by columns, each column corresponding to one
c  eigenvector.  the diagonalized matrix r=y'*p*y replaces p, while
c  q is replaced by a copy of the original p.
c
c  schmidt orthogonalization is used to transform to the normal
c  eigenvalue problem, which is then solved by the jacobi method.
c  the orthogonalizing transformation matrix is used as the initial
c  matrix of eigenvectors, so that the transformed eigenvectors (in
c  the original basis) are produced automatically.  jacobi plane
c  rotations are carried out until abs(p(i,j)).le.tol*abs(p(i,i)-p(j,j))
c  for each off-diagonal element p(i,j).
c
      implicit real*8 (a-h,o-z)
c
      dimension p(*),q(*),y(*)
c
      data ismax/50/
c          (ismax is the maximum number of jacobi sweeps allowed)
c
c  schmidt orthogonalization section
c
      y(1)=1.0d0/sqrt(q(1))
      if (k.eq.1) go to 7
      iz=1
      iy=k
      do 1 i=2,k
        do 11 j=1,k
 11       y(iy+j)=0.0d0
        y(iy+i)=1.0d0
        snorm=q(iz+i)
        jy=0
        do 13 j=1,(i-1)
          do 131 l=1,j
 131        y(iy+j)=y(iy+j)+q(iz+l)*y(jy+l)
          snorm=snorm-y(iy+j)**2
 13       jy=jy+k
        jy=0
        do 14 j=1,(i-1)
          s=0.0d0
          ly=jy
          do 141 l=j,(i-1)
            s=s-y(ly+j)*y(iy+l)
 141        ly=ly+k
          y(iy+j)=s
 14       jy=jy+k
        snorm=sqrt(snorm)
        do 15 j=1,i
 15       y(iy+j)=y(iy+j)/snorm
        iy=iy+k
 1      iz=iz+i
c
c  transformation section
c
      t=y(1)
c     (y(1) is stored temporarily in t in order to make room for a
c     work-space array in the first column of y)
c
      ij=0
      do 2 i=1,k
        do 21 j=1,i
          ij=ij+1
 21       q(ij)=p(ij)
 2      continue
      y(1)=t*q(1)
      p(1)=t*y(1)
      ij=1
      iz=k
      do 3 i=2,k
        lm=0
        do 31 l=1,i
          s=0.0d0
          do 311 m=1,l
            lm=lm+1
 311        s=s+y(iz+m)*q(lm)
          if (l.eq.i) go to 31
          ml=lm+l
          lp=l+1
          do 312 m=lp,i
            s=s+y(iz+m)*q(ml)
 312        ml=ml+m
 31       y(l)=s
        ij=ij+1
        p(ij)=t*y(1)
        jz=k
        do 32 j=2,i
          ij=ij+1
          p(ij)=0.0d0
          do 321 m=1,j
 321        p(ij)=p(ij)+y(jz+m)*y(m)
 32       jz=jz+k
 3      iz=iz+k
      y(1)=t
      do 4 i=2,k
 4      y(i)=0.0d0
c
c  jacobi procedure section
c
      nrot=0
      do 5 isweep=1,ismax
        nrots=0
        iy=k
        iz=1
        do 51 i=2,k
          ip=i+1
          ii=iz+i
          jy=0
          jz=0
          do 511 j=1,(i-1)
            jj=jz+j
            a=p(iz+j)
            h=p(ii)-p(jj)
            if (abs(a).le.tol*abs(h)) go to 5111
            nrots=nrots+1
            t=1.0d0
            b=h/(a+a)
            if (b.ne.0.0d0) t=1.0d0/(b+sign(sqrt(1.0d0+b**2),b))
            c=1.0d0/sqrt(1.0d0+t**2)
            s=c*t
            tau=s/(1.0d0+c)
            d=t*a
            p(ii)=p(ii)+d
            p(jj)=p(jj)-d
            p(iz+j)=0
              do 5113 l=1,j-1
                a=p(jz+l)
                b=p(iz+l)
                p(jz+l)=a-s*(b+a*tau)
 5113           p(iz+l)=b+s*(a-b*tau)
            if (j.ne.(i-1)) then
              jp=j+1
              lz=jj
              do 5115 l=jp,(i-1)
                a=p(lz+j)
                b=p(iz+l)
                p(lz+j)=a-s*(b+a*tau)
                p(iz+l)=b+s*(a-b*tau)
 5115           lz=lz+l
            endif
            if (i.ne.k) then
              lz=ii
              do 5117 l=ip,k
                a=p(lz+j)
                b=p(lz+i)
                p(lz+j)=a-s*(b+a*tau)
                p(lz+i)=b+s*(a-b*tau)
 5117           lz=lz+l
            endif
              do 5118 l=1,k
                a=y(jy+l)
                b=y(iy+l)
                y(jy+l)=a-s*(b+a*tau)
 5118           y(iy+l)=b+s*(a-b*tau)
 5111       jz=jj
 511        jy=jy+k
          iz=ii
 51       iy=iy+k
        if (nrots.eq.0) go to 8
        nrot=nrot+nrots
        if (k.eq.2) go to 8
 5      continue
c
c  if came through here, convergence has not been achieved
c
      write (iw,6) ismax,nrot
 6      format (/' ** no convergence in jacobi diagonalization of',
     x  ' the interaction matrix **'/' number of sweeps =',i8/
     x  ' number of rotations =',i5)
c
      stop
c
 7    a=p(1)/q(1)
      q(1)=p(1)
      p(1)=a
      return
c
c  order p and y by increasing eigenvalues
c
 8    do 81 i=1,k-1
        l=k-i
        jy=0
        jj=0
        do 811 j=1,l
          jp=j+1
          jj=jj+j
          if (p(jj+jp).ge.p(jj)) go to 811
c
c         exchange roots j and j+1
c
          jz=jj-j
            do 8112 m=1,j-1
              a=p(jz+m)
              p(jz+m)=p(jj+m)
 8112         p(jj+m)=a
          jjp=jj+jp
          a=p(jj)
          p(jj)=p(jjp)
          p(jjp)=a
          if (jp.ne.k) then
            mj=jjp+j
            jpp=jp+1
            do 8114 m=jpp,k
              a=p(mj)
              p(mj)=p(mj+1)
              p(mj+1)=a
 8114         mj=mj+m
          endif
            do 8115 m=1,k
              jm=jy+m
              a=y(jm)
              y(jm)=y(jm+k)
 8115         y(jm+k)=a
 811      jy=jy+k
 81     continue
      return
      end
