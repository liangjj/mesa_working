*deck comlr2
      subroutine comlr2 (nm, n, low, igh, int, hr, hi, wr, wi, zr, zi,
     +   ierr)
c***begin prologue  comlr2
c***purpose  compute the eigenvalues and eigenvectors of a complex upper
c            hessenberg matrix using the modified lr method.
c***library   slatec (eispack)
c***category  d4c2b
c***type      complex (comlr2-c)
c***keywords  eigenvalues, eigenvectors, eispack, lr method
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure comlr2,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the modified lr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  comhes  has been used to reduce
c     this general matrix to hessenberg form.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, hr, hi, zr and zi, as declared in the
c          calling program dimension statement.  nm is an integer
c          variable.
c
c        n is the order of the matrix h=(hr,hi).  n is an integer
c          variable.  n must be less than or equal to nm.
c
c        low and igh are two integer variables determined by the
c          balancing subroutine  cbal.  if  cbal  has not been used,
c          set low=1 and igh equal to the order of the matrix, n.
c
c        int contains information on the rows and columns
c          interchanged in the reduction by  comhes, if performed.
c          only elements low through igh are used.  if you want the
c          eigenvectors of a complex general matrix, leave int as it
c          came from  comhes.  if the eigenvectors of the hessenberg
c          matrix are desired, set int(j)=j for these elements.  int
c          is a one-dimensional integer array, dimensioned int(igh).
c
c        hr and hi contain the real and imaginary parts, respectively,
c          of the complex upper hessenberg matrix.  their lower
c          triangles below the subdiagonal contain the multipliers
c          which were used in the reduction by  comhes, if performed.
c          if the eigenvectors of a complex general matrix are
c          desired, leave these multipliers in the lower triangles.
c          if the eigenvectors of the hessenberg matrix are desired,
c          these elements must be set to zero.  hr and hi are
c          two-dimensional real arrays, dimensioned hr(nm,n) and
c          hi(nm,n).
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed, but the location hr(1,1) contains the norm
c          of the triangularized matrix.
c
c        wr and wi contain the real and imaginary parts, respectively,
c          of the eigenvalues of the upper hessenberg matrix.  if an
c          error exit is made, the eigenvalues should be correct for
c          indices ierr+1, ierr+2, ..., n.  wr and wi are one-
c          dimensional real arrays, dimensioned wr(n) and wi(n).
c
c        zr and zi contain the real and imaginary parts, respectively,
c          of the eigenvectors.  the eigenvectors are unnormalized.
c          if an error exit is made, none of the eigenvectors has been
c          found.  zr and zi are two-dimensional real arrays,
c          dimensioned zr(nm,n) and zi(nm,n).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after a total of 30*n iterations.
c                     the eigenvalues should be correct for indices
c                     ierr+1, ierr+2, ..., n, but no eigenvectors are
c                     computed.
c
c     calls csroot for complex square root.
c     calls cdiv for complex division.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  cdiv, csroot
c***revision history  (yymmdd)
c   760101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  comlr2
c
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,nm,nn,igh,im1,ip1
      integer itn,its,low,mp1,enm1,iend,ierr
      real hr(nm,*),hi(nm,*),wr(*),wi(*),zr(nm,*),zi(nm,*)
      real si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,s1,s2
      integer int(*)
c
c***first executable statement  comlr2
      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 100 i = 1, n
c
         do 100 j = 1, n
            zr(i,j) = 0.0e0
            zi(i,j) = 0.0e0
            if (i .eq. j) zr(i,j) = 1.0e0
  100 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by comhes ..........
      iend = igh - low - 1
      if (iend .le. 0) go to 180
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
      do 160 ii = 1, iend
         i = igh - ii
         ip1 = i + 1
c
         do 120 k = ip1, igh
            zr(k,i) = hr(k,i-1)
            zi(k,i) = hi(k,i-1)
  120    continue
c
         j = int(i)
         if (i .eq. j) go to 160
c
         do 140 k = i, igh
            zr(i,k) = zr(j,k)
            zi(i,k) = zi(j,k)
            zr(j,k) = 0.0e0
            zi(j,k) = 0.0e0
  140    continue
c
         zr(j,i) = 1.0e0
  160 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0e0
      ti = 0.0e0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         s1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))
     1             + abs(hr(l,l)) + abs(hi(l,l))
         s2 = s1 + abs(hr(l,l-1)) + abs(hi(l,l-1))
         if (s2 .eq. s1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
      xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0e0 .and. xi .eq. 0.0e0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0e0
      yi = (hi(enm1,enm1) - si) / 2.0e0
      call csroot(yr**2-yi**2+xr,2.0e0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0e0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = abs(hi(en,enm1)) + abs(hi(enm1,en-2))
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements ..........
      xr = abs(hr(enm1,enm1)) + abs(hi(enm1,enm1))
      yr = abs(hr(en,enm1)) + abs(hi(en,enm1))
      zzr = abs(hr(en,en)) + abs(hi(en,en))
c     .......... for m=en-1 step -1 until l do -- ..........
      do 380 mm = l, enm1
         m = enm1 + l - mm
         if (m .eq. l) go to 420
         yi = yr
         yr = abs(hr(m,m-1)) + abs(hi(m,m-1))
         xi = zzr
         zzr = xr
         xr = abs(hr(m-1,m-1)) + abs(hi(m-1,m-1))
         s1 = zzr / yi * (zzr + xr + xi)
         s2 = s1 + yr
         if (s2 .eq. s1) go to 420
  380 continue
c     .......... triangular decomposition h=l*r ..........
  420 mp1 = m + 1
c
      do 520 i = mp1, en
         im1 = i - 1
         xr = hr(im1,im1)
         xi = hi(im1,im1)
         yr = hr(i,im1)
         yi = hi(i,im1)
         if (abs(xr) + abs(xi) .ge. abs(yr) + abs(yi)) go to 460
c     .......... interchange rows of hr and hi ..........
         do 440 j = im1, n
            zzr = hr(im1,j)
            hr(im1,j) = hr(i,j)
            hr(i,j) = zzr
            zzi = hi(im1,j)
            hi(im1,j) = hi(i,j)
            hi(i,j) = zzi
  440    continue
c
         call cdiv(xr,xi,yr,yi,zzr,zzi)
         wr(i) = 1.0e0
         go to 480
  460    call cdiv(yr,yi,xr,xi,zzr,zzi)
         wr(i) = -1.0e0
  480    hr(i,im1) = zzr
         hi(i,im1) = zzi
c
         do 500 j = i, n
            hr(i,j) = hr(i,j) - zzr * hr(im1,j) + zzi * hi(im1,j)
            hi(i,j) = hi(i,j) - zzr * hi(im1,j) - zzi * hr(im1,j)
  500    continue
c
  520 continue
c     .......... composition r*l=h ..........
      do 640 j = mp1, en
         xr = hr(j,j-1)
         xi = hi(j,j-1)
         hr(j,j-1) = 0.0e0
         hi(j,j-1) = 0.0e0
c     .......... interchange columns of hr, hi, zr, and zi,
c                if necessary ..........
         if (wr(j) .le. 0.0e0) go to 580
c
         do 540 i = 1, j
            zzr = hr(i,j-1)
            hr(i,j-1) = hr(i,j)
            hr(i,j) = zzr
            zzi = hi(i,j-1)
            hi(i,j-1) = hi(i,j)
            hi(i,j) = zzi
  540    continue
c
         do 560 i = low, igh
            zzr = zr(i,j-1)
            zr(i,j-1) = zr(i,j)
            zr(i,j) = zzr
            zzi = zi(i,j-1)
            zi(i,j-1) = zi(i,j)
            zi(i,j) = zzi
  560    continue
c
  580    do 600 i = 1, j
            hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
            hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
  600    continue
c     .......... accumulate transformations ..........
         do 620 i = low, igh
            zr(i,j-1) = zr(i,j-1) + xr * zr(i,j) - xi * zi(i,j)
            zi(i,j-1) = zi(i,j-1) + xr * zi(i,j) + xi * zr(i,j)
  620    continue
c
  640 continue
c
      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.0e0
c
      do 720 i = 1, n
c
         do 720 j = i, n
            norm = norm + abs(hr(i,j)) + abs(hi(i,j))
  720 continue
c
      hr(1,1) = norm
      if (n .eq. 1 .or. norm .eq. 0.0e0) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = hr(i,en)
            zzi = hi(i,en)
            if (i .eq. enm1) go to 760
            ip1 = i + 1
c
            do 740 j = ip1, enm1
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
c
  760       yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0e0 .or. yi .ne. 0.0e0) go to 775
            yr = norm
  770       yr = 0.5e0*yr
            if (norm + yr .gt. norm) go to 770
            yr = 2.0e0*yr
  775       call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
  780    continue
c
  800 continue
c     .......... end backsubstitution ..........
      enm1 = n - 1
c     .......... vectors of isolated roots ..........
      do 840 i = 1, enm1
         if (i .ge. low .and. i .le. igh) go to 840
         ip1 = i + 1
c
         do 820 j = ip1, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low+1 do -- ..........
      do 880 jj = low, enm1
         j = n + low - jj
         m = min(j-1,igh)
c
         do 880 i = low, igh
            zzr = zr(i,j)
            zzi = zi(i,j)
c
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
c
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
