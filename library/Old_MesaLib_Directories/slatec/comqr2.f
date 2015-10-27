*deck comqr2
      subroutine comqr2 (nm, n, low, igh, ortr, orti, hr, hi, wr, wi,
     +   zr, zi, ierr)
c***begin prologue  comqr2
c***purpose  compute the eigenvalues and eigenvectors of a complex upper
c            hessenberg matrix.
c***library   slatec (eispack)
c***category  d4c2b
c***type      complex (hqr2-s, comqr2-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the qr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  corth  has been used to reduce
c     this general matrix to hessenberg form.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, hr, hi, zr, and zi, as declared in the
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
c        ortr and orti contain information about the unitary trans-
c          formations used in the reduction by  corth, if performed.
c          only elements low through igh are used.  if the eigenvectors
c          of the hessenberg matrix are desired, set ortr(j) and
c          orti(j) to 0.0e0 for these elements.  ortr and orti are
c          one-dimensional real arrays, dimensioned ortr(igh) and
c          orti(igh).
c
c        hr and hi contain the real and imaginary parts, respectively,
c          of the complex upper hessenberg matrix.  their lower
c          triangles below the subdiagonal contain information about
c          the unitary transformations used in the reduction by  corth,
c          if performed.  if the eigenvectors of the hessenberg matrix
c          are desired, these elements may be arbitrary.  hr and hi
c          are two-dimensional real arrays, dimensioned hr(nm,n) and
c          hi(nm,n).
c
c     on output
c
c        ortr, orti, and the upper hessenberg portions of hr and hi
c          have been destroyed.
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
c     calls pythag(a,b) for sqrt(a**2 + b**2).
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
c***routines called  cdiv, csroot, pythag
c***revision history  (yymmdd)
c   760101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  comqr2
c
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1
      integer itn,its,low,lp1,enm1,iend,ierr
      real hr(nm,*),hi(nm,*),wr(*),wi(*),zr(nm,*),zi(nm,*)
      real ortr(*),orti(*)
      real si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,s1,s2
      real pythag
c
c***first executable statement  comqr2
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
c                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0e0 .and. orti(i) .eq. 0.0e0) go to 140
         if (hr(i,i-1) .eq. 0.0e0 .and. hi(i,i-1) .eq. 0.0e0) go to 140
c     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
c
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
c
         do 130 j = i, igh
            sr = 0.0e0
            si = 0.0e0
c
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
c
            sr = sr / norm
            si = si / norm
c
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
c
  130    continue
c
  140 continue
c     .......... create real subdiagonal elements ..........
  150 l = low + 1
c
      do 170 i = l, igh
         ll = min(i+1,igh)
         if (hi(i,i-1) .eq. 0.0e0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0e0
c
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
c
  170 continue
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
     1             + abs(hr(l,l)) +abs(hi(l,l))
         s2 = s1 + abs(hr(l,l-1))
         if (s2 .eq. s1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
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
      si = 0.0e0
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
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0e0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0e0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0e0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0e0
      if (en .eq. n) go to 540
      ip1 = en + 1
c
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0e0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
c
  600 continue
c
      if (si .eq. 0.0e0) go to 240
c
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
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
      do  840 i = 1, enm1
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
