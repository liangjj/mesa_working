      subroutine cbabk3(nm,n,low,igh,scale,m,zr,zi,yr,yi,work,task)
c***begin prologue  cbabk3
c***date written   760101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  d4c4
c***keywords  eigenvalues,eigenvectors,eispack
c***author  smith, b. t., et al.
c***purpose  forms eigenvectors of complex general matrix from
c            eigenvectors of matrix output from cbal.
c***description
c
c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  cbal.
c
c        scale contains information determining the permutations
c          and scaling factors used by  cbal.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema, c. b. moler, *matrix eigen-
c                 system routines - eispack guide*, springer-verlag,
c                 1976.
c***routines called  (none)
c***end prologue  cbabk2
c
      implicit real *8 (a-h,o-z)
      integer i,j,k,m,n,ii,nm,igh,low
      real*8 scale(n), zr(nm,m), zi(nm,m), yr(nm,m),yi(nm,m),work(n)
      real*8 s
      character*(*) task
c
c***first executable statement  cbabk2
      if (m.eq.0) then
          return
       endif   
      if(task.eq.'right eigenvectors') then
         if (igh .eq. low) go to 120
c
         do 110 i = low, igh
            do 100 j = 1, m
               zr(i,j) = zr(i,j) * scale(i)
               zi(i,j) = zi(i,j) * scale(i)
  100       continue
c
  110    continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120    do 140 ii = 1, n
            i = ii
            if (i .ge. low .and. i .le. igh) go to 140
            if (i .lt. low) i = low - ii
            k = scale(i)
            if (k .eq. i) go to 140
c
            do 130 j = 1, m
               s = zr(i,j)
               zr(i,j) = zr(k,j)
               zr(k,j) = s
               s = zi(i,j)
               zi(i,j) = zi(k,j)
               zi(k,j) = s
  130       continue
c
  140    continue
c
      elseif(task.eq.'left eigenvectors') then
         call vinv(work,scale,n) 
         if (igh .eq. low) go to 220
c
         do 210 i = low, igh
            do 200 j = 1, m
               yr(i,j) = yr(i,j) * work(i)
               yi(i,j) = yi(i,j) * work(i)
  200       continue
c
  210    continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  220    do 240 ii = 1, n
            i = ii
            if (i .ge. low .and. i .le. igh) go to 240
            if (i .lt. low) i = low - ii
            k = scale(i)
            if (k .eq. i) go to 240
c
            do 230 j = 1, m
               s = yr(i,j)
               yr(i,j) = yr(k,j)
               yr(k,j) = s
               s = yi(i,j)
               yi(i,j) = yi(k,j)
               yi(k,j) = s
  230       continue
c
  240    continue
c
      elseif(task.eq.'left and right eigenvectors') then
         call vinv(work,scale,n)
         if (igh .eq. low) go to 320
c
         do 310 i = low, igh
            do 300 j = 1, m
               yr(i,j) = zr(i,j) * work(i)
               yi(i,j) = zi(i,j) * work(i)
               zr(i,j) = zr(i,j) * scale(i)
               zi(i,j) = zi(i,j) * scale(i)
  300       continue
c
  310    continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  320    do 340 ii = 1, n
            i = ii
            if (i .ge. low .and. i .le. igh) go to 340
            if (i .lt. low) i = low - ii
            k = scale(i)
            if (k .eq. i) go to 340
c
            do 330 j = 1, m
               s = zr(i,j)
               zr(i,j) = zr(k,j)
               zr(k,j) = s
               s = zi(i,j)
               zi(i,j) = zi(k,j)
               zi(k,j) = s
               s = yr(i,j)
               yr(i,j) = yr(k,j)
               yr(k,j) = s
               s = yi(i,j)
               yi(i,j) = yi(k,j)
               yi(k,j) = s
  330       continue
c
  340    continue
      endif  
c                      
      return
      end
