*deck dcscal
      subroutine dcscal (a, nrda, nrow, ncol, cols, colsav, rows,
     +   rowsav, anorm, scales, iscale, ic)
c***begin prologue  dcscal
c***subsidiary
c***purpose  subsidiary to dbvsup and dsuds
c***library   slatec
c***type      double precision (cscale-s, dcscal-d)
c***author  watts, h. a., (snla)
c***description
c
c     this routine scales the matrix a by columns when needed.
c
c***see also  dbvsup, dsuds
c***routines called  ddot
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dcscal
      double precision ddot
      integer ic, ip, iscale, j, k, ncol, nrda, nrow
      double precision a(nrda,*), alog2, anorm, ascale, cols(*),
     1     colsav(*), cs, p, rows(*), rowsav(*), s,
     2     scales(*), ten20, ten4
c
      save ten4, ten20
      data ten4,ten20 /1.0d4,1.0d20/
c
c     begin block permitting ...exits to 130
c        begin block permitting ...exits to 60
c***first executable statement  dcscal
            if (iscale .ne. (-1)) go to 40
c
               if (ic .eq. 0) go to 20
                  do 10 k = 1, ncol
                     cols(k) = ddot(nrow,a(1,k),1,a(1,k),1)
   10             continue
   20          continue
c
               ascale = anorm/ncol
               do 30 k = 1, ncol
                  cs = cols(k)
c        .........exit
                  if ((cs .gt. ten4*ascale) .or. (ten4*cs .lt. ascale))
     1               go to 60
c        .........exit
                  if ((cs .lt. 1.0d0/ten20) .or. (cs .gt. ten20))
     1               go to 60
   30          continue
   40       continue
c
            do 50 k = 1, ncol
               scales(k) = 1.0d0
   50       continue
c     ......exit
            go to 130
   60    continue
c
         alog2 = log(2.0d0)
         anorm = 0.0d0
         do 110 k = 1, ncol
            cs = cols(k)
            if (cs .ne. 0.0d0) go to 70
               scales(k) = 1.0d0
            go to 100
   70       continue
               p = log(cs)/alog2
               ip = -0.5d0*p
               s = 2.0d0**ip
               scales(k) = s
               if (ic .eq. 1) go to 80
                  cols(k) = s*s*cols(k)
                  anorm = anorm + cols(k)
                  colsav(k) = cols(k)
   80          continue
               do 90 j = 1, nrow
                  a(j,k) = s*a(j,k)
   90          continue
  100       continue
  110    continue
c
c     ...exit
         if (ic .eq. 0) go to 130
c
         do 120 k = 1, nrow
            rows(k) = ddot(ncol,a(k,1),nrda,a(k,1),nrda)
            rowsav(k) = rows(k)
            anorm = anorm + rows(k)
  120    continue
  130 continue
      return
      end
