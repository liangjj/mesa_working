*deck cscale
      subroutine cscale (a, nrda, nrow, ncol, cols, colsav, rows,
     +   rowsav, anorm, scales, iscale, ic)
c***begin prologue  cscale
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (cscale-s, dcscal-d)
c***author  watts, h. a., (snla)
c***description
c
c     this routine scales the matrix a by columns when needed
c
c***see also  bvsup
c***routines called  sdot
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  cscale
      dimension a(nrda,*),cols(*),colsav(*),scales(*),
     1          rows(*),rowsav(*)
c
      save ten4, ten20
      data ten4,ten20/1.e+4,1.e+20/
c
c***first executable statement  cscale
      if (iscale .ne. (-1)) go to 25
c
      if (ic .eq. 0) go to 10
      do 5 k=1,ncol
    5    cols(k)=sdot(nrow,a(1,k),1,a(1,k),1)
c
   10 ascale=anorm/ncol
      do 20 k=1,ncol
         cs=cols(k)
         if ((cs .gt. ten4*ascale) .or. (ten4*cs .lt. ascale)) go to 50
         if ((cs .lt. 1./ten20) .or. (cs .gt. ten20)) go to 50
   20 continue
c
   25 do 30 k=1,ncol
   30    scales(k)=1.
      return
c
   50 alog2=log(2.)
      anorm=0.
      do 100 k=1,ncol
         cs=cols(k)
         if (cs .ne. 0.) go to 60
         scales(k)=1.
         go to 100
   60    p=log(cs)/alog2
         ip=-0.5*p
         s=2.**ip
         scales(k)=s
         if (ic .eq. 1) go to 70
         cols(k)=s*s*cols(k)
         anorm=anorm+cols(k)
         colsav(k)=cols(k)
   70    do 80 j=1,nrow
   80       a(j,k)=s*a(j,k)
  100 continue
c
      if (ic .eq. 0) return
c
      do 200 k=1,nrow
         rows(k)=sdot(ncol,a(k,1),nrda,a(k,1),nrda)
         rowsav(k)=rows(k)
  200    anorm=anorm+rows(k)
      return
      end
