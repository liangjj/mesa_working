*deck dfulmt
      subroutine dfulmt (i, j, aij, indcat, prgopt, dattrv, iflag)
c***begin prologue  dfulmt
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (fulmat-s, dfulmt-d)
c***author  (unknown)
c***description
c
c     decodes a standard two-dimensional fortran array passed
c     in the array dattrv(ia,*).  the row dimension ia and the
c     matrix dimensions mrelas and nvars must simultaneously be
c     passed using the option array, prgopt(*).  it is an error
c     if this data is not passed to dfulmt( ).
c     example-- (for use together with dsplp().)
c      external dusrmt
c      dimension dattrv(ia,*)
c      prgopt(01)=7
c      prgopt(02)=68
c      prgopt(03)=1
c      prgopt(04)=ia
c      prgopt(05)=mrelas
c      prgopt(06)=nvars
c      prgopt(07)=1
c     call dsplp(  ... dfulmt instead of dusrmt...)
c
c***see also  dsplp
c***routines called  xermsg
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c***end prologue  dfulmt
      double precision aij,zero,dattrv(*),prgopt(*)
      integer iflag(10)
      save zero
c***first executable statement  dfulmt
      if (.not.(iflag(1).eq.1)) go to 50
c     initialize pointers to process full two-dimensional fortran
c     arrays.
      zero = 0.d0
      lp = 1
   10 next = prgopt(lp)
      if (.not.(next.le.1)) go to 20
      nerr = 29
      level = 1
      call xermsg ('slatec', 'dfulmt',
     +   'in dsplp, row dim., mrelas, nvars are missing from prgopt.',
     +   nerr, level)
      iflag(1) = 3
      go to 110
   20 key = prgopt(lp+1)
      if (.not.(key.ne.68)) go to 30
      lp = next
      go to 10
   30 if (.not.(prgopt(lp+2).eq.zero)) go to 40
      lp = next
      go to 10
   40 iflag(2) = 1
      iflag(3) = 1
      iflag(4) = prgopt(lp+3)
      iflag(5) = prgopt(lp+4)
      iflag(6) = prgopt(lp+5)
      go to 110
   50 if (.not.(iflag(1).eq.2)) go to 100
   60 i = iflag(2)
      j = iflag(3)
      if (.not.(j.gt.iflag(6))) go to 70
      iflag(1) = 3
      go to 110
   70 if (.not.(i.gt.iflag(5))) go to 80
      iflag(2) = 1
      iflag(3) = j + 1
      go to 60
   80 aij = dattrv(iflag(4)*(j-1)+i)
      iflag(2) = i + 1
      if (.not.(aij.eq.zero)) go to 90
      go to 60
   90 indcat = 0
      go to 110
  100 continue
  110 return
      end
