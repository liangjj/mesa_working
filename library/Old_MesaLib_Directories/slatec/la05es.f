*deck la05es
      subroutine la05es (a, irn, ip, n, iw, ia, reals)
c***begin prologue  la05es
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (la05es-s, la05ed-d)
c***author  (unknown)
c***description
c
c     this subprogram is a slight modification of a subprogram
c     from the c. 1979 aere harwell library.  the name of the
c     corresponding harwell code can be obtained by deleting
c     the final letter =s= in the names used here.
c     revised sep. 13, 1979.
c
c     royalties have been paid to aere-uk for use of their codes
c     in the package given here.  any primary usage of the harwell
c     subroutines requires a royalty agreement and payment between
c     the user and aere-uk.  any usage of the sandia written codes
c     splp( ) (which uses the harwell subroutines) is permitted.
c
c***see also  splp
c***routines called  (none)
c***common blocks    la05ds
c***revision history  (yymmdd)
c   811215  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  la05es
      logical reals
      real a(*)
      integer irn(*), iw(*)
      integer ip(*)
      common /la05ds/ small, lp, lenl, lenu, ncp, lrow, lcol
c***first executable statement  la05es
      ncp = ncp + 1
c     compress file of positive integers. entry j starts at irn(ip(j))
c  and contains iw(j) integers,j=1,n. other components of irn are zero.
c  length of compressed file placed in lrow if reals is .true. or lcol
c  otherwise.
c  if reals is .true. array a contains a real file associated with irn
c  and this is compressed too.
c  a,irn,ip,iw,ia are input/output variables.
c  n,reals are input/unchanged variables.
c
      do 10 j=1,n
c store the last element of entry j in iw(j) then overwrite it by -j.
         nz = iw(j)
         if (nz.le.0) go to 10
         k = ip(j) + nz - 1
         iw(j) = irn(k)
         irn(k) = -j
   10 continue
c kn is the position of next entry in compressed file.
      kn = 0
      ipi = 0
      kl = lcol
      if (reals) kl = lrow
c loop through the old file skipping zero (dummy) elements and
c     moving genuine elements forward. the entry number becomes
c     known only when its end is detected by the presence of a negative
c     integer.
      do 30 k=1,kl
         if (irn(k).eq.0) go to 30
         kn = kn + 1
         if (reals) a(kn) = a(k)
         if (irn(k).ge.0) go to 20
c end of entry. restore irn(k), set pointer to start of entry and
c     store current kn in ipi ready for use when next last entry
c     is detected.
         j = -irn(k)
         irn(k) = iw(j)
         ip(j) = ipi + 1
         iw(j) = kn - ipi
         ipi = kn
   20    irn(kn) = irn(k)
   30 continue
      if (reals) lrow = kn
      if (.not.reals) lcol = kn
      return
      end
