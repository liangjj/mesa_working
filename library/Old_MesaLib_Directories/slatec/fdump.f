*deck fdump
      subroutine fdump
c***begin prologue  fdump
c***purpose  symbolic dump (should be locally written).
c***library   slatec (xerror)
c***category  r3
c***type      all (fdump-a)
c***keywords  error, xermsg
c***author  jones, r. e., (snla)
c***description
c
c        ***note*** machine dependent routine
c        fdump is intended to be replaced by a locally written
c        version which produces a symbolic dump.  failing this,
c        it should be replaced by a version which prints the
c        subprogram nesting list.  note that this dump must be
c        printed on each of up to five files, as indicated by the
c        xgetua routine.  see xsetua and xgetua for details.
c
c     written by ron jones, with slatec common math library subcommittee
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  fdump
c***first executable statement  fdump
      return
      end
