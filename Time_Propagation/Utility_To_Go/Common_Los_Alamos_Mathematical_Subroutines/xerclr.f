      subroutine xerclr
c***begin prologue  xerclr
c***date written   790801   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  resets current error number to zero
c***description
c
c     abstract
c        this routine simply resets the current error number to zero.
c        this may be necessary to do in order to determine that
c        a certain error has occurred again since the last time
c        numxer was referenced.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  7 june 1978
c***references  jones r.e., *slatec common mathematical library error
c                 handling package*, sand78-1189, sandia laboratories,
c                 1978.
c***routines called  j4save
c***end prologue  xerclr
c
c***first executable statement  xerclr
      junk = j4save(1,0,.true.)
      return
      end
