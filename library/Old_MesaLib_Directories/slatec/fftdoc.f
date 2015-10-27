*deck fftdoc
      subroutine fftdoc
c***begin prologue  fftdoc
c***purpose  documentation for fftpack, a collection of fast fourier
c            transform routines.
c***library   slatec
c***category  j1, z
c***type      all (fftdoc-a)
c***keywords  documentation, fast fourier transform, fft
c***author  swarztrauber, p. n., (ncar)
c***description
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                       version 3  june 1979
c
c          a package of fortran subprograms for the fast fourier
c           transform of periodic and other symmetric sequences
c                              by
c                       paul n swarztrauber
c
c    national center for atmospheric research, boulder, colorado 80307
c        which is sponsored by the national science foundation
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     this package consists of programs which perform fast fourier
c     transforms for both complex and real periodic sequences and
c     certain other symmetric sequences that are listed below.
c
c     1.   rffti     initialize rfftf and rfftb
c     2.   rfftf     forward transform of a real periodic sequence
c     3.   rfftb     backward transform of a real coefficient array
c
c     4.   ezffti    initialize ezfftf and ezfftb
c     5.   ezfftf    a simplified real periodic forward transform
c     6.   ezfftb    a simplified real periodic backward transform
c
c     7.   sinti     initialize sint
c     8.   sint      sine transform of a real odd sequence
c
c     9.   costi     initialize cost
c     10.  cost      cosine transform of a real even sequence
c
c     11.  sinqi     initialize sinqf and sinqb
c     12.  sinqf     forward sine transform with odd wave numbers
c     13.  sinqb     unnormalized inverse of sinqf
c
c     14.  cosqi     initialize cosqf and cosqb
c     15.  cosqf     forward cosine transform with odd wave numbers
c     16.  cosqb     unnormalized inverse of cosqf
c
c     17.  cffti     initialize cfftf and cfftb
c     18.  cfftf     forward transform of a complex periodic sequence
c     19.  cfftb     unnormalized inverse of cfftf
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   780201  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900723  purpose section revised.  (wrb)
c***end prologue  fftdoc
c***first executable statement  fftdoc
       return
      end
