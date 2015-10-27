*deck dohtrl
      subroutine dohtrl (q, n, nrda, diag, irank, div, td)
c***begin prologue  dohtrl
c***subsidiary
c***purpose  subsidiary to dbvsup and dsuds
c***library   slatec
c***type      double precision (ohtrol-s, dohtrl-d)
c***author  watts, h. a., (snla)
c***description
c
c     for a rank deficient problem, additional orthogonal
c     householder transformations are applied to the left side
c     of q to further reduce the triangular form.
c     thus, after application of the routines dorthr and dohtrl
c     to the original matrix, the result is a nonsingular
c     triangular matrix while the remainder of the matrix
c     has been zeroed out.
c
c***see also  dbvsup, dsuds
c***routines called  ddot
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dohtrl
      double precision ddot
      integer irank, irp, j, k, kir, kirm, l, n, nmir, nrda
      double precision dd, diag(*), diagk, div(*), q(nrda,*), qs, sig,
     1     sqd, td(*), tdv
c***first executable statement  dohtrl
      nmir = n - irank
      irp = irank + 1
      do 40 k = 1, irank
         kir = irp - k
         diagk = diag(kir)
         sig = (diagk*diagk) + ddot(nmir,q(irp,kir),1,q(irp,kir),1)
         dd = sign(sqrt(sig),-diagk)
         div(kir) = dd
         tdv = diagk - dd
         td(kir) = tdv
         if (k .eq. irank) go to 30
            kirm = kir - 1
            sqd = dd*diagk - sig
            do 20 j = 1, kirm
               qs = ((tdv*q(kir,j))
     1               + ddot(nmir,q(irp,j),1,q(irp,kir),1))/sqd
               q(kir,j) = q(kir,j) + qs*tdv
               do 10 l = irp, n
                  q(l,j) = q(l,j) + qs*q(l,kir)
   10          continue
   20       continue
   30    continue
   40 continue
      return
      end
