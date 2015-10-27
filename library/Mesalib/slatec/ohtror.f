*deck ohtror
      subroutine ohtror (q, n, nrda, diag, irank, div, td)
c***begin prologue  ohtror
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (ohtror-s)
c***author  watts, h. a., (snla)
c***description
c
c     for a rank deficient problem, additional orthogonal
c     householder transformations are applied to the right side
c     of q to further reduce the triangular form.
c     thus, after application of the routines orthol and ohtror
c     to the original matrix, the result is a nonsingular
c     triangular matrix while the remainder of the matrix
c     has been zeroed out.
c
c***see also  bvsup
c***routines called  sdot
c***revision history  (yymmdd)
c   750601  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  ohtror
      dimension q(nrda,*),diag(*),div(*),td(*)
c***first executable statement  ohtror
      nmir=n-irank
      irp=irank+1
      do 30 k=1,irank
         kir=irp-k
         diagk=diag(kir)
         sig=(diagk*diagk)+sdot(nmir,q(kir,irp),nrda,q(kir,irp),nrda)
         dd=sign(sqrt(sig),-diagk)
         div(kir)=dd
         tdv=diagk-dd
         td(kir)=tdv
         if (k .eq. irank) go to 30
         kirm=kir-1
         sqd=dd*diagk-sig
         do 20 j=1,kirm
            qs=((tdv*q(j,kir))+sdot(nmir,q(j,irp),nrda,q(kir,irp),nrda))
     1               /sqd
            q(j,kir)=q(j,kir)+qs*tdv
            do 10 l=irp,n
   10          q(j,l)=q(j,l)+qs*q(kir,l)
   20    continue
   30 continue
      return
      end
