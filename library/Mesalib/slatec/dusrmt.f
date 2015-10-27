*deck dusrmt
      subroutine dusrmt (i, j, aij, indcat, prgopt, dattrv, iflag)
c***begin prologue  dusrmt
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (usrmat-s, dusrmt-d)
c***author  (unknown)
c***description
c
c   the user may supply this code
c
c***see also  dsplp
c***routines called  (none)
c***revision history  (yymmdd)
c   811215  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  dusrmt
      double precision prgopt(*),dattrv(*),aij
      integer iflag(*)
c
c***first executable statement  dusrmt
      if(iflag(1).eq.1) then
c
c     this is the initialization step.  the values of iflag(k),k=2,3,4,
c     are respectively the column index, the row index (or the next col.
c     index), and the pointer to the matrix entry's value within
c     dattrv(*).  also check (dattrv(1)=0.) signifying no data.
           if(dattrv(1).eq.0.d0) then
           i = 0
           j = 0
           iflag(1) = 3
           else
           iflag(2)=-dattrv(1)
           iflag(3)= dattrv(2)
           iflag(4)= 3
           endif
c
           return
      else
           j=iflag(2)
           i=iflag(3)
           l=iflag(4)
           if(i.eq.0) then
c
c     signal that all of the nonzero entries have been defined.
                iflag(1)=3
                return
           else if(i.lt.0) then
c
c     signal that a switch is made to a new column.
                j=-i
                i=dattrv(l)
                l=l+1
           endif
c
           aij=dattrv(l)
c
c     update the indices and pointers for the next entry.
           iflag(2)=j
           iflag(3)=dattrv(l+1)
           iflag(4)=l+2
c
c     indcat=0 denotes that entries of the matrix are assigned the
c     values from dattrv(*).  no accumulation is performed.
           indcat=0
           return
      endif
      end
