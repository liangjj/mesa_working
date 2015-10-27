*deck @(#)hinout.f	5.1  11/6/94
      subroutine hinout(row,nnci)
c
c  this routine reads in the hamiltonian matrix (hin) from iunt2a
c  and writes it out (hout) on iunt2a
c
      implicit real*8 (a-h,o-z)
      dimension row(nnci)
c
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
c
      entry hin(row,nnci)
      read (iunt2a) row
c
c      write(iw,*)'hin: nnci,row',nnci,row
c
      return
c
      entry hout(row,nnci)
      write (iunt2a) row
c
c      write(iw,*)'hout: nnci,row',nnci,row
c
      return
c
      end
