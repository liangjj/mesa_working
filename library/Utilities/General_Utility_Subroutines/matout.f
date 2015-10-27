*deck @(#)matout.f	5.1  11/6/94
      subroutine matout(a,nad,nbd,m,n,iout)
c***begin prologue     matout
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix output, matrix print
c***author             saxe, paul (lanl)
c***source             @(#)matout.f	5.1   11/6/94
c***purpose
c                      sends a matrix to the output file.
c***description
c                      call matout(a,nad,nbd,m,n,iout)
c                        a        input matrix dimensioned (nad,nbd).
c                        m        number of rows to print.
c                        n        number of columns to print.
c                        iout     unit number of output file.
c
c***references
c***routines called    (none)
c***end prologue       matout
      implicit real*8 (a-h,o-z)
      character*6 line
      dimension a(nad,nbd)
    1 format(2x,10(7x,i5))
    2 format(2x,21a6)
    3 format(2x,i2,2x,10f12.7)
    4 format(/)
      data line / '------' /
      save line
      ii=0
      jj=0
  200 ii=ii+1
      jj=jj+1
      kk=5*jj
      nn=n
      if(n.gt.kk) nn=kk
      ll=2*(nn-ii+1)+1
      write (iout,1) (i,i=ii,nn)
      write (iout,2) (line,i=1,ll)
      do 101 i=1,m
      write (iout,3) i,(a(i,j),j=ii,nn)
  101 continue
      if(n.le.kk) go to 201
      write (iout,4)
      ii=kk
      go to 200
  201 return
      end
