*deck @(#)print.f	5.1  11/6/94
      subroutine print(a,nad,m,iout)
c***begin prologue     print
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix output, matrix print
c***author             saxe, paul (lanl)
c***source             @(#)print.f	5.1   11/6/94
c***purpose
c                      sends a vector to the output file.
c***description
c                      call print(a,nad,m,iout)
c                        a        input vector of length nad.
c                        m        number of elements per line of output.
c                        iout     unit number of the output file.
c
c***references
c***routines called    (none)
c***end prologue       print
      implicit real*8 (a-h,o-z)
      character*6 line
      dimension a(nad)
1     format(2x,10(7x,i5))
2     format(2x,21a6)
3     format(2x,i2,2x,10f12.7)
4     format(/)
      data line / '------' /
      save line
      ii=0
      jj=0
  200 ii=ii+1
      jj=jj+1
      kk=6*jj
      nn=kk+kk*(kk-1)/2
      mm=m
      if(m.gt.kk) mm=kk
      ll=2*(mm-ii+1)+1
      write (iout,1) (i,i=ii,mm)
      write (iout,2) (line,i=1,ll)
      do 101 i=ii,m
      i1=i*(i-1)/2+ii
      i2=i+i*(i-1)/2
      if(i2.gt.nn) i2=i1+5
      write (iout,3) i,(a(j),j=i1,i2)
  101 continue
      if(m.le.kk) go to 201
      write (iout,4)
      ii=kk
      go to 200
  201 return
      end
