*deck @(#)vecout.f	5.1  11/6/94
      subroutine vecout(c,eigval,num)
c***begin prologue     vecout
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, output
c***author             saxe, paul (lanl)
c***source
c***purpose            prints eigenvectors and associated eigenvalues.
c***description
c                      call vecout(c,eigval,num)
c                        c         input eigenvector matrix (num,num).
c                        eigval    input eigenvalues (num).
c                        num       matrix dimension.
c
c***references
c***routines called    min
c***end prologue       vecout
      implicit integer (a-z)
c
      real*8 c(num,num),eigval(num)
      character*11 line
c
      common /io/     inp,iout
c
      data line /'-----------'/
      save line
c
      mx=0
    1 continue
      mn=mx+1
      mx=min(mx+7,num)
      write (iout,2) (iq,iq=mn,mx)
    2 format (/,1x,10(i6,5x))
      write (iout,3) (line,iq=mn,mx)
    3 format (1x,7a11)
      do 5 i=1,num
         write (iout,4) (c(i,iq),iq=mn,mx)
    4    format (1x,7f11.6)
    5 continue
      write (iout,6) (eigval(iq),iq=mn,mx)
    6 format (/,1x,7f11.6)
c
      if (mx.lt.num) go to 1
c
      return
      end
