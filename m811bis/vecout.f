*deck @(#)vecout.f	1.1  11/30/90
      subroutine vecout(c,eigval,num)
c
      implicit integer (a-z)
c
      real*8 c(num,num),eigval(num)
      character*11 line
c
      common /io/     inp,ioutpt
c
      data line /'-----------'/
c
      mx=0
    1 continue
      mn=mx+1
      mx=min(mx+7,num)
      write (ioutpt,2) (iq,iq=mn,mx)
    2 format (/,1x,10(i6,5x))
      write (ioutpt,3) (line,iq=mn,mx)
    3 format (1x,7a11)
      do 5 i=1,num
         write (ioutpt,4) (c(i,iq),iq=mn,mx)
    4    format (1x,7f11.6)
    5 continue
      write (ioutpt,6) (eigval(iq),iq=mn,mx)
    6 format (/,1x,7f11.6)
c
      if (mx.lt.num) go to 1
c
      return
      end
