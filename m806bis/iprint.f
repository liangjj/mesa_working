*deck @(#)iprint.f	1.2  7/30/91
      subroutine iprint(a,n,label)
      implicit integer (a-z)
      integer a(*)
      character*(*) label
      common /io/ inp,iout
      write (iout,1) (label(i:i),i=1,len(label))
    1 format (/,1x,79a1)
      mx=0
    2 continue
         mn=mx+1
         mx=min(n,mx+10)
         write (iout,3) (i,i=mn,mx)
    3    format(/,10x,10i6)
         write (iout,4) ('-',i=1,(mx-mn+1)*6)
    4    format(10x,60a1)
         do 6 i=mn,n
             ii=i*(i-1)/2
             write (iout,5) i,(a(ii+j),j=mn,min(i,mx))
    5        format (1x,i6,':   ',10i6)
    6    continue
      if (mx.lt.n) go to 2
c
      return
      end
