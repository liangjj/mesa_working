*deck @(#)getlin.f	5.1  11/6/94
      subroutine getlin
c
      implicit integer (a-z)
c
      character*1 image
      character*80 tmpbuf
c
      common /tapes/  out,errout,input,drttap
      common /buffer/ pos
      common /buffc/  image(80)
c
      pos=0
      read (input,1,end=2) tmpbuf
    1 format (a80)
      call locase(tmpbuf,tmpbuf)
      do 100 i=1,80
         image(i)=tmpbuf(i:i)
  100 continue
      return
    2 write (errout,3) input
    3 format (//,' encountered an end to tape',i2,' while scanning '
     #,          'for orbital codes')
      call lnkerr('found end of input while still looking')
      end
