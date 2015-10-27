*deck @(#)fixta.f	5.1  11/6/94
      subroutine fixta(ta,nob)
      implicit integer(a-z)
      real*8 ta(nob,nob)
c
      nob1=nob-1
      do 1 i=1,nob1
         i1=i+1
         do 2 j=i1,nob
            ta(j,i)=0.0d+00
   2     continue
  1   continue
c
      do 10 i=1,nob
         do 20 j=1,i
            ta(j,i)=-ta(j,i)
  20     continue
 10   continue
c
      do 30 i=1,nob
         ta(i,i)=ta(i,i)*0.5d+00
  30  continue
c
      return
      end
