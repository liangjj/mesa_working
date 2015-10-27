*deck @(#)setta.f	5.1  11/6/94
      subroutine setta(ta,nob)
      implicit integer(a-z)
      real*8 ta(nob,nob)
c
      nob1=nob-1
      do 10 i=1,nob1
         i1=i+1
         do 11 j=i1,nob
            ta(j,i)=0.0d+00
  11     continue
 10   continue
c
      do 15 i=1,nob
         ta(i,i)=(.5d+00)*ta(i,i)
  15  continue
c
      do 20 i=1,nob
         do 25 j=1,nob
            ta(j,i)=-ta(j,i)
  25     continue
  20  continue
c
c
      return
      end
