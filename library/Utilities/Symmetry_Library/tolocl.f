*deck @(#)tolocl.f	5.1  11/6/94
      subroutine tolocl(local,master,cosine,origin)
c
      implicit integer (a-z)
c
      real*8 local(3),master(3),cosine(3,3),origin(3),temp(3)
c
      do 1 i=1,3
         temp(i)=master(i)-origin(i)
    1 continue
c
      call ebc(local,temp,cosine,1,3,3)
c
c
      return
      end
