*deck @(#)tomast.f	5.1  11/6/94
      subroutine tomast(master,local,cosine,origin)
c
      implicit integer (a-z)
c
      real*8 master(3),local(3),cosine(3,3),origin(3)
c
      call vmove(master,origin,3)
      call apbc(master,cosine,local,3,3,1)
c
c
      return
      end
