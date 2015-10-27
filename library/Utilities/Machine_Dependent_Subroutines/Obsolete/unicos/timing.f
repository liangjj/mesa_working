*deck @(#)timing.f	5.1  11/6/94
      subroutine timing(t1,t2,t3)
c
      real*8 t1,t2,t3
c
      t1=second()
      t2=timef()/1000.0
c
      return
      end
