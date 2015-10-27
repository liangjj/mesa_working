*deck %W%  %G%
      subroutine timing(t1,t2,t3)
c
      call timeused(it1,it2,it3)
      t1=float(it1)/1.e6
      t2=float(it2)/1.e6
      t3=float(it3)/1.e6
      return
      end
