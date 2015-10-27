*deck @(#)constd.f	5.1  11/6/94
      block data constd
c
      implicit character*1 (a-z)
c
      common /const/  blank, digit(10)
c
      data blank /' '/
      data digit /'0','1','2','3','4','5','6','7','8','9'/
      end
