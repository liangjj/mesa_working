*deck @(#)skipbl.f	5.1  11/6/94
      subroutine skipbl
c
      implicit character*1 (a-z)
c
      common /const/  blank, digit(10)
c
    1 if (getchr(char).eq.blank) go to 1
      call bkspac
      return
      end
