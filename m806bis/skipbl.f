*deck @(#)skipbl.f	1.2  7/30/91
      subroutine skipbl
c
      implicit character*1 (a-z)
c
      common /const/  blank, digit(10)
c
    1 if (getc(char).eq.blank) go to 1
      call bkspac
      return
      end
