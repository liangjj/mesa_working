      subroutine test1
      real*8 a
#ifdef DECFLAG
      a=1.d0
      a=5.d0
#endif DECFLAG
      return
      end
