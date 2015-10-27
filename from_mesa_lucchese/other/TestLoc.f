      Program TestLoc

      implicit none

      integer iv(2), i1, i2
      real*8 ir(2)

      i1 = loc(iv(1))
      i2 = loc(iv(2))

      write (6, "('for integer')")
      write (6, "('i1 =', i20)") i1
      write (6, "('i2 =', i20)") i2
      write (6, "('diff =', i10)") i2-i1
      write (6, "('diff by 4 =', i10)") (i2-i1)/4

      i1 = loc(ir(1))
      i2 = loc(ir(2))

      write (6, "('for real*8')")
      write (6, "('i1 =', i20)") i1
      write (6, "('i2 =', i20)") i2
      write (6, "('diff =', i10)") i2-i1
      write (6, "('diff by 4 =', i10)") (i2-i1)/4


      end

      
