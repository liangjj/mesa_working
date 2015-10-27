      Program TestShift

      implicit none
      integer val, i, j, k, l, valnew
      integer ip, jp, kp, lp
      integer mask, nshift
      integer iout, jout, kout, lout
      integer shiftl, shiftr
      external shiftl, shiftr
      val = 7

      mask = 1023
      nshift = 10

      DO 10 i = 0, 20
         valnew = shiftl(val, 3)
         write (6, 20) i, valnew
         val = valnew
 20      format("lshift i", i5, '  valnew ', O25.24)
 10   continue

      val = 7*(2**27)
      DO 30 i = 0, 20
         valnew = shiftr(val, 3)
         write (6, 20) i, valnew
         val = valnew
 40      format("rshift i", i5, '  valnew ', O25.24)
 30   continue

c ** test packing two integers
      write (6, 45)
 45   format('Testing two packing')
      do 50 ip = 1, nshift
         i = 2**(ip-1)
         do 60 jp = 1, nshift
            j = 2**(jp-1)
            val = or(shiftl(i, nshift), j)
            jout=and(val, mask)
            iout=shiftr(val, nshift)
            if (jout .ne. j .or. iout .ne. i) THEN
               write (6, 70) i, j, val, iout, jout
 70            format('i', O5.4, '  j', O5.4, /,
     $         'val', O25.24, /,
     $         'iout', O5.4, '  jout', O5.4)
               stop
            end if
 60      continue
 50   continue

c ** test packing three integers
      write (6, 345)
 345  format('Testing three packing')
      do 350 ip = 1, nshift
         i = 2**(ip-1)
         do 360 jp = 1, nshift
            j = 2**(jp-1)
            do 380 kp = 1, nshift
               k = 2**(kp-1)
               val = or(shiftl(i, nshift), j)
               val = or(shiftl(val, nshift), k)
               kout=and(val, mask)
               jout=and(shiftr(val, nshift), mask)
               iout=shiftr(val, 2*nshift)
               if (kout .ne. k .or. jout .ne. j .or. iout .ne. i) THEN
                  write (6, 370) i, j, k, val, iout, jout, kout
 370              format('i', O5.4, '  j', O5.4, '  k', O5.4, /,
     $                 'val', O25.24, /,
     $                 'iout', O5.4, '  jout', O5.4, '  kout', O5.4)
                  stop
               end if
 380        continue
 360     continue
 350  continue


c ** test packing four integers
      write (6, 445)
 445  format('Testing four packing')
      do 450 ip = 1, nshift
         i = 2**(ip-1)
         do 460 jp = 1, nshift
            j = 2**(jp-1)
            do 480 kp = 1, nshift
               k = 2**(kp-1)
               do 490 lp = 1, nshift
                  l = 2**(lp-1)
                  val = or(shiftl(i, nshift), j)
                  val = or(shiftl(val, nshift), k)
                  val = or(shiftl(val, nshift), l)
                  lout=and(val, mask)
                  kout=and(shiftr(val, nshift), mask)
                  jout=and(shiftr(val, 2*nshift), mask)
                  iout=shiftr(val, 3*nshift)
                  if (kout .ne. k .or. jout .ne. j 
     $                 .or. iout .ne. i .or. lout .ne. l) THEN
                     write (6, 470) i, j, k, l, val, 
     $                    iout, jout, kout, lout
 470                 format('i', O5.4, '  j', O5.4, '  k', O5.4,
     $                    '  l', O5.4/,
     $                    'val', O25.24, /,
     $                    'iout', O5.4, '  jout', O5.4, '  kout', O5.4,
     $                    '  lout', O5.4)
                     stop
                  end if
 490           continue
 480        continue
 460     continue
 450  continue
      
      end
      function shiftl(i,j)
c
      implicit integer(a-z)
      integer shiftl
c
c
      shiftl=lshift(i,j)
c
      return
      end
      function shiftr(i,j)
c
      implicit integer(a-z)
      integer shiftr
c
c
      shiftr=rshift(i,j)
c
c
      return
      end
