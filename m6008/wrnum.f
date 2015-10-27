*deck @(#)wrnum.f	1.1 9/8/91
c***begin prologue     wrnum
c***date written       890528   (yymmdd)
c***revision date               (yymmdd)
c***keywords           write, bound
c***author             schneider, barry (lanl)
c***source             m6008
c***purpose            write out kohn numerator matrix
c*** 
c
c***references         none      
c
c***routines called    iosys
c***end prologue       wrnum
      subroutine wrnum(hambb,mxb)
      implicit integer (a-z)
      real *8 hambb
      dimension hambb(mxb,mxb)
      call iosys ('write real "bnd bnd num" to kohnint',mxb*mxb,hambb,
     1            0,' ')
      return
      end
