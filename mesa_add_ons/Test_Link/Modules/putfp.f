*deck @(#)putfp.f	5.1  11/6/94
      subroutine putfp(x,n,string,cursor)
c***begin prologue     putfp
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           character, string, floating point
c***author             martin, richard (lanl)
c***source             @(#)putfp.f	5.1   11/6/94
c***purpose            converts a floating point number into a string.
c***description
c                      call putfp(x,n,string,cursor)
c                        x       the floating point number.
c                        n       the number of places to keep after the
c                                decimal point.
c                        string  the string in which to insert the number.
c                        cursor  the number is inserted into string
c                                at index cursor+1.  the cursor is updated
c                                after the insertion.
c
c***references
c***routines called    fptoc(chr)
c***end prologue       putfp
      implicit integer(a-z)
      character*(*) string
      character*16 fptoc,cnum
      real*8 x,x1,ten
      parameter (ten=1.0d+01)
c
c
c     round the number off at the nth decimal point.
      x1=x+sign(0.5*ten**(-n),x)
c     convert the number.
      cnum=fptoc(x1)
c     find the decimal point.
      id=index(cnum,'.')
      bcur=cursor+1
      ecur=cursor+id+n
      string(bcur:ecur)=cnum(1:id+n)
      cursor=ecur
c
c
      return
      end
