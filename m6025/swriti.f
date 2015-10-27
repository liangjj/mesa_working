*deck swriti
      subroutine swriti ( String, Number )
c
c---------------------------------- SWritI ----------------------------------
c
c >>> SWritI	-- string-write-integer; inserts a number in the string
c
c----------------------------------------------------------------------------
c
      character*(*) String
      integer       Number
c
      integer		ic1, ic2
c
      character*8	IForm
c
c----------------------------------------------------------------------------
c
      ic1 = index(String,'{')
      ic2 = index(String,'}')
      if ( ic1 .eq. 0 ) RETURN
      if ( ic2 .lt. ic1 ) RETURN
c
      write ( IForm, '(''(i'',i2,'')'')' ) ic2-ic1+1
c
      write ( String(ic1:ic2), IForm ) Number
c
      return
      end
