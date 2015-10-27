*deck swritd
      subroutine swritd ( String, Number )
c
c---------------------------------- SWritD ----------------------------------
c
c >>> SWritD	-- string-write-double; inserts a number in the string
c
c----------------------------------------------------------------------------
c
      character*(*) String
      real*8        Number
c
      integer		ic1, ic2, ic3
      integer		nf, nd
c
      character*12	IForm
c
c----------------------------------------------------------------------------
c
      ic1 = index(String,'{')
      ic2 = index(String,'}')
      if ( ic1 .eq. 0 ) RETURN
      if ( ic2 .lt. ic1 ) RETURN
c
      nf = ic2-ic1+1
c
      ic3 = index(String(ic1:ic2),'.')
      IF ( ic3 .gt. 0 ) THEN
        ic3 = ic3 + ic1-1
        nd = max(0,ic2 - ic3 - 5)
      ELSE
        nd = max(0,nf-7)
      END IF
c
      IForm = ' '
      write ( IForm, '(''(1pe'',i2,''.'',i2,'')'')' ) nf,nd
c
      write ( String(ic1:ic2), IForm ) Number
      IForm = String(ic1:ic2)
c
      return
      end
