*deck @(#)ctofp.f	5.1  11/6/94
      function ctofp(c)
c***begin prologue     ctofp
c***date written       850601  (yymmdd)
c***revision date      910304  (yymmdd)
c
c    3 march,1991      rlm at lanl
c      the variable lenc was never defined -- changed it to len(c).
c***keywords           character, floating point, conversion
c***author             saxe, paul (lanl)
c***source             @(#)ctofp.f	5.1   11/6/94
c***purpose            converts a character string into a floating point
c                      number.
c***description
c                      ctofp is a real function used as:
c                        fpnum=ctofp(c)
c                          c  the input character string.
c
c                      ctofp converts a string of the form smmmm.nnnndseee,
c                      the signs and exponents need not be present.
c
c***references
c***routines called    (none)
c***end prologue       ctofp
c
      implicit integer (a-z)
c
      real*8 ctofp
      character*(*) c,expon*4
c
      data expon /'ddee'/
      save expon
c
c
c     ----- find decimal point and exponent marker if present
c
      locdec=index(c,'.')
      locexp=0
      do 1 i=1,4
         locexp=max(locexp,index(c,expon(i:i)))
    1 continue
c
c     ----- integral part of value -----
c
      if (locdec.gt.0) then
         loc=locdec-1
      else if (locexp.gt.0) then
         loc=locexp-1
      else
         loc=len(c)
      end if
      ctofp=abs(ctoi(c(1:loc)))
c
c     ----- fractional part of value -----
c
      if (locdec.gt.0) then
         if (locexp.gt.0) then
            n=locexp-locdec-1
         else
            n=cskipb(c,' ')-locdec
         end if
         if (n.gt.0) then
            ctofp=ctofp+ctoi(c(locdec+1:locdec+n))/10.0d+00**n
         end if
      end if
c
c     ----- exponent part -----
c
      if (locexp.gt.0) then
         ctofp=ctofp*10.0**(ctoi(c(locexp+1:)))
      end if
c
c     ----- and sign -----
c
      if (c(1:1).eq.'-') ctofp=-ctofp
c
c
      return
      end
