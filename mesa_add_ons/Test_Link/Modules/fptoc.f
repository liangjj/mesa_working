*deck @(#)fptoc.f	5.1  11/6/94
      function fptoc(fpnum)
c***begin prologue     fptoc
c***date written       850601  (yymmdd)
c***revision date      870131  (yymmdd)
c
c   31 january 1987  pws at lanl
c
c      changing call to amod to generic mod to suit the sun f77.
c
c***keywords           character, floating point, conversion
c***author             martin, richard (lanl)
c***source             @(#)fptoc.f	5.1   11/6/94
c***purpose            converts a floating point number into a string.
c***description
c                      fptoc is a character function used as:
c                        string=fptoc(fpnum)
c                          fpnum  the floating point number to convert.
c
c                      the string returned has the form smmm.nnnnn, there is
c                      no provision for exponential notation.
c                      the decimal portion past the 12th digit is padded
c                      with zeroes.
c
c***references
c***routines called    (none)
c***end prologue       fptoc
      implicit integer(a-z)
      character*16 fptoc
      character str*16
      character itoc*16
      real*8 fpnum,one,fract,ten,zero,absfp,round,rfp
      parameter (zero=0.0d+00,one=1.0d+00,ten=10.0d+00,round=0.5d-12)
      parameter (maxsiz=16,maxdec=12)
c
c
c     round to maximum precision.
      rfp=fpnum+sign(round,fpnum)
      absfp=abs(rfp)
c
c     get the portion before the decimal point, and its length.
      int=aint(rfp)
      if(rfp.le.zero) then
         str='-'//itoc(int)
      else
         str=itoc(int)
      endif
      lenstr=index(str,' ')
      if(lenstr.gt.0.and.lenstr.le.maxsiz) str(lenstr:lenstr)='.'
c
c     get the fractional portion.
      if(lenstr.lt.maxsiz) then
         fract=absfp
         id=0
   10    id=id+1
            fract=ten*mod(fract,one)
            if(id.gt.maxdec) then
               n=0
            else
               n=aint(fract)
            endif
            str(lenstr+1:)=itoc(n)
            lenstr=lenstr+1
            if(lenstr.lt.maxsiz) goto 10
      endif
c
c
      fptoc=str
c
c
      return
      end
