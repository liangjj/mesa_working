*deck @(#)dcaptl.f	4.1  7/7/93
      function dcaptl(in)
c
c***begin prologue     dcaptl
c***date written       861201  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           string, lower case, decapitlizing
c***author             saxe, paul (lanl)
c***source             @(#)dcaptl.f	4.1   7/7/93
c***purpose            converts uppercase characters in a string to lowercase
c***description
c                      dcaptl is a character function used as:
c                      lower=dcaptl(string)
c                        string  the input chracter string.
c                        dcaptl  the output character string.
c***references
c***routines called    (none)
c***end prologue       dcaptl
c
      implicit integer(a-z)
c
      character*(*) dcaptl
      character*(*) in
      character*26 lower,upper
c
      data lower /'abcdefghijklmnopqrstuvwxyz'/
      data upper /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      save lower,upper
c
      dcaptl=' '
      do 1 i=1,min(len(in),len(dcaptl))
         loc=index(upper,in(i:i))
         if (loc.le.0) then
            dcaptl(i:i)=in(i:i)
         else
            dcaptl(i:i)=lower(loc:loc)
         end if
    1 continue
c
c
      return
      end
