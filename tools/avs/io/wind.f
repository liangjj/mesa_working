*deck  @(#)wind.f	4.1   7/7/93
      subroutine wind(file)
c
c***begin prologue     wind
c***date written       870708   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           winding a file
c***author             saxe, paul (lanl)
c***source             @(#)wind.f	4.1   7/7/93
c
c***purpose            to position a formatted file at the end of
c     information.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       wind
c
      implicit integer (a-z)
c
      rewind (unit=file)
c
 1    continue
         read (file,2,end=3)
 2       format(/)
c
 
      go to 1
c
 3    continue
c
      backspace file
c
c
      return
      end
