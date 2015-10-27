*deck @(#)unitwt.f	4.1  7/7/93
      function unitwt(name)
c
c***begin prologue     unitwt
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)unitwt.f	4.1   7/7/93
c***purpose            to wait for asynchronous i/o to finish.
c
c***description        #
c
c
c***references
c
c***routines called    unit   (cftlib)
c
c   common blocks:     (none)
c
c***end prologue       unitwt
c
      implicit integer (a-z)
      integer unitwt
c
      unitwt=-1
cps      unitwt=unit(name)
c
c
      return
      end
