*deck cmpre.f 
c***begin prologue     cmpre
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             cmpre
c***purpose            expansion of functions and derivatives in coordinate functions
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       cmpre
      subroutine cmpre(pt,wpt,cp,dcp,ddcp,type,parity,n)
c
      implicit integer (a-z)
      character*(*) type, parity
      real*8 pt, wpt, cp, dcp, ddcp, fun
      dimension pt(n), wpt(n), cp(n,n), dcp(n,n), ddcp(n,n)
      common/io/inp, iout
      pointer (pc,fun(1))
c
      fex=1
      dfex=fex+n
      ddfex=dfex+n
      fapp=ddfex+n
      dfapp=fapp+n
      ddfapp=dfapp+n
      c=ddfapp+n
      scr=c+n
      need=wpadti(scr+n*n)
      call getmem(need,pc,ngot,'cmpre',0)
      call chkfun(cp,dcp,ddcp,pt,wpt,fun(fex),fun(dfex),fun(ddfex),
     1            fun(fapp),fun(dfapp),fun(ddfapp),
     3            fun(c),fun(scr),type,parity,n)
      return
      end


















