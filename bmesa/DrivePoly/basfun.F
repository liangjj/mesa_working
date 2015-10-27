*deck basfun.f 
c***begin prologue     basfun
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             basfun
c***purpose            calculate the basis set from the dvr functions.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       basfun
      subroutine basfun(pt,wpt,cp,dcp,ddcp,wtt,dwtt,ddwtt,
     1                  keep,n,nkeep,prn)
c
      implicit integer (a-z)
      real*8 pt, wpt, cp, dcp, ddcp, wtt, dwtt, ddwtt
      real*8 scr
      logical prn
      character*80 title
      dimension pt(n), wpt(n), cp(n,n), dcp(n,n), ddcp(n,n)
      dimension wtt(n), dwtt(n), ddwtt(n)
      dimension keep(nkeep)
      common/io/inp, iout
      pointer (pscr,scr(1))
      if(nkeep.ne.n) then
         need=wptoin(nkeep*nkeep)
         pt1=1
         wpt1=pt1+nkeep  
         cp1=wpt1+nkeep
         dcp1=cp1+nkeep*nkeep
         ddcp1=dcp1+nkeep*nkeep
         wtt1=ddcp1+nkeep*nkeep
         dwtt1=wtt1+nkeep
         ddwtt1=dwtt1+nkeep
         call getmem(need,pscr,ngot,'basfun',0)
         call trimv(pt,pt,scr,keep,nkeep,n)
         call trimv(wpt,wpt,scr,keep,nkeep,n)
         call trimm(cp,cp,scr,keep,nkeep,n)
         call trimm(dcp,dcp,scr,keep,nkeep,n)
         call trimm(ddcp,ddcp,scr,keep,nkeep,n)
         call trimv(wtt,wtt,scr,keep,nkeep,n)
         call trimv(dwtt,dwtt,scr,keep,nkeep,n)
         call trimv(ddwtt,ddwtt,scr,keep,nkeep,n)
         call getmem(-ngot,pscr,idum,'basfun',idum)
      endif 
      if(prn) then
         title='trimmed coordinates'
         call prntrm(title,pt,nkeep,1,nkeep,1,iout)
         title='trimmed weights'
         call prntrm(title,wpt,nkeep,1,nkeep,1,iout)
         title='trimmed coordinate functions'
         call prntrm(title,cp,nkeep,nkeep,nkeep,nkeep,iout)
         title='trimmed first derivative of coordinate functions'
         call prntrm(title,dcp,nkeep,nkeep,nkeep,nkeep,iout)
         title='trimmed second derivative of coordinate functions'
         call prntrm(title,ddcp,nkeep,nkeep,nkeep,nkeep,iout)
         title='trimmed weight function function'
         call prntrm(title,wtt,nkeep,1,nkeep,1,iout)
         title='trimmed first derivative of weight function function'
         call prntrm(title,dwtt,nkeep,1,nkeep,1,iout)
         title='trimmed second derivative of weight function function'
         call prntrm(title,ddwtt,nkeep,1,nkeep,1,iout)
      endif
      return
      end


















