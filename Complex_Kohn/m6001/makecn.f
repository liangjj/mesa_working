*deck @(#)makecn.f	1.1 9/7/91
c***begin prologue     makecn
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           makecn, link 6001, numerical, orbital
c***author             schneider, barry (lanl)
c***source             m6001
c***purpose            numerical orbital tabulation
c***description        calculates molecular orbitals on physical grid
c***                   for integral calculation
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       makecn
      subroutine makecn (fcon,grid,pre,rsq,anorm,ncon,nfmax,nptmx,
     1                   nwrite,nwds,reg,prnt)
      implicit integer (a-z)
      parameter (dimpr=300 , dimcen=10)
      real *8 grid, alf, cont, rloc, charge, fcon
      real *8  pre, anorm, rsq
      logical prnt
      common /io/ inp,iout
      common /aosi/ npr, ndum, nxyzc(dimpr,4), nprc(dimpr)
      common /aosr/ alf(dimpr), cont(dimpr)
      common /rloc/ rloc(3,dimcen), charge(dimcen)
      dimension fcon(nptmx,nfmax), grid(4,nptmx), anorm(*)
      dimension pre(nptmx), rsq(nptmx)
c-----------------------------------------------------------------------c
c                  set up file to hold                                  c
c                  contracted ao's on grid                              c
c                  nwbuf can hold at least one contracted               c
c                  function for one buffer load of points               c
c-----------------------------------------------------------------------c
c-----------------------------------------------------------------------c
c             calculate exponential factor associated                   c
c             with each primitive in this contracted                    c
c             function and use in all places                            c
c                       that it occurs                                  c 
c-----------------------------------------------------------------------c
c-----------------------------------------------------------------------c
c            loop over contracted functions                             c
c-----------------------------------------------------------------------c 
      wrdcnt=0
      call rzero(fcon,nptmx*nfmax)
      first=0
      last=0
      cntcon=0
      con1=0
      con2=0
      do 10 con=1,ncon
         cntcon=cntcon+1
         first=first+1
         last=last+nprc(con)
         cen=nxyzc(first,4)
c-----------------------------------------------------------------------c
c                   points loop for factors                             c
c                   depending only on contracted function               c
c-----------------------------------------------------------------------c
         do 20 i=1,nptmx
            pre(i)=1.d+00
   20    continue
         if (nxyzc(first,1).ne.0) then
             do 30 i=1,nptmx
                pre(i)=pre(i)*(grid(1,i)-rloc(1,cen))**nxyzc(first,1)
   30        continue
         endif
         if (nxyzc(first,2).ne.0) then
             do 40 i=1,nptmx
                pre(i)=pre(i)*(grid(2,i)-rloc(2,cen))**nxyzc(first,2)
   40        continue
         endif
         if (nxyzc(first,3).ne.0) then
             do 50 i=1,nptmx 
                pre(i)=pre(i)*(grid(3,i)-rloc(3,cen))**nxyzc(first,3) 
   50        continue
         endif
         do 60 i=1,nptmx
            rsq(i)=(grid(1,i)-rloc(1,cen))*(grid(1,i)-rloc(1,cen))+
     1          (grid(2,i)-rloc(2,cen))*(grid(2,i)-rloc(2,cen))+
     2          (grid(3,i)-rloc(3,cen))*(grid(3,i)-rloc(3,cen))
   60    continue
c-----------------------------------------------------------------------c
c              outer loop over primitives                               c
c              inner over points                                        c
c-----------------------------------------------------------------------c
         do 70 prim=first,last
            do 80 i=1,nptmx
               fcon(i,cntcon)=fcon(i,cntcon)+cont(prim)*exp(-alf(prim)*
     1                                                  rsq(i))
   80       continue         
   70    continue
         do 90 i=1,nptmx
            fcon(i,cntcon)=fcon(i,cntcon)*pre(i)
c-----------------------------------------------------------------------c
c                   calculate normalization as check                    c
c-----------------------------------------------------------------------c
            anorm(con)=anorm(con)+fcon(i,cntcon)*fcon(i,cntcon)*
     1                                           grid(4,i)
   90    continue
         wrdcnt=wrdcnt+nptmx
         if (cntcon.eq.nfmax) then
c-----------------------------------------------------------------------c
c                  print array if requested                             c
c-----------------------------------------------------------------------c
             if (prnt) then
                 con1=con1+1
                 con2=con2+nfmax
                 call wrtorb(fcon,reg,con1,con2,nfmax,nptmx)
                 con1=con2
             endif
c-----------------------------------------------------------------------c
c                the array is filled dump it out                        c
c-----------------------------------------------------------------------c
             nwrite=nwrite+1
             call iosys ('write real "con array" to orbs '//
     1                   'without rewinding',wrdcnt,fcon,0,' ')
             nwds=nwds+wrdcnt
             wrdcnt=0
             cntcon=0
             call rzero(fcon,nptmx*nfmax)
         endif
            first=last
   10 continue
c----------------------------------------------------------------------c
c                    get the dregs if any are left                     c
c----------------------------------------------------------------------c
      if (wrdcnt.ne.0) then
          nwrite=nwrite+1
          call iosys ('write real "con array" to orbs without '//
     1                'rewinding',wrdcnt,fcon,0,' ')
          nwds=nwds+wrdcnt
      endif
      if (prnt) then
          con1=con1+1
          con2=con2+nfmax
          con2=min(con2,ncon)
          call wrtorb(fcon,reg,con1,con2,nfmax,nptmx)
      endif
      return
      end
