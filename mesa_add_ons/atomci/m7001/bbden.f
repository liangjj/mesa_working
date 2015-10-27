*deck @(#)bbden.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            driver for bound-bound densities
c***                 
c***references       
c
c***routines called  
c***end prologue       m7001
      subroutine bbden(fns,rho,val,nsym,npts,nf,npair,dimsym,tst)
      implicit integer (a-z)
      character*8 typea, ylm, ylm1
      character*16 type
      real*8 fns, rho, val
      logical tst
      dimension fns(npts,nf), rho(*), val(*), nsym(dimsym)
      common /io/ inp, iout
      typea='bound'
      oldpr=0   
      call nprmax(nsym,nsym,typea,typea,dimsym,oldpr) 
      npr=nf*(nf+1)/2
      ijmax=npair
      filsiz=npr*npts
      call iosys ('create real "bound-bound densities" on atomci',
     1             filsiz,0,0,' ')
      ilen=length(typea)
      write(iout,*) typea(1:ilen),' ',typea(1:ilen),' density'
      prcnt=0
      ispt=1 
      do 10 is=1,dimsym
         if (nsym(is).ne.0) then
             type='non-diagonal'
             call lmtyp(is,l,m,ylm)
             islen=length(ylm)
             jspt=1
             do 20 js=1,is
                if (nsym(js).ne.0) then
                    if (is.eq.js) then
                        type='diagonal' 
                    endif
                    call lmtyp(js,l,m,ylm1)
                    jslen=length(ylm1)
                    call mkden(fns(1,ispt),fns(1,ispt),'bound',
     1                         fns(1,jspt),fns(1,jspt),'bound',rho,
     2                         rho,nsym(is),nsym(js),
     3                         ijmax,type,prcnt,npts)
                    jspt=jspt+nsym(js)
                    write(iout,*) '   finished ',ylm(1:islen),
     1                            '*',ylm1(1:jslen)
                endif  
   20        continue
             ispt=ispt+nsym(is)
         endif
   10 continue
      if (prcnt.ne.npr) then
          call lnkerr('error in pair count for bb case')
      endif
      call iosys ('rewind "bound-bound densities" on atomci '//
     1            'read-and-write',0,0,0,' ')
      if (tst) then
          do 30 is=1,dimsym
             if(nsym(is).ne.0) then
                nn=nsym(is)*(nsym(is)+1)/2
                call lmtyp(is,l,m,ylm)
                islen=length(ylm)
                do 40 js=1,is
                   if (nsym(js).ne.0) then
                       call lmtyp(js,l,m,ylm1)
                       jslen=length(ylm1)
                       if (is.eq.js) then
                           write(iout,*) '   testing ',ylm(1:islen),
     1                                   '*',ylm1(1:jslen),
     2                                   ' density integral'
                           call rhotst(rho,'bound',rho,'bound',
     1                                 val,val,nn,npts)
                       else
                           do 50 ip=1,nsym(is)
                              do 60 jp=1,nsym(js)
                                 call iosys ('read real "bound-bound '//
     1                                       'densities" from atomci '//
     2                                       'without rewinding',npts,
     3                                        rho,0,' ') 
   60                         continue
   50                      continue
                       endif
                   endif
   40           continue
             endif
   30     continue
      endif
      return
      end









