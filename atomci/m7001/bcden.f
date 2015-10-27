*deck @(#)bcden.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            driver for bound-complex densities
c***                 
c***references       
c
c***routines called  
c***end prologue       m7001
      subroutine bcden(fns,fnsc,rho,val,nsymb,nsymc,npts,nb,nc,
     1                 npair,dimsym,tst)
      implicit integer (a-z)
      character*8 typea, typeb, ylm, ylm1
      character*16 type
      real*8 fns, rho, val
      complex*16 fnsc
      logical tst
      dimension fns(npts,nb), rho(*), val(*), nsymb(dimsym)
      dimension fnsc(npts,nc), nsymc(dimsym)
      common /io/ inp, iout
      typea='bound'
      typeb='complex'
      oldpr=0   
      call nprmax(nsymb,nsymc,typea,typeb,dimsym,oldpr) 
      npr=nb*nc
      ijmax=npair/2
      filsiz=2*npr*npts
      call iosys ('create real "bound-complex densities" on atomci',
     1             filsiz,0,0,' ')
      ilen=length(typea)
      jlen=length(typeb)
      write(iout,*) typea(1:ilen),' ',typeb(1:jlen),' density'
      type='non-diagonal'
      prcnt=0
      ispt=1 
      do 10 is=1,dimsym
         if (nsymb(is).ne.0) then
             call lmtyp(is,l,m,ylm)
             islen=length(ylm)
             jspt=1
             do 20 js=1,dimsym
                if (nsymc(js).ne.0) then
                    call lmtyp(js,l,m,ylm1)
                    jslen=length(ylm1)
                    call mkden(fns(1,ispt),fns(1,ispt),'bound',
     1                         fnsc(1,jspt),fnsc(1,jspt),'complex',rho,
     2                         rho,nsymb(is),nsymc(js),
     3                         ijmax,type,prcnt,npts)
                    jspt=jspt+nsymc(js)
                    write(iout,*) '   finished ',ylm(1:islen),
     1                            '*',ylm1(1:jslen)
                endif  
   20        continue
             ispt=ispt+nsymb(is)
         endif
   10 continue
      if (prcnt.ne.npr) then
          call lnkerr('error in pair count for bc case')
      endif
      call iosys ('rewind "bound-complex densities" on atomci '//
     1            'read-and-write',0,0,0,' ')
      if (tst) then
          do 30 is=1,dimsym
             if(nsymb(is).ne.0) then
                call lmtyp(is,l,m,ylm)
                islen=length(ylm)
                do 40 js=1,dimsym
                   nn=nsymb(is)*nsymc(js)
                   if (is.eq.js) then
                       if (nsymc(js).ne.0) then
                           call lmtyp(js,l,m,ylm1)
                           jslen=length(ylm1)

                           write(iout,*) '   testing ',ylm(1:islen),
     1                                   '*',ylm1(1:jslen),
     2                                   ' density integral'
                           call rhotst(rho,'bound',rho,'complex',
     1                                     val,val,nn,npts)
                       endif
                   else
                       call iosys ('read real "bound-complex '//
     1                             'densities" from atomci '//
     2                             'without rewinding',2*nn*npts,rho,
     3                              0,' ')
                   endif
   40           continue
             endif
   30     continue
      endif
      return
      end









