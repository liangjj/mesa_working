*deck @(#)ccden.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            driver for complex-complex densities
c***                 
c***references       
c
c***routines called  
c***end prologue       m7001
      subroutine ccden(fnsc,rho,val,nsymc,npts,nc,npair,dimsym,tst)
      implicit integer (a-z)
      character*8 typea, ylm, ylm1
      character*16 type
      real*8 rho, val
      complex*16 fnsc
      logical tst
      dimension fnsc(npts,nc), rho(*), val(*), nsymc(dimsym)
      common /io/ inp, iout
      typea='complex'
      oldpr=0   
      call nprmax(nsymc,nsymc,typea,typea,dimsym,oldpr) 
      npr=nc*(nc+1)/2
      ijmax=npair/2
      filsiz=2*npr*npts
      call iosys ('create real "complex-complex densities" on atomci',
     1             filsiz,0,0,' ')
      ilen=length(typea)
      write(iout,*) typea(1:ilen),' ',typea(1:ilen),' density'
      prcnt=0
      ispt=1 
      do 10 is=1,dimsym
         if (nsymc(is).ne.0) then
             type='non-diagonal'
             call lmtyp(is,l,m,ylm)
             islen=length(ylm)
             jspt=1
             do 20 js=1,is
                if (nsymc(js).ne.0) then
                    if (is.eq.js) then
                        type='diagonal' 
                    endif
                    call lmtyp(js,l,m,ylm1)
                    jslen=length(ylm1)
                    call mkden(fnsc(1,ispt),fnsc(1,ispt),'complex',
     1                         fnsc(1,jspt),fnsc(1,jspt),'complex',rho,
     2                         rho,nsymc(is),nsymc(js),
     3                         ijmax,type,prcnt,npts)
                    jspt=jspt+nsymc(js)
                    write(iout,*) '   finished ',ylm(1:islen),
     1                            '*',ylm1(1:jslen)
                endif  
   20        continue
             ispt=ispt+nsymc(is)
         endif
   10 continue
      if (prcnt.ne.npr) then
          call lnkerr('error in pair count for cc case')
      endif
      call iosys ('rewind "complex-complex densities" on atomci '//
     1            'read-and-write',0,0,0,' ')
      if (tst) then
          do 30 is=1,dimsym
             if(nsymc(is).ne.0) then
                nn=nsymc(is)*(nsymc(is)+1)/2
                call lmtyp(is,l,m,ylm)
                islen=length(ylm)
                do 40 js=1,is
                   if (nsymc(js).ne.0) then
                       call lmtyp(js,l,m,ylm1)
                       jslen=length(ylm1)
                       if (is.eq.js) then
                           write(iout,*) '   testing ',ylm(1:islen),
     1                                   '*',ylm1(1:jslen),
     2                                   ' density integral'
                           call rhotst(rho,'complex',rho,'complex',
     1                                 val,val,nn,npts)
                       else
                           do 50 ip=1,nsymc(is)
                              do 60 jp=1,nsymc(js) 
                                 call iosys ('read real "complex-'//
     1                                       'complex densities" '//
     2                                       'from atomci without '//
     3                                       'rewinding',2*npts,
     4                                        rho,0,' ')
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









