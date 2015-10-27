*deck @(#)ffden.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            driver for free-free densities
c***                 
c***references       
c
c***routines called  
c***end prologue       m7001
      subroutine ffden(frea,freb,rho,val,nsymf,nlm,npts,nchan,sze,
     1                 npair,dimsym,typea,typeb,tst)
      implicit integer (a-z)
      character*16 type
      character*(*) typea, typeb
      character*8 ylm, ylm1
      real*8 rho, val
      complex*16 frea, freb
      logical tst, typlog, chnlog, symlog
      dimension frea(npts,sze), freb(npts,sze), nsymf(dimsym,nchan)
      dimension  rho(*), val(*), nlm(nchan)
      common /io/ inp, iout
      oldpr=0
      do 10 ic=1,nchan
         do 20 jc=1,ic
            call nprmax(nsymf(1,ic),nsymf(1,jc),typea,
     1                  typeb,dimsym,oldpr) 
   20    continue
   10 continue
      npr=sze*sze
      typlog=typea.eq.typeb
      if (typlog) then
          npr=sze*(sze+1)/2
      endif
      ijmax=npair/2
      filsiz=npr*2*npts
      call iosys ('create real "'//typea//'-'//typeb//' densities" '//
     1            'on atomci',filsiz,0,0,' ')
      write(iout,*) typea,' ',typeb,' density'
      prcnt=0
      cntic=1
      do 30 ic=1,nchan
         upprc=nchan
         if (typlog) then
             upprc=ic
         endif 
         cntjc=1
         do 40 jc=1,upprc
            chnlog=ic.eq.jc
            cnta=cntic 
            do 50 is=1,dimsym
               if (nsymf(is,ic).ne.0) then
                   cntb=cntjc
                   upprs=dimsym
                   if (typlog.and.chnlog) then
                        upprs=is
                   endif
                   call lmtyp(is,l,m,ylm)
                   islen=length(ylm)
                   do 60 js=1,upprs
                      type='non-diagonal'
                      if (nsymf(js,jc).ne.0) then
                          symlog=is.eq.js
                          if (typlog.and.chnlog.and.symlog) then
                              type='diagonal'
                          endif
                          call lmtyp(js,l,m,ylm1)
                          jslen=length(ylm1)
                          call mkden(frea(1,cnta),frea(1,cnta),
     1                               typea,freb(1,cntb),freb(1,cntb),
     2                               typeb,rho,rho,nsymf(is,ic),
     3                               nsymf(js,jc),ijmax,type,prcnt,npts)
                          write(iout,*) '   finished ',ylm(1:islen),
     1                                  '*',ylm1(1:jslen),' for '//
     2                                  'channels ic = ',ic,' jc = ',jc
                          cntb=cntb+nsymf(js,jc)
                      endif  
   60              continue
                   cnta=cnta+nsymf(is,ic) 
               endif
   50       continue
            cntjc=cntjc+nlm(jc)
   40    continue
         cntic=cntic+nlm(ic)
   30 continue
      if (prcnt.ne.npr) then
          call lnkerr('error in pair count for ff case')
      endif
      return
      end









