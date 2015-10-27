*deck @(#)cfden.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            driver for complex-free densities
c***                 
c***references       
c
c***routines called  
c***end prologue       m7001
      subroutine cfden(fnsc,fre,rho,val,nsymc,nsymf,npts,nb,
     1                 nchan,sze,npair,dimsym,typ,tst)
      implicit integer (a-z)
      character*8 typea, typeb, ylm, ylm1
      character*16 type
      character*(*) typ
      real*8 rho, val
      complex*16 fnsc, fre
      logical tst
      dimension fnsc(npts,nb), rho(*), val(*), nsymc(dimsym)
      dimension fre(npts,sze), nsymf(dimsym,nchan)
      common /io/ inp, iout
      typea='complex'
      ilen=length(typea)
      typeb=typ
      oldpr=0   
      do 10 nc=1,nchan
         call nprmax(nsymc,nsymf(1,nc),typea,typeb,dimsym,oldpr) 
   10 continue
      npr=sze*nb
      ijmax=npair/2
      filsiz=2*npr*npts
      type='non-diagonal'
      typeb=typ
      jlen=length(typeb)
      call iosys ('create real "'//typea(1:7)//'-'//typeb(1:5)//
     1            ' densities" on atomci',filsiz,0,0,' ')
      write(iout,*) typea(1:ilen),' ',typeb(1:jlen),' density'
      prcnt=0
      ispt=1 
      do 20 is=1,dimsym
         if (nsymc(is).ne.0) then
             call lmtyp(is,l,m,ylm)
             islen=length(ylm)
             jspt=1
             do 30 ic=1,nchan
                do 40 js=1,dimsym
                   if (nsymf(js,ic).ne.0) then
                        call lmtyp(js,l,m,ylm1)
                        jslen=length(ylm1)
                        call mkden(fnsc(1,ispt),fnsc(1,ispt),
     1                             typea(1:7),fre(1,jspt),
     2                             fre(1,jspt),typeb(1:5),rho,
     3                             rho,nsymc(is),nsymf(js,ic),
     4                             ijmax,type,prcnt,npts)
                        write(iout,*) '   finished ',ylm(1:islen),
     1                                '*',ylm1(1:jslen),
     2                                ' for channel ',ic
                        jspt=jspt+nsymf(js,ic)
                   endif
   40           continue
   30        continue
             ispt=ispt+nsymc(is)
         endif
   20 continue 
      if (prcnt.ne.npr) then
          call lnkerr('error in pair count for cf case')
      endif
      call iosys ('rewind "'//typea(1:7)//'-'//typeb(1:5)//
     1            ' densities" on atomci read-and-write',0,0,0,' ')
      if (tst) then
          do 50 is=1,dimsym
             if(nsymc(is).ne.0) then
                call lmtyp(is,l,m,ylm)
                islen=length(ylm)
                do 60 ic=1,nchan
                   do 70 js=1,dimsym
                      nn=nsymc(is)*nsymf(js,ic)
                      if (nsymf(js,ic).ne.0) then
                          if (is.eq.js) then
                              call lmtyp(js,l,m,ylm1)
                              jslen=length(ylm1)
                              write(iout,*) '   testing ',
     1                                      ylm(1:islen),'*',
     2                                      ylm1(1:jslen),
     3                                      ' density integral for '//
     4                                      'channel ',ic
                              call rhotst(rho,typea(1:7),rho,
     1                                    typeb(1:5),val,val,nn,npts)
                          else
                              call iosys ('read real "'//typea(1:7)
     1                                    //'-'//typeb(1:5)//
     2                                    ' densities" from '//
     3                                    'atomci without '//
     4                                    'rewinding',2*nn*npts,
     5                                     rho,0,' ')
                          endif
                      endif
   70              continue
   60           continue
             endif
   50     continue
      endif
      return
      end









