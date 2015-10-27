*deck @(#)mkden.f	1.1 9/8/91
c***begin prologue     mkden
c***date written       920614   (yymmdd)
c***revision date               (yymmdd)
c***keywords           density
c***author             schneider, barry(nsf)
c***source             @(#)m7000
c***purpose            calculate densities for numerically
c***                   tabulated orbitals.
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine mkden(fnar,fnac,typa,fnbr,fnbc,typb,rhor,rhoc,
     1                 na,nb,np,type,prtot,n)
      implicit integer (a-z)
      real *8 fnar, fnbr, rhor
      complex *16 fnac, fnbc, rhoc
      character *(*) typa, typb, type
      character *128 filnme
      dimension fnar(n,na), fnbr(n,nb), fnac(n,na), fnbc(n,nb)
      dimension rhor(n,np), rhoc(n,np)
      common /io/ inp, iout
c**********************************************************************c
c               fnar=fnac and rhor=rhoc in the calling routine.        c
c               for this routine to work correctly it is necessary     c
c               for rhor and rhoc to be equivalenced for the write     c
c**********************************************************************c  
      filnme='"'//typa//'-'//typb//' densities"'
      cntwd=0
      prcnt=0
      cnt=0
      if (typa.eq.'bound') then
c**********************************************************************c
c                   bound-bound                                        c
c**********************************************************************c
          if (typb.eq.'bound') then
              if (type.eq.'diagonal') then
                  nprtot=na*(na+1)/2
                  chkwd=nprtot*n
                  do 10 i=1,na
                     do 20 j=1,i
                        cnt=cnt+1
                        do 30 k=1,n
                           rhor(k,cnt)=fnar(k,i)*fnar(k,j)
   30                   continue
                        if (cnt.eq.np) then
                            cntwd=cntwd+np*n
                            prcnt=prcnt+np
                            call iosys ('write real '//filnme//' to '//
     1                                  'atomci without rewinding',
     2                                   np*n,rhor,0,' ')
                            cnt=0
                        endif
   20                continue
   10             continue
              else
                  nprtot=na*nb
                  chkwd=nprtot*n
                  do 40 i=1,na
                     do 50 j=1,nb
                        cnt=cnt+1
                        do 60 k=1,n
                           rhor(k,cnt)=fnar(k,i)*fnbr(k,j)
   60                   continue
                        if (cnt.eq.np) then
                            cntwd=cntwd+np*n
                            prcnt=prcnt+np
                            call iosys ('write real '//filnme//' to '//
     1                                  'atomci without rewinding',
     2                                   np*n,rhor,0,' ')
                            cnt=0
                        endif
   50                continue
   40             continue
              endif
              if (cnt.ne.0) then
                  cntwd=cntwd+cnt*n
                  prcnt=prcnt+cnt
                  call iosys ('write real '//filnme//' to atomci '//
     1                        'without rewinding',cnt*n,rhor,0,' ')
              endif                       
              if (cntwd.ne.chkwd) then
                  write(iout,1) cntwd, chkwd
                  call lnkerr('error in write of density file')
              endif
c**********************************************************************c
c                   bound-other                                        c
c**********************************************************************c
          elseif (typb.eq.'complex'.or.typb.eq.'free0'.or.
     1                                 typb.eq.'free1') then
              nprtot=na*nb               
              chkwd=na*nb*2*n
              do 70 i=1,na
                 do 80 j=1,nb
                    cnt=cnt+1
                    do 90 k=1,n
                       rhoc(k,cnt)=fnar(k,i)*fnbc(k,j)
   90               continue
                    if (cnt.eq.np) then
                        prcnt=prcnt+np
                        cntwd=cntwd+2*n*np   
                        call iosys ('write real '//filnme//' to '//
     1                              'atomci without rewinding',
     1                               2*np*n,rhor,0,' ')
                        cnt=0
                    endif
   80            continue
   70         continue
              if (cnt.ne.0) then
                  prcnt=prcnt+cnt
                  cntwd=cntwd+cnt*2*n
                  call iosys ('write real '//filnme//' to atomci '//
     1                        'without rewinding',2*cnt*n,rhor,0,' ')
              endif                            
              if (cntwd.ne.chkwd) then
                  write(iout,1) cntwd, chkwd
                  call lnkerr('error in write of density file')
              endif
          endif
      elseif (typa.eq.'complex') then
c**********************************************************************c
c                    complex-complex                                   c
c**********************************************************************c
          if (typb.eq.'complex') then
              if (type.eq.'diagonal') then
                  nprtot=na*(na+1)/2               
                  chkwd=nprtot*n*2
                  do 100 i=1,na
                     do 200 j=1,i
                        cnt=cnt+1
                        do 300 k=1,n
                           rhoc(k,cnt)=fnac(k,i)*fnac(k,j)
  300                   continue
                        if (cnt.eq.np) then
                            prcnt=prcnt+np
                            cntwd=cntwd+2*np*n
                            call iosys ('write real '//filnme//' to '//
     1                                  'atomci without rewinding',
     2                                   2*np*n,rhor,0,' ')
                            cnt=0
                        endif
  200                continue
  100             continue
              else
                  nprtot=na*nb               
                  chkwd=nprtot*n*2
                  do 400 i=1,na
                     do 500 j=1,nb
                        cnt=cnt+1
                        do 600 k=1,n
                           rhoc(k,cnt)=fnac(k,i)*fnbc(k,j)
  600                   continue
                        if (cnt.eq.np) then
                            prcnt=prcnt+np
                            cntwd=cntwd+2*np*n
                            call iosys ('write real '//filnme//' to '//
     1                                  'atomci without rewinding',
     2                                   2*np*n,rhor,0,' ')
                            cnt=0
                        endif
  500                continue
  400             continue
              endif
              if (cnt.ne.0) then
                  prcnt=prcnt+cnt
                  cntwd=cntwd+cnt*2*n
                  call iosys ('write real '//filnme//' to atomci '//
     1                        'without rewinding',2*cnt*n,rhor,0,' ')
              endif                            
              if (cntwd.ne.chkwd) then
                  write(iout,1) cntwd, chkwd
                  call lnkerr('error in write of density file')
              endif
c**********************************************************************c
c                    complex-free                                      c
c**********************************************************************c
          elseif (typb.eq.'free0'.or.typb.eq.'free1') then
              nprtot=na*nb
              chkwd=na*nb*2*n
              do 700 i=1,na
                 do 800 j=1,nb
                    cnt=cnt+1
                    do 900 k=1,n
                       rhoc(k,cnt)=fnac(k,i)*fnbc(k,j)
  900               continue
                    if (cnt.eq.np) then
                        prcnt=prcnt+np
                        cntwd=cntwd+2*np*n
                        call iosys ('write real '//filnme//' to '//
     1                              'atomci without rewinding',
     2                               2*np*n,rhor,0,' ')
                        cnt=0
                    endif
  800            continue
  700         continue
              if (cnt.ne.0) then
                  prcnt=prcnt+cnt
                  cntwd=cntwd+2*cnt*n
                  call iosys ('write real '//filnme//' to atomci '//
     1                        'without rewinding',2*cnt*n,rhor,0,' ')
              endif                            
              if (cntwd.ne.chkwd) then
                  write(iout,1) cntwd, chkwd
                  call lnkerr('error in write of density file')
              endif
          endif
      elseif (typa.eq.'free0'.or.typa.eq.'free1') then
              if (type.eq.'diagonal') then
                  nprtot=na*(na+1)/2
                  chkwd=nprtot*n*2
                  do 1000 i=1,na
                     do 2000 j=1,i
                        cnt=cnt+1
                        do 3000 k=1,n
                           rhoc(k,cnt)=fnac(k,i)*fnac(k,j)
 3000                   continue
                       if (cnt.eq.np) then
                           prcnt=prcnt+np
                           cntwd=cntwd+2*np*n
                           call iosys ('write real '//filnme//' to '//
     1                                 'atomci without rewinding',
     2                                  2*np*n,rhor,0,' ')
                           cnt=0
                       endif
 2000            continue
 1000         continue
              else
                  nprtot=na*nb
                  chkwd=nprtot*n*2
                  do 5000 i=1,na
                     do 6000 j=1,nb
                        cnt=cnt+1
                        do 7000 k=1,n
                           rhoc(k,cnt)=fnac(k,i)*fnbc(k,j)
 7000                   continue
                        if (cnt.eq.np) then
                            prcnt=prcnt+np
                            cntwd=cntwd+2*np*n
                            call iosys ('write real '//filnme//' to '//
     1                                  'atomci without rewinding',
     2                                   2*np*n,rhor,0,' ')
                            cnt=0
                        endif
 6000                continue
 5000             continue
              endif   
              if (cnt.ne.0) then
                  prcnt=prcnt+cnt
                  cntwd=cntwd+2*cnt*n
                  call iosys ('write real '//filnme//' to atomci '//
     1                        'without rewinding',2*cnt*n,rhor,0,' ')
              endif                            
              if (cntwd.ne.chkwd) then
                  write(iout,1) cntwd, chkwd
                  call lnkerr('error in write of density file')
              endif
      else    
              call lnkerr('wrong type passed into mkden')
      endif
      prtot=prtot+prcnt
      if (prcnt.ne.nprtot) then
          write(iout,2) prcnt, nprtot 
      endif
    1 format(//,5x,'##### error in word count #####',/,5x,
     1             'words written',1x,i6,1x,'should be',1x,i6)    
    2 format(//,'*no. pairs written',1x,i5,1x,'should be',1x,i5)   
      return
      end


