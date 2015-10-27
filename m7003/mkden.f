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
     1                 na,nb,np,n)
      implicit integer (a-z)
      real *8 fnar, fnbr, rhor
      complex *16 fnac, fnbc, rhoc
      character *(*) typa, typb
      character *128 filnme
      dimension fnar(n,na), fnbr(n,nb), fnac(n,na), fnbc(n,nb)
      dimension rhor(n,np), rhoc(n,np)
      common /io/ inp, iout
      filnme='"'//typa//'-'//typb//' densities"'
      if (typa.eq.'bound') then
          if (typb.eq.'bound') then
              nprtot=na*(na+1)/2
              chkwd=na*(na+1)*n/2
              cntwd=0
              cnt=0
              prcnt=0
              do 10 i=1,na
                 do 20 j=1,i
                    cnt=cnt+1
                    do 30 k=1,n
                       rhor(k,cnt)=fnar(k,i)*fnar(k,j)
   30               continue
                    if (cnt.eq.np) then
                        cntwd=cntwd+np*n
                        prcnt=prcnt+np
                        call iosys ('write real '//filnme//' to '//
     1                              'atomci without rewinding',
     2                               np*n,rhor,0,' ')
                        cnt=0
                    endif
   20            continue
   10         continue
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
          elseif (typb.eq.'complex'.or.typb.eq.'free0'.or.
     1                                 typb.eq.'free1') then
              nprtot=na*nb               
              chkwd=na*nb*2*n
              cnt=0
              cntwd=0
              prcnt=0
              do 40 i=1,na
                 do 50 j=1,nb
                    cnt=cnt+1
                    do 60 k=1,n
                       rhoc(k,cnt)=fnar(k,i)*fnbc(k,j)
   60               continue
                    if (cnt.eq.np) then
                        prcnt=prcnt+np
                        cntwd=cntwd+2*n*np   
                        call iosys ('write real '//filnme//' to '//
     1                              'atomci without rewinding',
     1                               2*np*n,rhoc,0,' ')
                        cnt=0
                    endif
   50            continue
   40         continue
              if (cnt.ne.0) then
                  prcnt=prcnt+cnt
                  cntwd=cntwd+cnt*2*n
                  call iosys ('write real '//filnme//' to atomci '//
     1                        'without rewinding',2*cnt*n,rhoc,0,' ')
              endif                            
              if (cntwd.ne.chkwd) then
                  write(iout,1) cntwd, chkwd
                  call lnkerr('error in write of density file')
              endif
          endif
      elseif (typa.eq.'complex') then
          if (typb.eq.'complex') then
              nprtot=na*(na+1)/2               
              chkwd=na*(na+1)*n
              cntwd=0 
              cnt=0
              prcnt=0
              do 70 i=1,na
                 do 80 j=1,i
                    cnt=cnt+1
                    do 90 k=1,n
                       rhoc(k,cnt)=fnac(k,i)*fnac(k,j)
   90               continue
                    if (cnt.eq.np) then
                        prcnt=prcnt+np
                        cntwd=cntwd+2*np*n
                        call iosys ('write real '//filnme//' to '//
     1                              'atomci without rewinding',
     2                               2*np*n,rhoc,0,' ')
                        cnt=0
                    endif
   80            continue
   70         continue
              if (cnt.ne.0) then
                  prcnt=prcnt+cnt
                  cntwd=cntwd+cnt*2*n
                  call iosys ('write real '//filnme//' to atomci '//
     1                        'without rewinding',2*cnt*n,rhoc,0,' ')
              endif                            
              if (cntwd.ne.chkwd) then
                  write(iout,1) cntwd, chkwd
                  call lnkerr('error in write of density file')
              endif
          elseif (typb.eq.'free0'.or.typb.eq.'free1') then
              nprtot=na*nb
              chkwd=na*nb*2*n
              cntwd=0
              cnt=0
              prcnt=0
              do 100 i=1,na
                 do 200 j=1,nb
                    cnt=cnt+1
                    do 300 k=1,n
                       rhoc(k,cnt)=fnac(k,i)*fnbc(k,j)
  300               continue
                    if (cnt.eq.np) then
                        prcnt=prcnt+np
                        cntwd=cntwd+2*np*n
                        call iosys ('write real '//filnme//' to '//
     1                              'atomci without rewinding',
     2                               2*np*n,rhoc,0,' ')
                        cnt=0
                    endif
  200            continue
  100         continue
              if (cnt.ne.0) then
                  prcnt=prcnt+cnt
                  cntwd=cntwd+2*cnt*n
                  call iosys ('write real '//filnme//' to atomci '//
     1                        'without rewinding',2*cnt*n,rhoc,0,' ')
              endif                            
              if (cntwd.ne.chkwd) then
                  write(iout,1) cntwd, chkwd
                  call lnkerr('error in write of density file')
              endif
          endif
      elseif (typa.eq.'free0'.or.typa.eq.'free1') then
              nprtot=na*(na+1)/2
              chkwd=na*(na+1)*n
              cntwd=0
              cnt=0
              prcnt=0
              do 400 i=1,na
                 do 500 j=1,i
                    cnt=cnt+1
                    do 600 k=1,n
                       rhoc(k,cnt)=fnac(k,i)*fnac(k,j)
  600               continue
                    if (cnt.eq.np) then
                        prcnt=prcnt+np
                        cntwd=cntwd+2*np*n
                        call iosys ('write real '//filnme//' to '//
     1                              'atomci without rewinding',
     2                               2*np*n,rhoc,0,' ')
                        cnt=0
                    endif
  500             continue
  400         continue
              if (cnt.ne.0) then
                  prcnt=prcnt+cnt
                  cntwd=cntwd+2*cnt*n
                  call iosys ('write real '//filnme//' to atomci '//
     1                        'without rewinding',2*cnt*n,rhoc,0,' ')
              endif                            
              if (cntwd.ne.chkwd) then
                  write(iout,1) cntwd, chkwd
                  call lnkerr('error in write of density file')
              endif
      else    
              call lnkerr('wrong type paseed into mkden')
      endif
      write(iout,2) prcnt, nprtot 
    1 format(//,5x,'##### error in word count #####',/,5x,
     1             'words written',1x,i6,1x,'should be',1x,i6)    
    2 format(//,'*no. pairs written',1x,i5,1x,'should be',1x,i5)   
      return
      end
