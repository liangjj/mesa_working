c $Header$
*deck genleg.f
c***begin prologue     genleg
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legndre, link 6203
c***author             schneider, barry (lanl)
c***source             m6203
c***purpose            generate real spherical harmonics from m = 0 to
c***                   mmax and l = m to lmax for each shell.
c***                   both m values and all l values are stored.      
c***references         
c***routines called
c***end prologue       genleg
      subroutine genleg (plm,phim,cphi,sphi,phi,wphi,theta,wthet,
     1                   dfct,ddfct,mpnt,lpnt,ns,nth,nph,lmax,
     2                   mmax,maxfac,buffer,dirin,dirout)
      implicit integer(a-z)
      dimension plm(buffer), phim(*), cphi(*), sphi(*), phi(*)
      dimension theta(*), wphi(*), wthet(*), mpnt(ns), lpnt(0:mmax,ns)
      dimension nth(ns), nph(ns)
      real*8 plm, phim, cphi, sphi, phi, theta, dfct, ddfct
      real*8 wphi, wthet
      character*2 itoc
      character*(*) dirin, dirout
      common /io/ inp, iout
      phiwdc=1
      thewdc=1
      do 10 is=1,ns
         call iosys ('read real "phi points shell-'//itoc(is)//
     1               '" from lamdat',nph(is),phi(phiwdc),0,' ')
         call iosys ('read real "phi weights shell-'//itoc(is)//
     1               '" from lamdat',nph(is),wphi(phiwdc),0,' ')
         phiwdc=phiwdc+nph(is) 
         thewdc=thewdc+nth(is)
   10 continue
      if (dirin.eq.'create') then
          words=0
          words1=0
          do 20 is=1,ns
             words=words+(2*mmax+1)*nph(is)
             do 30 m=0,mmax
                words1=words1+(lmax-m+1)*nth(is)
   30        continue
   20     continue
          call iosys ('create real phim on lamdat',words,0,0,' ')
          call iosys ('create real plm on lamdat',words1,0,0,' ')
c
c         mpnt is a pointer to the first word for m = 0 in each shell.
c         lpnt is a pointer to the first l value for a given m in each shell
          plmcnt=0
          phicnt=0
          totplm=0
          phiwdc=1
          thewdc=1       
          do 40 is=1,ns
             mpnt(is)=phicnt
             wd1=nph(is)+nph(is)
             wd2=nph(is) 
             call miscfn(phi(phiwdc),cphi,sphi,nph(is))
             do 50 m=0,mmax
                wds=(lmax+1-m)*nth(is)
                wds1=wd1
                if (m.eq.0) then
                    wds1=wd2
                endif  
                lpnt(m,is)=totplm+1
                call scmphi(cphi,sphi,phim(phicnt+1),nph(is),m)
                testwd=plmcnt+wds
                if (testwd.gt.buffer) then
                    call iosys ('write real plm on lamdat without '//
     1                          'rewinding',plmcnt,plm,0,' ')
                    totplm=totplm+plmcnt
                    plmcnt=0
                endif
                call legend(plm(plmcnt+1),theta(thewdc),dfct,ddfct,
     1                      nth(is),lmax,m,maxfac)
                plmcnt=plmcnt+wds
                phicnt=phicnt+wds1
   50        continue
             phiwdc=phiwdc+nph(is) 
             thewdc=thewdc+nth(is)
   40     continue
          if (phicnt.ne.words) then
              call lnkerr('error in words written to phim file')         
          endif
          call iosys ('write real phim on lamdat',phicnt,phim,0,' ')
          if(plmcnt.ne.0) then
             call iosys ('write real plm on lamdat without '//
     1                   'rewinding',plmcnt,plm,0,' ')
             totplm=totplm+plmcnt
          endif
          if (totplm.ne.words1) then
              call lnkerr('error in words written to legendre file')
          endif
          call iosys ('write integer "number of phim words" to lamdat',
     1                 1,phicnt,0,' ')
          call iosys ('write integer "number of plm words" to lamdat',
     1                 1,words1,0,' ')
          call iosys ('write integer "legendre m pointer" to lamdat',
     1                 ns,mpnt,0,' ')
          call iosys ('write integer "legendre l pointer" to lamdat',
     1                 (mmax+1)*ns,lpnt,0,' ')
          call iosys ('rewind phim on lamdat read-and-write',0,0,0,' ')
          call iosys ('rewind plm on lamdat read-and-write',0,0,0,' ')
          dirout='on-disk'
          if (words1.le.buffer) then
              dirout='in-core'
          endif
      elseif (dirin.eq.'read-all') then
          call iosys ('read integer "number of phim words" from lamdat',
     1                 1,phicnt,0,' ')
          call iosys ('read integer "number of plm words" from lamdat',
     1                 1,plmcnt,0,' ')
          call iosys ('read integer "legendre m pointer" from lamdat',
     1                 ns,mpnt,0,' ')
          call iosys ('read integer "legendre l pointer" from lamdat',
     1                 (mmax+1)*ns,lpnt,0,' ')
          call iosys ('read real phim from lamdat',phicnt,phim,0,' ')
          call iosys ('read real plm from lamdat',plmcnt,plm,0,' ')
          call iosys ('rewind phim on lamdat read-and-write',0,0,0,' ')
          call iosys ('rewind plm on lamdat read-and-write',0,0,0,' ')
          dirout='in-core'
      else
          dirout='on-disk'
      endif
      return
      end
