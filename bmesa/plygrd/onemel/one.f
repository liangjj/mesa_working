*deck one.f 
c***begin prologue     one
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           kinetic, potential, matrix elements
c***author             schneider, b. i.(nsf)
c***source             one
c***purpose            kinetic energy plus one body potential
c**                    matrix elements in dvr basis.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       one
      program one
c
      implicit integer (a-z)
      parameter ( maxgrd=20 ) 
      character*4096 ops
      character*2 itoc, ic, jc
      character*80 cpass, title, chrkey, prtit
      character*24 coord, comp, typfn, timh0
      character*32 qdtyp, titphr
      character*800 card
      character*128 filbec
      character*8 qtyp, phrse, atom
      character*3 notim, nospac
      logical dollar, logkey, toau, totrap, useau, prnhmo
      real*8 z, y
      real*8 fpkey
      real*8 pi, hbar, massau, lenau, timau, amass, omegat, bigomg
      real*8 lenscl, escal
      dimension npt(maxgrd), nmax(maxgrd)
      dimension q(maxgrd), wt(maxgrd) 
      dimension p(maxgrd,maxgrd), dp(maxgrd,maxgrd), ddp(maxgrd,maxgrd)
      dimension omegat(3)
      common/io/inp, iout      
      pointer (p1,z(1)), (p1,iz(1))
      pointer (p2,y(1)), (p2,iy(1))
      data pi/3.1415926535897932384d0/
      data massau, lenau, timau / 9.109558d-31, 5.291771d-11, 
     1                            2.41888d-17 /
c     hbar in joule-sec      
      data hbar/1.054592d-34/
c
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      call iosys ('open bec as old',0,0,0,filbec)
      call iosys('read integer "no. space dimensions" from bec',1,
     1            nosdim,0,' ')
      call iosys('read character "space dimension" from bec',
     1            0,0,0,nospac)
      call iosys('read character "time dimension" from bec',0,0,0,notim)
c
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      totrap=logkey(ops,'use-trap-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
      prnhmo=logkey(ops,'print=m6291=hamiltonian',.false.,' ')
      timh0=chrkey(ops,'pure-time-perturbation','none',' ')
      nospac=chrkey(ops,'m6291=space-dimension','yes',' ')
      notim=chrkey(ops,'m6291=time-dimension','no',' ')
      write(iout,1)
      if(toau) then
         hbar=1.d0
      endif
      if(useau) then
         hbar=1.d0
      endif
      amass=1.d0
      atom=chrkey(ops,'atom','cs',' ')
      if( dollar('$trap',card,cpass,inp) ) then
          call atmdat(atom,amass,omegat)
          if(toau) then
              write(iout,*) '          converting to atomic units'
              omegat(1)=omegat(1)*timau
              omegat(2)=omegat(2)*timau
              omegat(3)=omegat(3)*timau
              amass=amass/massau
          endif
          if(useau) then
             write(iout,*) '          assuming atomic units'
             amass=1.d0
             omegat(1)=1.d0
             omegat(2)=1.d0
             omegat(3)=1.d0
          endif
          if(totrap) then
             bigomg=max(omegat(1),omegat(2),omegat(3))
             omegat(1)=omegat(1)/bigomg
             omegat(2)=omegat(2)/bigomg
             omegat(3)=omegat(3)/bigomg
             lenscl=sqrt(hbar/(amass*bigomg))
             escal=hbar*bigomg
          endif
      endif
      if(nospac.eq.'yes') then
         write(iout,*)  'processing space information'
c
         do 10 i=1,nosdim
            ic=itoc(i)
            qtyp=chrkey(ops,'space-dimension-'//ic,'x',' ')
            phrse=qtyp
            call iosys('read integer "no. subgrids for '
     1                 //phrse//'" from bec',1,nsubg,0,' ')
            call iosys('read integer "left fn value for '
     1                 //phrse//'" from bec',1,fl,0,' ')
            call iosys('read integer "right fn value for '
     1                 //phrse//'" from bec',1,fr,0,' ')
            call iosys('read integer "no. points for '
     1                 //phrse//'" from bec',nsubg,npt,0,' ')
            call iosys('read integer "mod. no. points '//
     1                 'for '//phrse//'" from bec',nsubg,nmax,
     2                  0,' ')
            titphr=qtyp//' pointer'
            call iosys('read character "'//titphr//'" from bec',
     1                  0,0,0,title)
            call iosys('read integer '//title//' from bec '//
     1                 'without rewinding',nsubg,q,2*nsubg,' ')
            call iosys('read integer '//title//' from bec without '//
     1                 'rewinding',nsubg,wt,0,' ')
            call iosys('read integer '//title//' from bec without '//
     1                 'rewinding',maxgrd*maxgrd,p,3*maxgrd*maxgrd,' ')
            call iosys('read integer '//title//' from bec without '//
     1                 'rewinding',maxgrd*maxgrd,dp,0,' ')
            call iosys('read integer '//title//' from bec without '//
     1                 'rewinding',maxgrd*maxgrd,ddp,0,' ')
            titphr=qtyp//' title'
            call iosys('read character "'//titphr//'" from bec',
     1                  0,0,0,title)
            do 20 j=1,nsubg
               write(iout,2) qtyp, j, npt(j), nmax(j)
               jc=itoc(j)
               zq=1
               zwt=zq+npt(j)
               zp=zwt+npt(j)
               zdp=zp+npt(j)*npt(j)
               zddp=zdp+npt(j)*npt(j)
               need=wpadti(zddp+npt(j)*npt(j))
               call memory(need,p1,ngot1,'functions',0)
               call iosys('read real '//title//' from bec',npt(j),
     1                     z(zq),q(j)-1,' ')
c               prtit='points'
c               call prntrm(prtit,z(zq),nmax(j),1,nmax(j),1,iout)
                call iosys('read real '//title//' from bec',npt(j),
     1                     z(zwt),wt(j)-1,' ')
c               prtit='weights'
c               call prntrm(prtit,z(zwt),nmax(j),1,nmax(j),1,iout)
                call iosys('read real '//title//' from bec',
     1                     npt(j)*npt(j),z(zp),p(j,j)-1,' ')
                call iosys('read real '//title//' from bec',
     1                     npt(j)*npt(j),z(zdp),dp(j,j)-1,' ')
                call iosys('read real '//title//' from bec',
     1                     npt(j)*npt(j),z(zddp),ddp(j,j)-1,' ')
c               prtit='p'
c               call prntrm(prtit,z(zp),nmax(j),nmax(j),
c     1                     nmax(j),nmax(j),iout)
c               prtit='dp'
c               call prntrm(prtit,z(zdp),nmax(j),nmax(j),
c     1                     nmax(j),nmax(j),iout)
c               prtit='ddp'
c               call prntrm(prtit,z(zddp),nmax(j),nmax(j),
c     1                     nmax(j),nmax(j),iout)
               v=1
               ham=v+nmax(j)
               grad=ham+nmax(j)*nmax(j)
               need=wpadti(grad+nmax(j)*nmax(j))
               call memory(need,p2,ngot2,'ham',0)
               call vmat1d(z(zq),y(v),hbar,amass,lenscl,escal,
     1                     nmax(j),phrse,j,prnhmo)
               call iosys('write real "v mtrx-'//jc//' for '//
     1                     qtyp//'" to bec',nmax(j),y(v),0,' ')
               call h0(z(zp),z(zdp),z(zddp),z(zq),z(zwt),
     1                 y(ham),y(v),y(grad),hbar,amass,nmax(j),fl,fr,
     2                 qtyp,jc,.false.)
               call memory(-ngot2,p2,idum,'ham',idum)
               call memory(-ngot1,p1,idum,'functions',idum)
 20         continue   
 10      continue
      endif   
      if(notim.eq.'yes') then
         call iosys('read integer "no. time regions" from bec',1,
     1               ntreg,0,' ')
         write(iout,*) 'processing time information'
         do 30 i=1,ntreg
            ic=itoc(i)
            phrse='t-'//ic
            call iosys('read integer "no. subgrids for '//
     1                  phrse//'" from bec',1,nsubg,0,' ')
            call iosys('read integer "no. points for '//
     1                  phrse//'" from bec',nsubg,npt,0,' ')
            call iosys('read integer "mod. no. points for '
     1                      //phrse//'" from bec',nsubg,nmax,0,' ')
            titphr='t pointer region '//ic
            call iosys('read character "'//titphr//'" from bec',
     1                  0,0,0,title)
            call iosys('read integer '//title//' from bec without '//
     1                 'rewinding',nsubg,q,2*nsubg,' ')
            call iosys('read integer '//title//' from bec without '//
     1                 'rewinding',nsubg,wt,0,' ')
            call iosys('read integer '//title//' from bec without '//
     1                 'rewinding',maxgrd*maxgrd,p,3*maxgrd*maxgrd,' ')
            call iosys('read integer '//title//' from bec without '//
     1                 'rewinding',maxgrd*maxgrd,dp,0,' ')
            call iosys('read integer '//title//' from bec without '//
     1                 'rewinding',maxgrd*maxgrd,ddp,0,' ')
            titphr='t title region '//ic
            call iosys('read character "'//titphr//'" from bec',
     1                  0,0,0,title)
            do 40 j=1,nsubg
               write(iout,2) phrse, j, npt(j), nmax(j)
               jc=itoc(j)
               zq=1
               zwt=zq+npt(j)
               zp=zwt+npt(j)
               zdp=zp+npt(j)*npt(j)
               zddp=zdp+npt(j)*npt(j)
               need=wpadti(zddp+npt(j)*npt(j))
               call memory(need,p1,ngot1,'time',0)
               call iosys('read real '//title//' from bec',npt(j),
     1                     z(zq),q(j)-1,' ')
               call iosys('read real '//title//' from bec',npt(j),
     1                     z(zwt),wt(j)-1,' ')
               call iosys('read real '//title//' from bec',
     1                     npt(j)*npt(j),z(zp),p(j,j)-1,' ')
               call iosys('read real '//title//' from bec',
     1                    npt(j)*npt(j),z(zdp),dp(j,j)-1,' ')
               call iosys('read real '//title//' from bec',
     1                     npt(j)*npt(j),z(zddp),ddp(j,j)-1,' ')
c               prtit='p'
c               call prntrm(prtit,z(zp),nmax(j),nmax(j),
c     1                                 nmax(j),nmax(j),iout)  
c               prtit='dp'
c               call prntrm(prtit,z(zdp),nmax(j),nmax(j),
c     1                                 nmax(j),nmax(j),iout)  
c               prtit='ddp'
c               call prntrm(prtit,z(zddp),nmax(j),nmax(j),
c     1                                 nmax(j),nmax(j),iout)  
               v=1
               ham=v+nmax(j)
               need=wpadti(ham+2*nmax(j)*nmax(j))
               call memory(need,p2,ngot2,'hamt',0)
               call vt0(z(zq),y(v),nmax(j),timh0,phrse,j)      
               call iosys('write real "vt mtrx-'//jc//' for '//
     1                     phrse//'" to bec',nmax(j),y(v),0,' ')
               call h1(z(zp),z(zdp),z(zq),z(zwt),y(ham),hbar,
     1                 phrse,jc,nmax(j),prnhmo)
               call memory(-ngot1,p1,idum,'time',idum)
               call memory(-ngot2,p2,idum,'hamt',idum)
 40         continue   
 30      continue   
      endif
      call iosys('rewind all on bec read-and-write',0,0,0,' ')
      call chainx(0)               
      stop
 1    format(/,20x,'kinetic energy and/or time derivative code',//)
 2    format(/,5x,'computing matrix for coordinate = ',a8,/,5x,
     1            'grid number                     = ',i3,/,5x,
     2            'number of points                = ',i3,/,5x,
     3            'modified number of points       = ',i3)
      end


