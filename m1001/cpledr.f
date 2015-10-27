*deck @(#)cpledr.f	1.4  8/19/91
      subroutine cpledr(nsym,nbf,nob,nfob,ncob,naob,cv,
     1      nmix,ncsf,mdim,
     4      noci,istopj,stopj,
     5      b, sqcdf,
     6      cr,icr,ncore,lbufso,nograd,ndf,locsym,locp,len,mix,
     7      prtena)
cc
c revision  880322  (yymmdd)
c 22 march 1988  bhl at llnl
c  read ci vector from "mc root" instead of "ci root"
cc
c
      implicit real*8(a-h,o-z)
      integer wpadti
      character*4 itoc
      character*8 mcscf,citype,prtflg
      logical keyval,prtena
      logical debug
cc
      dimension nbf(2),nob(2),ncob(2),naob(2),cv(2)
      dimension cr(2),icr(2),nfob(2)
      dimension b(2)
      integer locp(2),len(2),mix(2),locsym(2)
      data kount/0/
      data debug/.false./
c
      common /io/ inp,iout
c
      common /number/zero,pt5,one,two,four,eight
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
c
      call iosys('read character mcscf_type from mcscr',-1,0,0,mcscf)
c
c..bhl..scf
c
      if(noci.eq.0) then
         call iosys('read character "mcscf: ci used" from rwf',
     $              -1,0,0,citype)
         if(kount.eq.0) then
            write(iout,12345) citype
12345       format(' mcscf type = ',a8)
            kount=1
         end if
      endif
c
c      if (.not.keyval(ops,'nroots',nroots))     nroots=1
c
      nroots=1
c     if (mcroot.ne.0) nroots=mcroot
c
c----------------------------------------
c     storage allocation for jacobi
c----------------------------------------
c
      nroots=1
      nao=naob(1)
      nnao=nao*(nao+1)/2
      ntnao=nnao*nnao
      npassj=1
      lenb=lbufso
c
c     allocate core
      call getscm(0,cr,maxcor,'size',0)
      maxcor=iadtwp(maxcor-100000)
c
      ibufix=1
      lgrad=ibufix+lbufso
      ldiag = lgrad + mdim
      lt = ldiag + mdim
      lc = lt + mdim*(istopj+1)
      lc0 = lc + mdim
      lhess = lc0 + ncsf
c
      nhess = istopj * (istopj + 1) / 2
c
      lthess = lhess + nhess
      lgl = lthess + nhess
      lolap = lgl + istopj
      lsg = lolap + istopj
c
      naa=naob(1)*(naob(1)+1)/2
      nai=naob(1)*nob(1)
      nci=ncob(1)*nob(1)
      naa2=naa*naob(1)*naob(1)
c
      ltg=lsg+naa2
      labcx=ltg+naa2
      labix=labcx+naa*nai
      laibx=labix+naa*nci
      lstor=laibx+1
c
      incor=1
      ihss=lstor
      if(lstor.gt.maxcor) then
         lstor=labcx
         labix=lstor
         laibx=lstor
         incor=0
      else
         lhss=nmix*nmix
         if(lstor+lhss+lenb.lt.maxcor) then
            ihss=lstor
            ibf=ihss+lhss
            lstor=ibf+lenb
            incor=2
         end if
       end if
c
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(wpadti(lstor),cr,ncore,'mc linear equations',0)
c
      ncor = iadtwp(ncore)- lstor
c
      if(incor.eq.2) then
         call gethes(cr(ihss),lhss,cr(ibf),lenb)
      end if
c
      if(incor.ne.0) then
          call getrab(cr(labcx),cr(ibufix),lbufso,'abcx',
     #                naob(1),nob(1),naob(1))
          if(ncob(1).ne.0) then
             call getrab(cr(labix),cr(ibufix),lbufso,'abix',
     #                   ncob(1),nob(1),naob(1))
          endif
      endif
c
c
c--------------------------------
c     read ci vector
c--------------------------------
      if(noci.eq.0) then
         call iosys('read real "mc root '//itoc(nroots)//'" from rwf',
     #               ncsf,cr(lc0),0,' ')
c
c---------------------------------
c     calculate ci gradient
c---------------------------------
c
      zero=0.d0
      junk=0
      if (citype.eq.'m902') then
         call mn820(cr(lstor),cr(lstor),ncor,'mcscr','mcscr',
     $              'mcscf','ci')
         call mn902(cr(lstor),cr(lstor),ncor,'h times c',cr(lc0),
     $        cr(lgrad+nmix),'mcscr','mcscr','ci',junk,'mc')
      else if (citype.eq.'m903') then
         call mn903(cr(lstor),cr(lstor),ncor,'h times c',
     $        cr(lc0),cr(lgrad+nmix),'mcscr','mcscr','ci',junk,'mc')
      end if
cccccc
c
c
c-------------------------
c     calculate ci energy
c-------------------------
      energy = 0.d0
      lc00 = lc0 - 1
      lg0 = lgrad + nmix - 1
      do 310 i = 1, ncsf
  310 energy = energy + cr(lc00+i) * cr(lg0+i)
c

      if(prtena.or.debug) then
         write (iout,3001) energy
 3001    format(/,' cpledr: energy = ',1x,f24.16)
      end if
c
c-----------------------------------
c     read ci diagonal elements
c-----------------------------------
c
      ld0 = ldiag + nmix - 1
      call iosys('read real "h diagonals" from rwf',ncsf,
     #           cr(ld0+1),0,' ')
c
      endif
c
ccc
c----------------------------------------------
c     modify ci gradient and diagonal elements
c----------------------------------------------
cc
      do 510 i = 1, ncsf
         cr(ld0+i) = cr(ld0+i) - energy
  510 continue
c
      if(nograd.eq.0) then
         do 511 i = 1, ncsf
         cr(lg0+i) = cr(lg0+i) - energy * cr(lc00+i)
  511    continue
      end if
c
cps      write(6,91929)(cr(lg0+i),i=1,ncsf)
cps91929 format(/,' ci gradient ',6(1x,f15.10))
c
c---------------------------------
c     initialize mcmxvc
c---------------------------------
      lbuf = lstor + 1
      lstor2 = lbuf + lenb
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(wpadti(lstor2),cr,ncore,'mc linear equations 2',0)
c
      ncor2 = iadtwp(ncore) - lstor2
c
      naobt = 0
      do 610 i = 1, nsym
  610 naobt = naobt + naob(i)
      lij = naobt * (naobt + 1) / 2
      lijkl=(lij*(lij+1))/2
c
c     call enmxvc(nsym,nbf,nob,nfob,ncob,naob,cv,cr(lc0),
c    1      locsym,len,locp,mix,
c    2      lij,lijkl,
c    3      nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
c    4      nfsm,ndshd1,ndshd2,nhd,
c    5      cr(lbuf),cr(lstor2),icr(lstor2),
c    6      ncor2)
c
c-------------------------------------------------------------
c     read diagonal elements of hessian and orbital gradient
c-------------------------------------------------------------
cc

      call rdiagh(cr(ldiag),cr(lbuf),lenb,nmix)
cc
c
      if(nograd.eq.0) then
         call iosys('rewind mcscf_gradient on rwf',0,0,0,' ')
         call iosys('read real mcscf_gradient from rwf'
     $              //' without rewinding',nmix,cr(lgrad),0,' ')
      else
         lg0=lgrad-1
         do 821 k=1,mdim
            cr(lg0+k)=b(k)
  821    continue
      end if
c
c
c
c---------------------------------
c     solve linear equations
c---------------------------------
c
      call mcjcbi(cr(lgrad),cr(ldiag),cr(lt),b,cr(lc0),
     $     cr(lgl),cr(lolap),cr(lhess),cr(lthess),cr(lc),
     $     nmix,ncsf,mdim,energy,
     $     ndstv,nf45,nf46,lb46,noci,npassj,istopj,stopj,sqcdf,
     $     cr(ibufix),lbufso,cr(labcx),cr(labix),cr(laibx),
     $     cr(ihss),incor,cr(lsg),cr(ltg),
     $     nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,locp,mix,
     $     lij,lijkl,
     $     nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
     $     nfsm,ndshd1,ndshd2,nhd,
     $     cr(lbuf),cr(lstor2),icr(wpadti(lstor2)),
     $     ncor2)

c
c     write (6,603) (b(i),i=1,mdim)
c 603 format('0*mcledr solution '//4(1x,f16.8))
c
c
      return
      end
