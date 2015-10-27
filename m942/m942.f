*deck @(#)m942.f	1.4  8/3/91
      program m942
c
c***begin prologue     m942
c***date written       000811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci, hamiltonian matrix, contracted
c***                   basis
c***author             schneider, barry(nsf)
c***source             @(#)m942.f	1.4   8/3/91
c
c***purpose            to construct matrix elements of the 
c***                   Hamiltonian in the contracted CI basis. 
c***description        raw Hamiltonian matrix elements are
c***                   read in which have been contructed in m940.
c***                   these raw matrix elements are over a 
c***                   primitive CI basis which has been divided
c***                   into a P and Q space.  the Q space part of
c***                   the Hamiltonian always proceeds the P space
c***                   part by construction and the Hamiltonian
c**                    matrix elements are stored on a file
c***                   containing only the non-zero elements with
c***                   indices.  note that there may be more than
c***                   one non-zero entry corresponding to a given
c***                   (I,J) index so that an entire pass of the
c***                   set is required to do the contraction.
c***                   the contracted H_PP and H_QP matrices are
c***                   also stored in compact form with only the
c***                   non-zero elements and their indices on the
c***                   file.  the current implimentation does NOT
c***                   allow for a contraction in Q space.  at a
c***                   later date one could, repartition the space
c***                   and move configurations from P to Q.  I
c***                   have also allowed for the complement of
c***                   the contracted P space to be constructed.
c***                   one could move this back into Q space if
c***                   desired.
c
c***references
c
c***routines called    (none)
c
c***end prologue       m942
c
c
      implicit integer (a-z)
c
      parameter (maxeng=300)
      character*4096 ops
      logical logkey, pcomp, noopt, incore, outcore
      logical nodsk, prnt, prndvd
      character*8 prtflg
      character*8 dsk, scatyp, filtyp
      character*4 itoc
      character*128 namham, nmcnfg, tmp, nmtarg, namchk, nmfile
      character*3 answer
      character*16 chrkey, unit, key
      character*32 header(10)
      character*80 title
      real*8 ci, pvec, hpp, hqp, buf, bufq, cvec, rjunk, opt
      real*8 energy, relene, rep, refeng
      real*8 fzcore
      real*8 fpkey
      real*8 thresh, cnverg
      real*4 time(2), delta(10), secnds
      dimension energy(maxeng), relene(maxeng), prnt(5), prndvd(11)
      pointer (pci,ci(1))
      pointer (ppvec,pvec(1)), (ppvec,apvec(1))
      pointer (pcvec,cvec(1)), (pcvec,acvec(1))
      pointer (pbuf,buf(1)), (pbuf,abuf(1))
      pointer (qbuf,bufq(1)), (qbuf,abufq(1))
      pointer (phpp,hpp(1)), (phpp,ahpp(1))
      pointer (phqp,hqp(1)), (phqp,ahqp(1))
      pointer(pjunk,rjunk(1))
      pointer (phopt,opt(1)), (phopt,aopt(1))
c
      common /io/ inp,iout
      data header /'buffers', 
     #             '"H_PP diag"',
     #             '"H_PP off diag"',
     #             '"number H_PP elements"',
     #             '"H_QP off diag"',
     #             '"number H_QP elements"',
     #             '"optical potential"',
     #             '"H_QQ diag"',
     #             '"H_QQ off diag"',
     #             '"number H_QQ elements"' /
      call drum
      time(1)=secnds(0.0)
      call iosys('read character options from rwf',-1,0,0,ops)
      scatyp=chrkey(ops,'scattering','none',' ')
      if(scatyp(1:4).eq.'none') then
         call chainx(0)
         stop
      endif
      if(scatyp(1:4).eq.'kohn') then
         call iosys('read character "kohn filename" from rwf',
     #               -1,0,0,nmfile)
         filtyp='kohn'
         filtyp=filtyp(1:4)
         len=4
      else if(scatyp(1:8).eq.'r-matrix') then
         call iosys('read character "r-matrix filename" from rwf',
     #               -1,0,0,nmfile)
         filtyp='rmtrx'
         filtyp=filtyp(1:5)
         len=8
      endif
      call iosys('open '//filtyp//' as old',0,0,0,nmfile)
      prnt(1)=logkey(ops,'print=ci=target-energies',.false.,' ')
      prnt(2)=logkey(ops,'print=ci=p-space-vectors',.false.,' ')
      prnt(3)=logkey(ops,'print=ci=hpp',.false.,' ')
      prnt(4)=logkey(ops,'print=ci=hqp',.false.,' ')
      prnt(5)=logkey(ops,'print=ci=hqq',.false.,' ')
      prnt(6)=logkey(ops,'print=ci=optical-potential',.false.,' ')
      prndvd(1)=logkey(ops,'print=ci=davidson',.false.,' ')      
      if(prndvd(1)) then
         do 10 i=2,11
            prndvd(i)=.true.
 10      continue
      endif   
      mcorci='ci'
      citype='m942'
c
      write(iout,*) 'm942: construct contracted CI matrix elements'
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
c     
c
c     open the file holding the target information
c     and read the information needed for the target CI
c     vectors
c
      call getenv('MESA_TMP',tmp)
      ltmp=cskipb(tmp,' ')
      nmtarg=tmp(1:ltmp)//'/target'
      call iosys('open target as old',0,0,0,nmtarg)
      write(iout,*) '          reading target information '
      call iosys('rewind all on target read-and-write',
     $            0,0,0,' ')
      call iosys('read integer "target walks" '//
     1           'from target',1,nwks,0,' ')
      call iosys('read integer "target orbitals" '//
     #           'from target',1,nsmall,0,' ')
      call iosys('read integer "target roots" from target',
     1                    1,roots,0,' ')
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,rep,0,' ')
      call iosys('read real "reference energy" from '//filtyp,
     $            1,refeng,0,' ')
      nen=intkey(ops,'scattering='//scatyp(1:len)
     #                //'=number-of-energies',1,' ')
      call fparr(ops,'scattering='//scatyp(1:len)//'=energy',relene,
     #                nen,' ')
      pcomp=logkey(ops,'scattering='//scatyp(1:len)//'=p-complement',
     #             .false.,' ')  
      thresh=fpkey(ops,'ci=tolerance',1.0d-10,' ')
      cnverg=fpkey(ops,'ci=convergence',1.0d-08,' ')
      mxiter=intkey(ops,'ci=iterations',15,' ')
      mxdvd=intkey(ops,'ci=davidson-vectors',5*npvec,' ')
      do 20 i=1,nen
         energy(i)=relene(i)+refeng
 20   continue   
c
c     open the main configuration file.  its name has been 
c     written to namchk.
c
      call iosys('read character "checkpoint filename" from rwf',
     $            0,0,0,namchk)
      call iosys('open chk as old',0,0,0,namchk)
      key=chrkey(ops,'int=drt=key','drt',' ')
      call pakstr(key,lenkey)
      call iosys('read character "drt file name '//key(1:lenkey)
     #            //'" from chk',0,0,0,dsk)
      call iosys('read character "drt unit file name '//dsk
     #               //'" from chk',0,0,0,nmcnfg)
      call iosys('open '//dsk//' as unknown',0,0,0,nmcnfg)
      write(iout,*) '          reading information from '//dsk
c
c    ----- read the needed information on the drt file
c
      call iosys('rewind all on '//dsk//' read-and-write',0,0,0,' ')
c
c     get the number of P and Q space configurations.
c
      call iosys('read integer nwksp from '//dsk,1,nwksp,0,' ')
      call iosys('read integer nwksq from '//dsk,1,nwksq,0,' ')
c
c     find out the size of the contracted P space.
c
      nconf=nwksp+nwksq
      nl2=nwksp/nwks
      npvec=roots*nl2
      write(iout,1) nsmall, roots, nl2, nwks, nwksp, nwksq, 
     #              nconf, npvec, refeng
      write(iout,2)(relene(i),i=1,nen)
      write(iout,3)(energy(i),i=1,nen)
      call getmem(0,pci,nwci,'first',0)
      call iosys('read integer mxcore from rwf',1,maxcor,0,' ')
c
c     build the contracted p-space
c
      ciene=1
      civec=ciene + roots
      need=wpadti(civec+nwks*roots)
      call getmem(need,pci,nwci,'civec',0)
      maxcor=maxcor-nwci
      vec=1
      need=wpadti(vec+npvec*nconf)
      call getmem(need,ppvec,nwpvec,'pvec',0)
      maxcor=maxcor-nwpvec
      ie=ciene
      iv=civec
      do 30 i=1,roots
         call iosys('read real "ci energy '//itoc(i)//'" from target',
     #              1,ci(ie),0,' ')
         call iosys('read real "ci root '//itoc(i)//'" from target',
     $               nwks,ci(iv),0,' ')
         ie=ie+1
         iv=iv+nwks
 30   continue
c
c     make the required contracted CI vectors from the
c     target vectors and the number of L**2 orbitals.
c
      if(prnt(1)) then
         write(iout,4) ( ci(ie), ie=ciene, ciene+roots-1)
      endif
      if(logkey(ops,'scattering=split',.false.,' ')) then
         nsplit=intkey(ops,'scattering=nsplit',0,' ')
         if(nsplit.eq.0 .or. nsplit.gt.20) then
            write(iout,*)'m942:  nsplit = 0 or nsplit gt 20 ',nsplit
            call lnkerr(' m942: nsplit ')
         end if
         call intarr(ops,'scattering=split',split,nsplit,' ')
         write(iout,*)' expanding ci vector with partitioned'
     $              //' virtual space'
         call svectr(ci(civec),pvec,roots,nwks,nconf,nl2,
     #               nsplit,split)
      else
         call pvectr(ci(civec),pvec,roots,nwks,nconf,nl2)
      end if
      call getmem(-nwci,pci,idum,'civec',idum)
      maxcor=maxcor+nwci
      call iosys('read character "hamiltonian filename" '//
     $           'from rwf',0,0,0,namham)
      call iosys('open hamiltonian as old',0,0,0,namham)
      call iosys('rewind all on hamiltonian',0,0,0,' ')
c
c     build the compliment of the contracted P space.
c
      call iosys('create real "P-space vectors" on hamiltonian',
     #            nconf*nconf,0,0,' ')
      call iosys('write real "P-space vectors" to hamiltonian '//
     #           'without rewinding',nconf*npvec,pvec,0,' ')
      if(pcomp) then
         ncomp = nwksp-npvec
         if(ncomp.ne.0) then
            lcvec = 1
            need=wpadti(lcvec+nconf*ncomp)
            call getmem(need,pcvec,nwcvec,'cvec',0)
            maxcor=maxcor-nwcvec
            s = 1
            eig=s+nwksp*nwksp
            scr=eig+nwksp
            work=scr+nconf*nwksp
            need=wpadti(work+5*nwksp)
            call getmem(need,pjunk,nwjunk,'junk',0)             
            maxcor=maxcor-nwjunk
            call vcomp(cvec,pvec,rjunk(s),rjunk(eig),rjunk(work),
     #                 rjunk(scr),nconf,nwksp,nwksq,
     #                 npvec,ncomp,prnt(2))
            call iosys('write real "P-space vectors" to hamiltonian '//
     #                 'without rewinding',nconf*ncomp,cvec,0,' ')
            call getmem(-nwjunk,pjunk,idum,'junk',idum)             
            call getmem(-nwcvec,pcvec,idum,'cvec',idum)
            maxcor=maxcor+nwjunk+nwcvec
         endif
         call getmem(-nwpvec,ppvec,idum,'pvec',idum)
         maxcor=maxcor+nwpvec
         npvec=npvec+ncomp
         vec=1
         need=wpadti(vec+npvec*nconf)
         call getmem(need,ppvec,nwpvec,'pvec',0)
      endif
      call iosys('write real "number of P-space vectors" to '//
     #           'hamiltonian',1,npvec,0,' ')
      call iosys('read real "P-space vectors" from hamiltonian',
     #            nconf*npvec,pvec,0,' ')
      call iosys('does "frozen core energy" exist on rwf',0,0,0,
     $            answer)
      if (answer.eq.'no') then
         fzcore=0.0d+00
      else
         call iosys('read real "frozen core energy" from rwf',1,
     #               fzcore,0,' ')
      end if
      call iosys('read integer "buffer size" from hamiltonian',
     $            1,lnbuf,0,' ')
      call iosys('read integer "number of elements" from '//
     #           'hamiltonian',1,nonz,0,' ')
      nodsk=.false.
      if(nonz.le.lnbuf) then
         nodsk=.true.
         lnbuf=nonz
      end if
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)
      rbuf=1
      ibuf=wpadti(rbuf+lnbuf)
      need=ibuf+lnbuf+lnbuf
      call getmem(need,pbuf,nwbuf,'buffers',0) 
      diag=1
      spp=diag+nconf
      rhpp=spp+nwksp*npvec 
      need=wpadti(rhpp+npvec*npvec) 
      if(need.gt.maxcor) then
         call lnkerr('not enough core to continue')
      endif
      call getmem(need,phpp,nwpp,'hpp',0) 
      call iosys('read real diagonals from hamiltonian',nconf,
     1            hpp(diag),0,' ')
      call sadd(hpp(diag),hpp(diag),rep,nconf)
      noopt=logkey(ops,'scattering='//scatyp(1:len)
     #                //'=no-optical-potential',.false.,' ')
      if(noopt) then
         call hampp(abuf(ibuf),buf(rbuf),hpp(diag),hpp(rhpp),
     #              pvec,hpp(spp),energy,header,nwksp,nwksq,
     #              nconf,npvec,lnbuf,nonz,nodsk,nen,prnt(3))

         call getmem(-nwpp,phpp,idum,'hpp',idum)
      else
         call iosys('write real '//header(8)//' to hamiltonian',
     #               nwksq,hpp(diag),0,' ')
c
c        if all of the non zero Hamiltonian matrix elements
c         can be held in core, read them in.
c

         maxcor=maxcor-nwpp
c
c        see if we can also do hqp on this pass
c
         rhqp=1
         need=wpadti(rhqp+nwksq*npvec)
         if(need.lt.maxcor) then
            call getmem(need,phqp,nqp,'hqp',0)
c
c           this is easy.  we can do the whole thing in one shot.
c
            write(iout,5)
            call ham(abuf(ibuf),buf(rbuf),hpp(diag),hpp(rhpp),
     #               hqp(rhqp),pvec,hpp(spp),energy,header,nwksp,
     #               nwksq,nconf,npvec,lnbuf,nonz,nodsk,nen,prnt(3))
            call getmem(-nqp,phqp,idum,'hqp',idum)
            call getmem(-nwpp,phpp,idum,'hpp',idum) 
            call getmem(-nwpvec,ppvec,idum,'pvec',idum)
            maxcor = maxcor + nwpvec + nwpp + nqp
         else
c
c           we need to make more than one pass.  first we will
c           calculate hpp and then after memory is reduced do
c           hpq.
c         
            call hampp(abuf(ibuf),buf(rbuf),hpp(diag),hpp(rhpp),
     #                 pvec,hpp(spp),energy,header,nwksp,
     #                 nwksq,nconf,npvec,lnbuf,nonz,nodsk,nen,prnt(3))
            call getmem(-nwpp,phpp,idum,'hpp',idum)
            maxcor = maxcor + nwpp
c
c           how many vectors can we get at a pass
c
            number=iadtwp(maxcor)/nwksq
            if(number.ge.npvec) then
               trips=1
               number=npvec
            else
               trips = npvec/number
               left=npvec-trips*number
            endif
            need=number*nwksq
            rhqp=1
            need=wpadti(rhqp+need)
c
            call getmem(need,phpq,nqp,'hqp',0)
            locat=1
            do 40 trp=1,trips
               call hamqp(abuf(ibuf),buf(rbuf),hqp(rhqp),pvec(locat),
     #                    header,nwksp,nwksq,nconf,number,lnbuf,nonz,
     #                    nodsk,prnt(4))
               locat=locat+number*nwksq
 40         continue
            if(left.ne.0) then   
               call hamqp(abuf(ibuf),buf(rbuf),hqp(rhqp),pvec(locat),
     #                    header,nwksp,nwksq,nconf,left,lnbuf,
     #                    nonz,nodsk,prnt(4))
            end if
            call getmem(-nqp,phqp,idum,'hqp',idum)
            call getmem(-nwpvec,ppvec,idum,'pvec',idum)
            maxcor = maxcor + nwpvec + nqp
         endif
c
c     
         rbufq=1
         ibufq=wpadti(rbufq+lnbuf)
         diagqq=iadtwp(ibufq+lnbuf+lnbuf) 
         need=wpadti(diagqq+nwksq)
         call getmem(need,qbuf,nwqq,'hqq',0)
         call hamqq(abuf(ibuf),buf(rbuf),abufq(ibufq),bufq(rbufq),
     #              bufq(diagqq),header,nwksq,lnbuf,nonz,nodsk,prnt(5))
         call getmem(-nwqq,qbuf,idum,'hqq',idum)      
         call getmem(-nwbuf,pbuf,idum,'buffers',idum) 
         maxcor=maxcor+nwqq+nwbuf
c
c        allocate memory for optical potential calculation
c
c        force an in-core inversion.  if not possible quit.
c
         incore=logkey(ops,'scattering='//scatyp(1:len)
     #                 //'=force-in-core',.false.,' ')
         outcore=logkey(ops,'scattering='//scatyp(1:len)
     #                 //'=force-out-of-core',.false.,' ')
         if(incore) then         
c
c           calculate the memory needed for an incore inversion
c
            call incmem(rbufq,ibufq,diagqq,rhqq,rhqp,xqp,rhopt,
     #                  rtmp,ipvt,header,lnbuf,nwksq,npvec,ntri,
     #                  nonzqq,nonzqp,nodsk,need,more,maxcor)
            if(more.gt.maxcor) then
c
c              cannot get enough memory for an forced incore inversion
c           
c
               write(iout,6)
               call lnkerr('not enough core:quit')
            else
               call getmem(need,qbuf,nwbuf,'buffers',0) 
               call getmem(more,phopt,nwopt,'hopt',0)
c
c                   in core inversion
c
               call optdir(bufq(diagqq),opt(rhqq),opt(rhqp),opt(xqp),
     #                     opt(rhopt),opt(rtmp),
     #                     bufq(rbufq),abufq(ibufq),
     #                     energy,aopt(ipvt),header,lnbuf,nwksq,npvec,
     #                     ntri,nen,nonzqq,nonzqp,nodsk,prnt(4)) 
             endif
         else if(outcore) then
               write(iout,7)
               call outmem(rbufq,ibufq,diagqq,rhqp,rhopt,vec,hvec,
     #                     resid,t,mat,mattmp,b,btmp,list,header,
     #                     lnbuf,nwksq,npvec,nonzqq,nonzqp,nodsk,
     #                     mxdvd,need,more,maxcor)
               call getmem(need,qbuf,nwbuf,'buffers',0)
               call getmem(more,phopt,nwopt,'hopt',0)
               call optit(bufq(diagqq),opt(rhqp),opt(rhopt),
     #                    bufq(rbufq),abufq(ibufq),energy,
     #                    opt(vec),opt(hvec),opt(resid),opt(mat),
     #                    opt(mattmp),opt(b),opt(btmp),opt(t),
     #                    aopt(list),cnverg,thresh,nwksq,npvec,nen,
     #                    nonzqq,nonzqp,mxiter,mxdvd,lnbuf,nodsk,
     #                    header,prnt(4),prndvd)
         else 
c
c                   let the code decide
c
            call incmem(rbufq,ibufq,diagqq,rhqq,rhqp,xqp,rhopt,rtmp,
     #                  ipvt,header,lnbuf,nwksq,npvec,ntri,
     #                  nonzqq,nonzqp,nodsk,need,more,maxcor)
            if(more.lt.maxcor) then
               call getmem(need,qbuf,nwbuf,'buffers',0) 
               call getmem(more,phopt,nwopt,'hopt',0)
c
c                   in core inversion
c
               call optdir(bufq(diagqq),opt(rhqq),opt(rhqp),opt(xqp),
     #                     opt(rhopt),opt(rtmp),
     #                     bufq(rbufq),abufq(ibufq),
     #                     energy,aopt(ipvt),header,lnbuf,nwksq,npvec,
     #                     ntri,nen,nonzqq,nonzqp,nodsk,prndvd) 
            else 
               write(iout,7)
               call outmem(rbufq,ibufq,diagqq,rhqp,rhopt,vec,hvec,
     #                     resid,t,mat,mattmp,b,btmp,list,header,
     #                     lnbuf,nwksq,npvec,nonzqq,nonzqp,nodsk,
     #                     mxdvd,need,more,maxcor)
               call getmem(need,qbuf,nwbuf,'buffers',0)
               call getmem(more,phopt,nwopt,'hopt',0)
               call optit(bufq(diagqq),opt(rhqp),opt(rhopt),
     #                    bufq(rbufq),abufq(ibufq),energy,
     #                    opt(vec),opt(hvec),opt(resid),opt(mat),
     #                    opt(mattmp),opt(b),opt(btmp),opt(t),
     #                    aopt(list),cnverg,thresh,nwksq,npvec,nen,
     #                    nonzqq,nonzqp,mxiter,mxdvd,lnbuf,nodsk,
     #                    header,prnt(4),prndvd)
            endif
         endif
      endif
c                           END
      call iosys('rewind all on hamiltonian read-and-write',
     #            0,0,0,' ')
      call iosys('close hamiltonian',0,0,0,' ')         
      call chainx(0)
c
      stop
 1    format(//,25x,'optical potential calculation',/,5x,
     1            'no. target orbitals                   = ',i7,/,5x,
     2            'no. target roots                      = ',i7,/,5x,
     3            'no. L**2 p space orbitals             = ',i7,/,5x,
     4            'no. target configurations             = ',i7,/,5x,
     5            'no. p-space configurations            = ',i7,/,5x,
     6            'no. q-space configurations            = ',i7,/,5x,
     7            'total no. configurations              = ',i7,/,5x,
     8            'no. contracted p-space vectors        = ',i7,/,5x,
     9            'reference energy(Ground State Target) = ',e15.8)
 2    format(/,10x,'relative scattering energies',(/,5x,5(e15.8,1x)))
 3    format(/,10x,'absolute scattering energies',(/,5x,5(e15.8,1x)))
 4    format(/,10x,'reading in target ci vectors and energies',/,10x,
     1             'energies = ',(/,10x,5(e15.8,1x)))
 5    format(/,10x,'entire calculation can be done in one pass'
     #       /,10x,'everything will fit in available core')   
 6    format(/,10x,'incore inversion of HQQ forced with insufficient '
     #             'core. quit')
 7    format(/,10x,'out of core inversion of HQQ')
      end
