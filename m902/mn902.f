*deck @(#)mn902.f	5.1  11/6/94
      subroutine mn902(a,z,maxcor,oper,cvec,svec,tunit,tdunit,
     $                 calc,mcroot,vn)
c
c***begin prologue     mn902
c***date written       860000   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci, hamiltonian matrix
c***author             saxe, paul (lanl)
c***source             @(#)mn902.f	5.1   11/6/94
c
c***purpose            to construct the hamiltonian matrix,
c                      or portions thereof.
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       mn902
c
c
      implicit integer (a-z)
c
      character*8 prtflg
      character*4 itoc
      character*128 namham
      character*4096 ops
      character*(*) oper
      character*(*) tunit
      character*(*) tdunit
      character*(*) calc
      character*(*) vn
      character*3 answer
      character*16 unit
      character*16 dunit
      integer a(maxcor)
      logical logkey
      logical doci
      real*8 cvec(*)
      real*8 svec(*)
      real*8 z(*)
      real*8 rep
      real*8 fzcore
      real*8 cutoff
      real*8 fpkey
c
      common /io/ inp,iout
c
c
      unit=tunit
      dunit=tdunit
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
c
c
c     ----- get the constants for dividing up core -----
c
      call iosys('read integer ngroup from rwf',1,ngroup,0,' ')
      call iosys('read integer "symmetries in ci" from rwf',
     $            1,nsym,0,' ')
      call iosys('read integer norbs from rwf',1,norbs,0,' ')
      call iosys('read integer nrows from rwf',1,nrows,0,' ')
      call iosys('read integer nlevs from rwf',1,nlevs,0,' ')
      call iosys('read integer orbfrm from rwf',1,orbfrm,0,' ')
      call iosys('read integer numij from rwf',1,nnp,0,' ')
      call iosys('read integer nmax from rwf',1,nmax,0,' ')
      call iosys('read integer nijvir from rwf',1,nnpvir,0,' ')
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,rep,0,' ')
c
c
      if (oper.eq.'ci') then
         if(logkey(ops,'ci=space=p',.false.,' ')) then
            doci=.false.
            sirow=2
            sjrow=2
            call iosys('read integer nwksp from rwf',1,nwks,0,' ')
            if (prtflg.ne.'minimum') then
               write (iout,101) nwks
 101           format (5x,'forming the p-space hamiltonian of ',i5,
     #              ' configurations')
            end if
         else if(logkey(ops,'ci=space=q',.false.,' ')) then
            doci=.false.
            sirow=1
            sjrow=1
            call iosys('read integer nwksq from rwf',1,nwks,0,' ')
            if (prtflg.ne.'minimum') then
               write (iout,102) nwks
 102           format (5x,'forming the q-space hamiltonian of ',i5,
     #              ' configurations')
            end if
         else if (logkey(ops,'ci=space=pq',.false.,' ')) then
            doci=.false.
            sirow=1
            sjrow=2
            call iosys('read integer nwksp from rwf',1,nwksp,0,' ')
            call iosys('read integer nwksq from rwf',1,nwksq,0,' ')
            nwks=1
            if (prtflg.ne.'minimum') then
               write (iout,202) nwksp,nwksq
 202           format (5x,'forming the p-q block of the hamiltonian',/,
     $              10x,i4,' x ',i5,' configurations')
            end if
         else if(logkey(ops,'ci=all',.false.,' ')) then
            doci=.true.
            sirow=0
            sjrow=0
            call iosys('read integer nwksp from rwf',1,nwksp,0,' ')
            call iosys('read integer nwksq from rwf',1,nwksq,0,' ')
            nwks=nwksp+nwksq
            call iosys('write integer nwks to rwf',1,nwks,0,' ')
            if (prtflg.ne.'minimum') then
               write (iout,192) nwks
 192           format (5x,'diagonalizing the p+q hamiltonian of ',i5,
     #              ' configurations')
            end if
         else
            doci=.true.
            sirow=1
            sjrow=1
            call iosys('read integer nwks from rwf',1,nwks,0,' ')
            if (prtflg.ne.'minimum') then
               write (iout,103) nwks
 103           format (5x,'forming the entire hamiltonian of ',i5,
     #              ' configurations')
            end if
         end if
c
c        ----- buffer length for storing hamiltonian elements -----
c
         lnbuf=intkey(ops,'ci=buffer-length',20000,' ')
c
         call iosys('write integer "ci buffer length" to rwf',
     $               1,lnbuf,0,' ')
c
         cutoff=fpkey(ops,'ci=cutoff',1.0d-9,' ')
      else if (oper.eq.'density') then
         sirow=1
         sjrow=1
         call iosys('read integer nwks from rwf',1,nwks,0,' ')
         if (prtflg.ne.'minimum') then
            write (iout,104) nwks
 104        format (5x,'forming the density matrix from ',i5,
     $           ' configurations')
         end if
      else if (oper.eq.'transition density') then
         sirow=1
         sjrow=1
         call iosys('read integer nwks from rwf',1,nwks,0,' ')
         if (prtflg.ne.'minimum') then
            write (iout,105) nwks
 105        format (5x,'forming the averaged transition density matrix',
     $           ' from ',i5,' configurations')
         end if
      else if (oper.eq.'h times c') then
         sirow=1
         sjrow=1
         call iosys('read integer nwks from rwf',1,nwks,0,' ')
         if (prtflg.ne.'minimum') then
            write (iout,106) nwks
 106        format (5x,'forming h.c for ',i5,' configurations')
         end if
      else
         write(iout,*) 'oper:',oper
         call lnkerr('bad operation passed into mn902')
      end if
c
c     determine how many roots are desired.
      if(calc.eq.'mcscf') then
         nroots=mcroot
      else
         nroots=intkey(ops,'ci=nroots',1,' ')
         nroots=min(nwks,nroots)
      endif
c
c     ----- allocate core for drt arrays, integrals, vectors, etc. -----
c
      imngrp=1
      imxgrp=imngrp+ngroup
      jmngrp=imxgrp+ngroup
      jmxgrp=jmngrp+ngroup
      arc=jmxgrp+ngroup
      wt=arc+4*nrows
      nlwks=wt+4*nrows
      ijadd=nlwks+nrows
      ijgrp=ijadd+nnp
      kadd=ijgrp+nnp
      ladd=kadd+norbs*nsym
      orbsym=ladd+norbs*nsym
      b=orbsym+norbs
      refwt=b+nrows
      refarc=refwt+nlevs
      refb=refarc+nlevs
      irowsv=refb+nlevs
      jrowsv=irowsv+nlevs
      segsv=jrowsv+nlevs
      pagesv=segsv+nlevs
      iwtsv=pagesv+nlevs
      jwtsv=iwtsv+nlevs
      traksv=jwtsv+nlevs
c
      ints=iadtwp(traksv+nlevs)
      acoef=ints+nmax
      bcoef=acoef+nlevs
      if (oper.eq.'ci') then
         c=bcoef+nlevs
         rbuf=c+nwks
         ibuf=wpadti(rbuf+lnbuf)
         need=ibuf+2*lnbuf
      else if (oper.eq.'density') then
         c=bcoef+nlevs
         need=wpadti(c+nwks)
      else if (oper.eq.'transition density') then
         need=wpadti(bcoef+nlevs)
      else if (oper.eq.'h times c') then
         need=wpadti(bcoef+nlevs)
      end if
c
c     ----- get core, and then read in arrays -----
c
      call getscm(need,a,maxcor,'mn902',0)
c

      call iosys('read integer arc from rwf',4*nrows,a(arc),0,' ')
      call iosys('read integer weight from rwf',4*nrows,a(wt),0,' ')
      call iosys('read integer nlwks from rwf',nrows,a(nlwks),0,' ')
      call iosys('read integer ijadd from rwf',nnp,a(ijadd),0,' ')
      call iosys('read integer ijgrp from rwf',nnp,a(ijgrp),0,' ')
      call iosys('read integer kadd from rwf',
     $            norbs*nsym,a(kadd),0,' ')
      call iosys('read integer ladd from rwf',norbs*nsym,
     $            a(ladd),0,' ')
      call iosys('read integer orbsym from rwf',norbs,a(orbsym),0,' ')
      call iosys('read integer b from rwf',nrows,a(b),0,' ')
c
c     ----- fix up the symmetry array and add group mins and maxs -----
c
      call fixdrt(a(orbsym),norbs,a(ijgrp),nnp,a(imngrp),a(imxgrp),
     $            a(jmngrp),a(jmxgrp),ngroup)
c
      call iosys('does "frozen core energy" exist on rwf',0,0,0,
     $            answer)
      if (answer.eq.'no') then
         fzcore=0.0d+00
      else
         call iosys('read real "frozen core energy" from rwf',1,
     #               fzcore,0,' ')
      end if
c
      if (oper.eq.'ci') then
c
c        ----- open a scratch file for the hamiltonian -----
c
         junk=min(1600000,2*nwks**2)
         junk=max(junk,100000)
         if (logkey(ops,'ci=keep',.false.,' ').or.
     $       logkey(ops,'ci=space=q',.false.,' ').or.
     $       logkey(ops,'kohn',.false.,' ')) then
c
            call iosys('read character "hamiltonian filename" '//
     $                 'from rwf',0,0,0,namham)
            call iosys('open hamiltonian as new on ssd',
     $                 junk,0,0,namham)
         else
            call iosys('open hamiltonian as scratch on ssd',
     $                 junk,0,0,' ')
         end if
c
c        ----- form and write out the packed matrix -----
c
         call hamilt(a(arc),a(wt),a(nlwks),a(ijadd),a(kadd),a(ladd),
     #        a(orbsym),a(refwt),a(refarc),a(b),a(refb),z(ints),z(c),
     #        nrows,norbs,nlevs,orbfrm,nsym,nmax,nwks,nnp,
     #        a(irowsv),a(jrowsv),a(segsv),a(pagesv),
     $        a(iwtsv),a(jwtsv),a(traksv),z(acoef),z(bcoef),
     $        lnbuf,a(ibuf),z(rbuf),cutoff,ngroup,a(imngrp),a(imxgrp),
     $        a(jmngrp),a(jmxgrp),sirow,sjrow,ops,unit,prtflg)
c
c
         if (logkey(ops,'ci=space=p',.false.,' ')) then
            ibuf=1
            t1=1
            rbuf=iadtwp(ibuf+2*lnbuf)
            h=max(rbuf+lnbuf,t1+nwks)
            need=wpadti(h+nwks**2)
c
            call getscm(need,a,maxcor,'hpp',0)
c
            call fmh(lnbuf,a(ibuf),z(rbuf),z(h),nwks,z(t1),ops)
c
            call iosys('write real h(pp) to rwf',nwks**2,z(h),0,' ')
         end if
c
         if (logkey(ops,'ci=space=q=form',.false.,' ')) then
            ibuf=1
            t1=1
            rbuf=iadtwp(ibuf+2*lnbuf)
            h=max(rbuf+lnbuf,t1+nwks)
            need=wpadti(h+nwks**2)
c
            call getscm(need,a,maxcor,'hqq',0)
c
            call fmh(lnbuf,a(ibuf),z(rbuf),z(h),nwks,z(t1),ops)
c
            call iosys('write real h(qq) to rwf',nwks**2,z(h),0,' ')
         end if
c
         if (logkey(ops,'ci=space=pq',.false.,' ')) then
c
c            ----- read in and form hpq -----
c
            ibuf=1
            rbuf=iadtwp(ibuf+2*lnbuf)
            hpq=rbuf+lnbuf
            need=wpadti(hpq+nwksp*nwksq)
c
            call getscm(need,a,maxcor,'hpq',0)
c
            call fmhpq(lnbuf,a(ibuf),z(rbuf),z(hpq),nwksp,nwksq,ops)
c
            call iosys('write real h(pq) to rwf',
     $                  nwksp*nwksq,z(hpq),0,' ')
c
         else if (logkey(ops,'ci=diagonalize',doci,' ')) then
c
c           ----- reallocate core -----
c
            if(logkey(ops,'ci=kohn',.false.,' ')) return
c
            ibuf=1
            rbuf=iadtwp(ibuf+2*lnbuf)
c
c           ----- try for in-core diagonalization if possible -----
c
            nnpwks=nwks*(nwks+1)/2
c
            eigvec=1
            eigval=eigvec+nwks**2
            t1=eigval+nwks
            t2=t1+nwks
            h=max(t2+nwks,rbuf+lnbuf)
            need=wpadti(h+nnpwks)
c
            call getscm(0,a,canget,'mn902',0)
c
            if (logkey(ops,'ci=givens-householder',.true.,' ').and.
     $           need.le.canget) then
               call getscm(need,a,maxcor,'mn902',0)
c
c              ----- n**3 diagonalization -----
c
               call incore(lnbuf,a(ibuf),z(rbuf),nwks,nroots,z(eigvec),
     $              z(eigval),z(t1),z(t2),z(h),nnpwks,ops,rep,fzcore,
     $              prtflg,calc)
            else
c
c              ----- out-of-core davidson diagonalization -----
c
               if(logkey(ops,'mcscf',.false.,' ')) then
                  mxiter=intkey(ops,'mcscf=ci=iterations',20,' ')
                  junk=mxiter
                  mxiter=intkey(ops,'mcscf=ci=cycles',junk,' ')
               else
                  mxiter=intkey(ops,'ci=iterations',15,' ')
                  junk=mxiter
                  mxiter=intkey(ops,'ci=cycles',junk,' ')
               end if
c
               nwksg=intkey(ops,'ci=guess-size',0,' ')
               nattim=intkey(ops,'ci=nroots-at-a-time',nroots,' ')
               call iosys('read integer "number of elements"'//
     $                    'from hamiltonian',1,ntotal,0,' ')
c
               root=rbuf+lnbuf
               dvdvec=root+mxiter
               dvdmat=dvdvec+mxiter**2
               c=dvdmat+mxiter*(mxiter+1)/2
c
               call getscm(0,a,canget,'mn902',0)
c
               left=iadtwp(canget)-c
               maxvec=min(nattim,left/2/nwks)
               if (maxvec.lt.1) then
                  call lnkerr('not enough memory for two vectors !')
               end if
c
               s=c+nwks*maxvec
               need=wpadti(s+nwks*maxvec)
c
               call getscm(need,a,maxcor,'mn902',0)
c
c              ----- core allocation for guess routine -----
c
               nwksg=min(nwks,nwksg)
               if (nwksg.le.0) then
                  do 40 nwksg=10,500,10
                     if (wpadti(rbuf+lnbuf+nwks+nwks+nwksg+nwksg**2+
     $                    nwksg*(nwksg+1)/2)+nwks.gt.maxcor) go to 45
 40               continue
                  nwksg=500
                  go to 50
 45               continue
                  nwksg=nwksg-10
 50               continue
               end if
c
               nwksg=min(nwks,nwksg)
               nnpg=nwksg*(nwksg+1)/2
c
               cguess=rbuf+lnbuf
               sguess=cguess+nwks
               hguess=sguess+nwks
               smlvec=hguess+nnpg
               smlval=smlvec+nwksg**2
               ptgues=wpadti(smlval+nwksg)
               need=ptgues+nwks
c
               call getscm(need,a,maxcor,'guess',0)
c
               nguess=intkey(ops,'ci=nguess',nroots,' ')
               junk=min(nguess*nwks,16000000)
               call iosys('open guess as scratch',junk,0,0,' ')
c
               call guess2(a(ptgues),nwks,a(ibuf),z(rbuf),lnbuf,
     $              z(hguess),
     $              nwksg,nnpg,z(cguess),z(sguess),z(smlvec),z(smlval),
     $              ops,ntotal,rep+fzcore,prtflg)
c
               call david(z(c),z(s),nwks,maxvec,z(root),nroots,
     $              z(dvdvec),
     $              z(dvdmat),mxiter,mxiter*(mxiter+1)/2,rep,fzcore,ops,
     $              a(ibuf),z(rbuf),lnbuf,prtflg,calc)
c
               call iosys('destroy guess',0,0,0,' ')
            end if
         end if
c
         call iosys('close hamiltonian',0,0,0,' ')
c
      else if (oper.eq.'density') then
c
c        ----- density matrix section -----
c
         froot=max(mcroot,1)
         if(calc.eq.'mcscf') then
            call iosys('read real "mc root '//itoc(froot)//'" from rwf',
     $                  nwks,z(c),0,' ')
         else
            call iosys('read real "ci root '//itoc(froot)//'" from rwf',
     $                  nwks,z(c),0,' ')
         endif
c
         call dm02(a(arc),a(wt),a(nlwks),a(ijadd),a(kadd),a(ladd),
     #        a(orbsym),a(refwt),a(refarc),a(b),a(refb),z(ints),z(c),
     #        nrows,norbs,nlevs,orbfrm,nsym,nmax,nwks,nnp,
     #        a(irowsv),a(jrowsv),a(segsv),a(pagesv),
     $        a(iwtsv),a(jwtsv),a(traksv),z(acoef),z(bcoef),
     $        lnbuf,a(ibuf),z(rbuf),cutoff,ngroup,a(imngrp),a(imxgrp),
     $        a(jmngrp),a(jmxgrp),sirow,sjrow,ops,dunit)
c
c
      else if (oper.eq.'transition density') then
c
c        ----- transition density matrix section -----
c
         call tdm02(a(arc),a(wt),a(nlwks),a(ijadd),a(kadd),a(ladd),
     #        a(orbsym),a(refwt),a(refarc),a(b),a(refb),z(ints),
     $        cvec,svec,
     #        nrows,norbs,nlevs,orbfrm,nsym,nmax,nwks,nnp,
     #        a(irowsv),a(jrowsv),a(segsv),a(pagesv),
     $        a(iwtsv),a(jwtsv),a(traksv),z(acoef),z(bcoef),
     $        lnbuf,a(ibuf),z(rbuf),cutoff,ngroup,a(imngrp),a(imxgrp),
     $        a(jmngrp),a(jmxgrp),sirow,sjrow,ops,dunit)
c
      else if (oper.eq.'h times c') then
c
c        ----- product of a vector times hamiltonian section -----
c
         call hc02(a(arc),a(wt),a(nlwks),a(ijadd),a(kadd),a(ladd),
     #        a(orbsym),a(refwt),a(refarc),a(b),a(refb),z(ints),
     $        cvec,svec,
     #        nrows,norbs,nlevs,orbfrm,nsym,nmax,nwks,nnp,
     #        a(irowsv),a(jrowsv),a(segsv),a(pagesv),
     $        a(iwtsv),a(jwtsv),a(traksv),z(acoef),z(bcoef),
     $        lnbuf,a(ibuf),z(rbuf),cutoff,ngroup,a(imngrp),a(imxgrp),
     $        a(jmngrp),a(jmxgrp),sirow,sjrow,ops,unit)
c
      end if
c
c
      return
      end
