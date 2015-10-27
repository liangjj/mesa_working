*deck @(#)mn940.f	1.4  8/3/91
      subroutine mn940(ops)
c
c***begin prologue     mn940
c***date written       000811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci, hamiltonian matrix
c***author             saxe, paul (lanl)
c***source             @(#)m940.f	1.4   8/3/91
c
c***purpose            to construct the hamiltonian matrix,
c                      or portions thereof.
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       mn940
c
c
      implicit integer (a-z)
c
      character*(*) ops
      character*24 mcorci
      character*8 citype
      character*128 gints
      character*32 title
      logical logkey, ppon, qqon, pqon, pandqon, normal, partit
      logical conham, hamon, contrq
      logical diaham, prnt
      logical target
      character*8 prtflg
      character*5 dsk
      character*4 itoc
      character*128 namham, nmpart
      character*3 answer
      character*16 unit
      integer a
      logical dvd, dagtyp, pfile
      real*8 z
      real*8 rep
      real*8 fzcore
      real*8 cutoff
      real*8 fpkey
      dimension title(3)
      pointer(p,a(1)),(p,z(1))
      data ppon/.false./
      data pqon/.false./
      data qqon/.false./
      data pandqon/.false./
      data normal/.false./
      data conham/.false./
      data diaham/.false./ 
      data target/.false./
      data partit/.false./
c
      common /io/ inp,iout
c
c
      mcorci='ci'
      citype='m940'
c
c     ----- open  the integral file ----
c
      write(iout,*) 'm940: general ci'
      call iosys ('read character "guga integral filename" '//
     $            'from rwf',-1,0,0,gints)
      call iosys('open gints as old',0,0,0,gints)
      unit='gints'
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      call iosys('read integer "number of basis functions" from rwf',
     1            1,nbf,0,' ')
c
      pfile=.true.
c
c     decide on the type of ci you are doing
c
c
c     determine if this is a target state calculation
c
      dsk='rwf'
      call iosys('read character "partitioning filename" '//
     $           'from rwf',0,0,0,nmpart)
      call iosys('open hpart as unknown',0,0,0,nmpart)
      if ( logkey(ops,'ci=p',.false.,' ').or.
     1     logkey(ops,'ci=q',.false.,' ').or.
     2     logkey(ops,'ci=pq',.false.,' ').or.
     3     logkey(ops,'ci=p-and-q',.false.,' ').or.
     4     logkey(ops,'ci=all',.false.,' ') ) then
           dsk='hpart'
      endif
      write(iout,*) '          reading information from '//dsk
      citype='m940'
      if( logkey(ops,'ci=target',.false.,' ') ) then      
          target=.true.
          normal=.true.
          if(logkey(ops,'ci=form',.false.,' ')) then
             conham=.true.
          endif
          if(logkey(ops,'ci=unformatted',.false.,' ')) then
             hamon=.true.
          endif
          if(logkey(ops,'ci=diagonalize',.false.,' ')) then
             diaham=.true.
          endif
c
c         get the necessary input from ops
c
          call iosys('does "internal orbitals" exist on '//dsk,0,0,
     1                0,answer)
          if(answer.eq.'yes') then
             call iosys('read integer "internal orbitals" from '//dsk,
     1                   1,trgorb,0,' ')
          else
             trgorb=intkey(ops,'internal-orbitals',0,' ')
          endif
          if(trgorb.eq.0) then
             write(iout,*) ' error. number of internal or target '//
     1                     'orbitals not set' 
             call lnkerr('nsmall not set')
          endif
          call iosys('read integer nwks from '//dsk,1,nwks,0,' ')
          call iosys('read integer nwksp from '//dsk,1,nwksp,0,' ')
          call iosys('read integer nwksq from '//dsk,1,nwksq,0,' ')
          roots=intkey(ops,'ci=target-roots',1,' ')
          twalks=nwksp
          porb=nwksp/nwks
          call iosys('rewind all on hpart read-and-write',0,0,0,' ')
          call iosys('write integer "target roots" to hpart',1,
     1                roots,0,' ')
          call iosys('write integer "target orbitals" to hpart',1,
     1                trgorb,0,' ')
          call iosys('write integer "target walks" to '//
     1               'hpart',1,nwks,0,' ')
          call iosys('write integer nwksp to hpart',1,nwksp,0,' ') 
          call iosys('write integer nwksq to hpart',1,nwksq,0,' ') 
          call iosys('write integer roots to rwf',1,roots,0,' ')
          write(iout,1) trgorb, porb, nwks, twalks, roots
      elseif( logkey(ops,'ci=p',.false.,' ') ) then
          mcorci='p-ci'
          ppon=.true.
          if(logkey(ops,'ci=p-explicit',.false.,' ')) then
             partit=.true.
          endif
          if(logkey(ops,'ci=form',.false.,' ')) then
             conham=.true.
          endif
          if(logkey(ops,'ci=diagonalize',.false.,' ')) then
             diaham=.true.
          endif
          call iosys('read integer nwksp from '//dsk,1,nwksp,0,' ')
          call iosys('read integer nwksq from '//dsk,1,nwksq,0,' ')
          proots=intkey(ops,'ci=p-space-roots',1,' ')
          proots=min(proots,nwksp)
          nwks=nwksp+nwksq
          write (iout,2) nwks, nwksp, proots
          call iosys('write integer "p space roots" to '//
     1               'rwf',1,proots,0,' ')
      elseif(logkey(ops,'ci=q',.false.,' ')) then
          mcorci='q-ci'
          qqon=.true.
          if(logkey(ops,'ci=q-explicit',.false.,' ')) then
             partit=.true.
          endif
          if(logkey(ops,'ci=form',.false.,' ')) then
            conham=.true.
          endif
          if(logkey(ops,'ci=diagonalize',.false.,' ')) then
             diaham=.true.
             contrq=.false.
             if(logkey(ops,'ci=contract-q',.false.,' ')) then             
                contrq=.true.
             endif
          endif
          call iosys('read integer nwksp from '//dsk,1,nwksp,0,' ')
          call iosys('read integer nwksq from '//dsk,1,nwksq,0,' ')
          qroots=intkey(ops,'ci=q-space-roots',1,' ')
          qroots=min(qroots,nwksq)
          call iosys('write integer "q space roots" to '//
     1               'rwf',1,qroots,0,' ')
          nwks=nwksp+nwksq
          write (iout,3) nwks, nwksq, qroots
      elseif (logkey(ops,'ci=pq',.false.,' ')) then
          pqon=.true.
          partit=.true.
          if(logkey(ops,'ci=form',.false.,' ')) then
             conham=.true.
          endif
          call iosys('read integer nwksp from '//dsk,1,nwksp,0,' ')
          call iosys('read integer nwksq from '//dsk,1,nwksq,0,' ')
          nwks=nwksp+nwksq
          write (iout,4) nwks, nwksp, nwksq
      elseif (logkey(ops,'ci=p-and-q',.false.,' ')) then
          ppon=.true.
          pqon=.true.
          qqon=.true.
          partit=.true.
          if(logkey(ops,'ci=form',.false.,' ')) then
             conham=.true.
          endif
          if(logkey(ops,'ci=diagonalize',.false.,' ')) then
             diaham=.true.
             if(logkey(ops,'ci=contract-q',.false.,' ')) then             
                contrq=.true.
             endif
          endif
          call iosys('read integer nwksp from '//dsk,1,nwksp,0,' ')
          call iosys('read integer nwksq from '//dsk,1,nwksq,0,' ')
          nwks=nwksp+nwksq
          proots=intkey(ops,'ci=p-space-roots',1,' ')
          qroots=intkey(ops,'ci=q-space-roots',1,' ')
          proots=min(proots,nwksp)
          qroots=min(qroots,nwksq)
          write (iout,5) nwks, nwksp, proots, nwksq, qroots
      elseif(logkey(ops,'ci=all',.false.,' ')) then
          if(.not.logkey(ops,'int=drt=key=p-and-q',.false.,' ')) then
              write(iout,*) 'ci=all and p-and-q options'//
     1                      ' must both be on'
              call lnkerr('quit')
          endif
          mcorci='ci'
          pandqon=.true.
          if(logkey(ops,'ci=p-and-q-explicit',.false.,' ')) then
             partit=.true.
          endif
          if(logkey(ops,'ci=form',.false.,' ')) then
             conham=.true.
          endif
          if(logkey(ops,'ci=unformatted',.false.,' ')) then
             hamon=.true.
          endif
          if(logkey(ops,'ci=diagonalize',.false.,' ')) then
             diaham=.true.
          endif
          nvec=intkey(ops,'ci=number-of-test-vectors',0,' ')
          call iosys('read integer nwksp from '//dsk,1,nwksp,0,' ')
          call iosys('read integer nwksq from '//dsk,1,nwksq,0,' ')
          nwks=nwksp+nwksq 
          roots=intkey(ops,'ci=roots',1,' ')
          roots=min(roots,nwks)
          call iosys('write integer nwks to rwf',1,nwks,0,' ')
          call iosys('write integer roots to rwf',1,roots,0,' ')
          write (iout,6) nwks, nwksp, nwksq, roots
      else
          normal=.true.
          nwksp=0
          nwksq=0
          call iosys('read integer nwks from '//dsk,1,nwks,0,' ')
          roots=intkey(ops,'ci=roots',1,' ')
          roots=min(roots,nwks)
          if(logkey(ops,'ci=form',.false.,' ')) then
             conham=.true.
          endif
          if(logkey(ops,'ci=diagonalize',.false.,' ')) then
             diaham=.true.
          endif
          write(iout,7) nwks, roots
      endif
      call iosys('write character mcorci to rwf',0,0,0,mcorci)
      call iosys('write character "ci used" to rwf',0,0,0,citype)
c
c
c     ----- get the constants for dividing up core -----
c
      call iosys('read integer ngroup from '//dsk,1,ngroup,0,' ')
      call iosys('read integer "symmetries in ci" from '//dsk,
     $            1,nsym,0,' ')
      call iosys('read integer norbs from '//dsk,1,norbs,0,' ')
      call iosys('read integer nrows from '//dsk,1,nrows,0,' ')
      call iosys('read integer nlevs from '//dsk,1,nlevs,0,' ')
      call iosys('read integer orbfrm from '//dsk,1,orbfrm,0,' ')
      call iosys('read integer numij from '//dsk,1,nnp,0,' ')
      call iosys('read integer nmax from '//dsk,1,nmax,0,' ')
      call iosys('read integer nijvir from '//dsk,1,nnpvir,0,' ')
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,rep,0,' ')
c
c
c
c        ----- buffer length for storing hamiltonian elements -----
c
      lnbuf=intkey(ops,'ci=buffer-length',20000,' ')
c
      call iosys('write integer "ci buffer length" to rwf',
     $            1,lnbuf,0,' ')
c
      cutoff=fpkey(ops,'ci=cutoff',1.0d-9,' ')
c
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
      c=bcoef+nlevs
      rbuf=c+max(nwks,nwksp,nwksq)
      ibuf=wpadti(rbuf+lnbuf)
      need=ibuf+2*lnbuf
c
c            Lee, If the previous thing does not work, uncomment these
c            following lines, remake and we will rerun
c
c      write(iout,*) ngroup, nrows, nnp, nsym
c      write(iout,*) norbs, nlevs, traksv, nmax
c      write(iout,*) nwks, nwksp, nwksq, lnbuf
c      write(iout,*) need
c      call lnkerr('quit')
c
c     ----- get core for main arrays and then read them in -----
c
      call getmem(need,p,ngot,'m940',0)
c
      call iosys('read integer arc from '//dsk,4*nrows,a(arc),0,' ')
      call iosys('read integer weight from '//dsk,4*nrows,a(wt),0,' ')
      call iosys('read integer nlwks from '//dsk,nrows,a(nlwks),0,' ')
      call iosys('read integer ijadd from '//dsk,nnp,a(ijadd),0,' ')
      call iosys('read integer ijgrp from '//dsk,nnp,a(ijgrp),0,' ')
      call iosys('read integer kadd from '//dsk,
     $            norbs*nsym,a(kadd),0,' ')
      call iosys('read integer ladd from '//dsk,norbs*nsym,
     $            a(ladd),0,' ')
      call iosys('read integer orbsym from '//dsk,norbs,a(orbsym),0,' ')
      call iosys('read integer b from '//dsk,nrows,a(b),0,' ')
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
c
c        ----- open either a scratch or permanent file for the 
c        -----               hamiltonian             -----
c
      if(pfile) then
         call iosys('read character "hamiltonian filename" '//
     $              'from rwf',0,0,0,namham)
         call iosys('open hamiltonian as new',0,0,0,namham)
      else
         call iosys('open hamiltonian as scratch',0,0,0,' ')
      end if
      call iosys('write integer "buffer size" to hamiltonian',
     $            1,lnbuf,0,' ')
      prnt=logkey(ops,'print=ci=hmatrix',.false.,' ')
c----------------------------------------------------------------------c
c                           BEGIN
c           section to form and write out the packed hamiltonian matrices
c
      title(1)='buffers'
      title(2)='"number of elements"'
      title(3)='diagonals'
      if(.not.normal) then
         sirow=0
         sjrow=0
         write(iout,*) 'calculating h for p and q spaces'
      else
         write(iout,*) 'calculating h'
         sirow=1
         sjrow=1
      endif
c
c         note that the output from hamilt is an unsorted
c         list of non-zero matrix elements with their indices.
c         in addition, more than one piece of a matrix element
c         may appear.  the final matrix element must be summed over
c         all the pieces.  if these are needed explicitly, then
c         either one must do a canonical sort and add the pieces or
c         zero a hamiltonian array of the appropriate size, and
c         accumulate the matrix elements as one reads the list.
c
      call hamilt(a(arc),a(wt),a(nlwks),a(ijadd),a(kadd),a(ladd),
     #            a(orbsym),a(refwt),a(refarc),a(b),a(refb),
     #            z(ints),z(c),nrows,norbs,nlevs,orbfrm,nsym,
     #            nmax,nwks,nnp,a(irowsv),a(jrowsv),
     #            a(segsv),a(pagesv),a(iwtsv),a(jwtsv),
     #            a(traksv),z(acoef),z(bcoef),lnbuf,a(ibuf),
     #            z(rbuf),cutoff,ngroup,a(imngrp),
     #            a(imxgrp),a(jmngrp),a(jmxgrp),sirow,sjrow,
     #            ops,unit,prtflg,title,prnt)
c
c                           END
c
c
c----------------------------------------------------------------------c
c
c
c               free all memory
c  
      call getmem(-ngot,p,idum,'m940',idum)
c
c----------------------------------------------------------------------c
c
c                          BEGIN
c       section to explicitly construct the non-zero hamiltonian 
c       matrix elements of a partitioned hamiltonian.  the same statement
c       made about the hamilt subroutine above applies here.
c
c
      if ( ppon ) then
          ibuf=1
          rbuf=iadtwp(ibuf+2*lnbuf)
          diag=rbuf+lnbuf
          ipbuf=wpadti(diag+nwks)
          hpbuf=iadtwp(ipbuf+2*lnbuf)
          diagp=hpbuf+lnbuf
          need=wpadti(diagp+nwksp)
          if(conham) then
             hpp=diagp+nwksp
             need=wpadti(hpp+nwksp*nwksp)
          endif
c
          call getmem(need,p,ngot,'hpp',0)
c
          title(1)='"p space buffers"'
          title(2)='"number of p space elements"'
          title(3)='"p space diagonals"'
          call hampp(a(ibuf),z(rbuf),z(diag),a(ipbuf),z(hpbuf),
     1               z(diagp),z(hpp),lnbuf,
     2               nwksp,nwksq,conham,title,prnt)
          call getmem(-ngot,p,idum,'hpp',idum)
          if(diaham) then
             call diagh('ci',rep,fzcore,lnbuf,proots,nwksp,
     1                   ops,title)
          endif
      endif
      if ( qqon ) then
          ibuf=1
          rbuf=iadtwp(ibuf+2*lnbuf)
          diag=rbuf+lnbuf
          iqbuf=wpadti(diag+nwks)
          hqbuf=iadtwp(iqbuf+2*lnbuf)
          diagq=hqbuf+lnbuf
          need=wpadti(diagq+nwksq)
          if(conham) then
c
c            see remarks about h(pp)
c
             hqq=diagq+nwksp
             need=wpadti(hqq+nwksq*nwksq)
          endif
c
          call getmem(need,p,ngot,'hqq',0)            
c
          title(1)='"q space buffers"'
          title(2)='"number of q space elements"'
          title(3)='"q space diagonals"'
          call hamqq(a(ibuf),z(rbuf),z(diag),a(iqbuf),z(hqbuf),
     1               z(diagq),z(hqq),lnbuf,
     2               nwksp,nwksq,conham,title,prnt)
          call getmem(-ngot,p,idum,'hqq',idum)            
          if(diaham) then
             call diagh('ci',rep,fzcore,lnbuf,qroots,nwksq,
     1                   ops,title)
          endif
      endif
      if ( pqon ) then
          ibuf=1
          rbuf=iadtwp(ibuf+2*lnbuf)
          ipqbuf=wpadti(rbuf+lnbuf)
          hpqbuf=iadtwp(ipqbuf+2*lnbuf)
          need=wpadti(hpqbuf+lnbuf)
          if(conham) then
             hpq=hpqbuf+lnbuf
             need=wpadti(hpq+nwksp*nwksq)
          endif
c
          call getmem(need,p,ngot,'hpq',0)
          title(1)='"pq space buffers"'
          title(2)='"number of pq space elements"'
          call hampq(a(ibuf),z(rbuf),a(ipqbuf),z(hpqbuf),
     1               z(hpq),lnbuf,nwksp,nwksq,
     2               conham,title,prnt)
          call getmem(-ngot,p,idum,'hpq',idum)
      endif
      if (pandqon.or.normal) then
         if(conham.or.hamon) then
            title(1)='buffers'
            title(2)='"number of elements"'
            title(3)='diagonals'
            ibuf=1
            rbuf=iadtwp(ibuf+2*lnbuf)
            diag=rbuf+lnbuf
            h=diag+nwks
            if(conham) then
               need=wpadti(h+nwks*nwks)
               call getmem(need,p,ngot,'h',0)
c
c              see remarks about h(pp) above
c
               call ham(a(ibuf),z(rbuf),z(diag),z(h),
     1                  lnbuf,nwks,title,prnt)
            elseif(hamon) then
               need=wpadti(h)
               call getmem(need,p,ngot,'h',0)
               call hamunf(a(ibuf),z(rbuf),z(diag),
     1                     lnbuf,nwks,title)         
            endif
            call getmem(-ngot,p,idum,'h',idum)            
         endif
         if(diaham) then
            call diagh('ci',rep,fzcore,lnbuf,roots,nwks,
     1                  ops,title)
            if(nvec.gt.0) then
               nvec=min(nvec,roots)
               ibuf=1
               rbuf=iadtwp(ibuf+2*lnbuf)
               diag=rbuf+lnbuf
               vecs=diag+nwks
               temp=vecs+nvec*nwks            
               need=wpadti(temp+nvec*nwks)            
               title(1)='buffers'
               title(2)='"number of elements"'
               title(3)='diagonals'
               call getmem(need,p,ngot,'test',0)
               call tvmult(z(rbuf),a(ibuf),z(diag),z(vecs),
     1                     z(temp),rep,fzcore,lnbuf,nwks,nvec,title)
               call getmem(-ngot,p,idum,'test',idum)
            endif
         endif
      endif
c
c
      if( partit ) then
c
          if( logkey (ops,'ci=data=ops',0,' ') ) then          
c
              trgorb=intkey(ops,'ci=target-orbitals',0,' ')
              nwkst=intkey(ops,'ci=target-configurations',0,' ')
              porb=intkey(ops,'ci=p-space-orbitals',nwksp,' ')
              troots=intkey(ops,'ci=target-roots',1,' ')
              twalks=nwkst*porb
              qwalks=nwks-twalks
              npvec=troots*porb
c
          else
c                 get data from partitioning file
c
              call iosys('rewind all on hpart read-and-write',0,0,0,' ')
              call iosys('read integer "target walks" '//
     1                   'from hpart',1,nwkst,0,' ')
              call iosys('read integer "target orbitals" from hpart',1,
     1                    trgorb,0,' ')
              call iosys('read integer nwksp from hpart',1,
     1                    twalks,0,' ')
              porb=twalks/nwkst
              call iosys('read integer nwksq from hpart',
     1                    1,qwalks,0,' ')
              call iosys('read integer "target roots" from hpart',1,
     1                    troots,0,' ')
              npvec=troots*porb
          endif
          write(iout,8) trgorb, porb, nwks, twalks, qwalks, 
     1                  troots, npvec
          if(twalks.ne.nwksp) then
             write(iout,9) twalks, nwksp
             call lnkerr('error in p space walks')
          endif
          if(qwalks.ne.nwksq) then
             write(iout,11) qwalks, nwksq
             call lnkerr('error in q space walks')
          endif
c
c             get the p space vectors

          ibuf=1
          rbuf=iadtwp(ibuf+2*lnbuf)
          diag=rbuf+lnbuf
          call iosys('read integer mxcore from rwf',1,canget,0,' ')          
          pvec=diag+nwksp
          trpvec=pvec+max(troots*porb*nwksp,nwks)
          scr=trpvec+troots*porb*max(nwksp,nwksq)
          ntot=troots*porb
          need=wpadti(scr+ntot*nwksp)
          if(need.gt.canget) then
             call lnkerr('cannot get enough memory for tranpp')
          endif
          title(1)='"p space buffers"'
          title(2)='"number of p space elements"'
          title(3)='"p space diagonals"'
          call getmem(need,p,ngot,'tranpp',0)
c
c         transform the hpp to the contracted p-space
c
          call tranpp(a(ibuf),z(rbuf),z(diag),z(pvec),z(trpvec),z(scr),
     1                lnbuf,nwks,nwksp,ntot,title,prnt)
c
c         transform hpq to the contracted p-space
c
          title(1)='"pq space buffers"'
          title(2)='"number of pq space elements"'
          call tranpq(a(ibuf),z(rbuf),z(pvec),z(trpvec),lnbuf,
     1                nwksp,nwksq,ntot,title,prnt)
          if(contrq) then
             call getmem(-ngot,p,idum,'tranpp',idum)
             ibuf=1
             rbuf=iadtwp(ibuf+2*lnbuf)
             qvec=rbuf+lnbuf
             trqvec=qvec+nwksq*qroots
             need=wpadti(trqvec+ntot*nwksq)
             if(need.gt.canget) then
                call lnkerr('cannot get enough getmem for '//
     1                      'final transformation')
             endif
             call getmem(need,p,ngot,'fnltrn',0)
             call fnltrn(a(ibuf),z(rbuf),z(qvec),z(trqvec),
     1                   lnbuf,nwksp,nwksq,qroots,ntot,prnt)
             call getmem(-ngot,p,idum,'fnltrn',idum)
          else
             call scr2ham(a(ibuf),z(rbuf),lnbuf)
             call getmem(-ngot,p,idum,'tranpp',idum)
          endif
      endif
      call iosys('close hamiltonian',0,0,0,' ')         
      return
 1    format(//,25x,'target calculation',/,5x,
     1            'no. target orbitals        = ',i7,/,5x,
     2            'no. L**2 p space orbitals  = ',i7,/,5x,
     3            'total no. configurations   = ',i7,/,5x,
     4            'no. p-space configurations = ',i7,/5x,
     5            'no. target roots           = ',i7)
 2    format (//,25x,'forming the p-space hamiltonian',/,5x,
     1               'total no. configurations   = ',i7,/,5x,
     2               'no. p-space configurations = ',i7,/,5x,
     3               'no. p-space roots          = ',i7)
 3    format (//,25x,'forming the q-space hamiltonian',/,5x,
     1               'total no. configurations   = ',i7,/,5x,
     2               'no. q-space configurations = ',i7,/,5x,
     3               'no. q-space roots          = ',i7)
 4    format (//,25x,'forming the p-q block of the hamiltonian',
     1               'total no. configurations   = ',i7,/,5x, 
     2               'no. p-space configurations = ',i7,/,5x,
     3               'no. q-space configurations = ',i7)
 5    format (//,25x,'forming partitioned  hamiltonian',/,5x,
     1               'total no. configurations   = ',i7,/,5x, 
     2               'no. p-space configurations = ',i7,/,5x,
     3               'no. p-space roots          = ',i7,/,5x,
     4               'no. q-space configurations = ',i7,/,5x,
     5               'no. q-space roots          = ',i7)
 6    format (//,25x,'forming p+q hamiltonian',/,5x,
     1               'total no. configurations   = ',i7,/,5x,
     2               'no. p-space configurations = ',i7,/,5x,
     3               'no. q-space configurations = ',i7,/,5x,
     4               'no. roots                  = ',i7)
 7    format (//,25x,'forming hamiltonian',/,5x,
     1               'no. configurations = ',i7,/,5x,
     2               'no. roots          = ',i7)
 8    format(//,25x,'transforming hamiltonian to contracted space',/,5x,
     1            'no. target orbitals       = ',i7,/,5x,
     2            'no. L**2 p space orbitals = ',i7,/,5x,
     3            'no. walks                 = ',i7,/,5x,
     4            'no. target walks          = ',i7,/,5x,
     5            'no. q space walks         = ',i7,/,5x,
     6            'no. target roots          = ',i7,/,5x,
     7            'size of p space           = ',i7)
 9    format(/,25x,'target walks = ',i7,1x,'p space walks = ',i7)
 11   format(/,25x,'q walks = ',i7,1x,'q space walks = ',i7)
      end











