c*deck @(#)pm806.f	1.2  7/30/91
      program m806
c
c***begin prologue     pm806
c***date written       yymmdd   (yymmdd)
c***revision date      871014   (yymmdd)
c
c 15 november 2004     bis at nsf
c      some significant changes in the way drt's are generated and
c      stored.  now more than one drt can be generated in a single
c      run.  each one is stored on a separate physical file.
c      the filenames are generated from the input drt 
c      information and their names put on the chk file.
c
c 14 october 1987      pws at lanl
c      merging m804 and m802 together so that m806 can generate
c      drt's for intermediate states for full ci's
c
c  8 december 1986   pws at lanl
c      modifying iosys open call to character unit name.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)pm806.f	1.2   7/30/91
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit integer (a-z)
c
c
      character*1 multrf,valenc,c(1000)
      character*1 itoc
      character*1 setno
      character*3 symchr
      character*3 codes,words*18
      character*16 cjunk
      character*4096 ops
      character*8 scatyp
      integer numel(50),symtp(50),spin(50)
      integer z
      pointer(p,z(1))
      logical qspace,pandq,symchk
      character*8 inpflg
      character*16 chrkey
      character*16 key
      logical logkey
      character*80 line
      character*128 tmp, nmcnfg, namchk
      character*8 dsk
      logical posinp
      logical noinp
c
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
      common /code/   fzc, fzv, cor, vir, doc, uoc, alp, bet, spe
      common /tapes/  out,input,drttap
      common /io/     inp,ioutq
      call drum
c
c
c     ----- recover the options string -----
c
c      call manmem(0,idum,idum,'m806',idum)
      write(ioutq,*) ' M806:Generate Direct Row Table for '
      write(ioutq,*) '      Scattering Codes'
      call iosys('read character options from rwf',-1,0,0,ops)
      if( logkey(ops,'drt=reuse',.false.,' ') ) then
         write(ioutq,*) 'skipping drt generation'
         call chainx(0)
         stop
      endif 
      call getmem(0,p,ngot,'first',0)
      call iosys('read integer mxcore from rwf',1,top,0,' ')
      call iosys('read integer mncore from rwf',1,bottom,0,' ')
      write(ioutq,*) 'minimum core available = ',bottom
      write(ioutq,*) 'maximum core available = ',top
      scatyp=chrkey(ops,'scattering','none',' ')
      if(chrkey(ops,'core-use=m806','maximum',' ').eq.'minimum') then
          top=bottom
      endif      
      write(ioutq,*) 'using ', top, ' words of memory'
      call getmem(top,p,maxcor,'m806',0)
      call getenv('MESA_TMP',tmp)
      ltmp=cskipb(tmp,' ')
      dsets=intkey(ops,'drt=sets',1,' ')
c
c             open the checkpoint file
c
      call iosys('read character "checkpoint filename" from rwf',
     $            0,0,0,namchk)
      call iosys('open chk as old',0,0,0,namchk)
c
c     ----- transfer global unit names to local ones -----
c
      input=inp
      out=ioutq
c
c
c
      do 5000 sets=1,dsets
c
         setno=itoc(sets)         
c        check for input flag replacement.
c
c        from the input generate the unique file name.
c
         key=chrkey(ops,'drt=key-'//setno,'drt-'//setno,' ')
         lkey=length(key)         
         nmcnfg=tmp(1:ltmp)//'/'//key(1:lkey)
         inpflg='$'//key(1:lkey)
         lflg=length(inpflg)
         write(out,5001) key(1:lkey)
 5001    format(/,15x,'generating configurations for set = ',a8)
c
c        read the options string.
c
         exv=intkey(ops,key(1:lkey)//'=vexcite',2,' ')
         exref=intkey(ops,key(1:lkey)//'=refexcite',0,' ')
c        out=chrkey(ops,key//'=out','out',' ')
         blksiz=intkey(ops,key(1:lkey)//'=blksiz',1000000,' ')
         absmax=blksiz
         levfrm=intkey(ops,key(1:lkey)//'=fermilevel',0,' ')
         nsym=intkey(ops,key(1:lkey)//'=nsym',0,' ')
         symchk=logkey(ops,key(1:lkey)//'=symmetry=rwf',
     1                 0,' ')
         nbf=intkey(ops,key(1:lkey)//'=nbf',0,' ')
         nrefs=intkey(ops,key(1:lkey)//'=nrefs',1,' ')
         ngroup=intkey(ops,key(1:lkey)//'=ngroups',0,' ')
c
c
         noinp=.not.posinp(inpflg(1:lflg),cjunk)
c
c
c     ----- see if doing p and q space drt's -----
c
         pandq=.false.
         dsk=key(1:lkey)
c
c        open the file for this drt set.
c
         call iosys('open '//dsk//' as new',0,0,0,nmcnfg)
         if(key(1:lkey).eq.'p-and-q') then
            pandq=.true.
         endif
         write(out,5002) dsk
c
c        write the names to the chk file so they can be used
c        later.  this is needed to pick off the correct file.
c        in other links.
c
         call iosys('write character "drt unit file name '//dsk
     #               //'" to chk',0,0,0,nmcnfg)
         call iosys('write character "drt file name '
     1              //key(1:lkey)//'" to chk',0,0,0,dsk)
 5002    format(/,15x,'configuration information put on file = ',a8)
c
c     ----- fill in defaults if needed -----
c
         if (nsym.eq.0) nsym=1
         if (nrefs.eq.0) nrefs=1
         if (nbf.eq.0) call iosys('read integer "number of basis '//
     $                            'functions" from rwf',1,nbf,0,' ')
         call iosys('write integer "number of drt functions" to '//
     1                 dsk,1,nbf,0,' ')
         call iosys('write integer "symmetries in ci" to rwf',
     1               1,nsym,0,' ')
         call iosys('write integer "symmetries in ci" to '//dsk,
     1               1,nsym,0,' ')
c
         write (out,2) nsym,nbf,nrefs,ngroup
    2    format (/,1x,'nsym = ',i3,1x,'nbf = ',i3,1x,'nrefs = ',i3,
     1             1x,'ngroup = ',i3)
c
c     ----- can only do p and q space from groups -----
c
         if (pandq.and.ngroup.lt.1) then
             call lnkerr('can only handle p-and-q space with '//
     1                   'group input')
         end if
c
         bfkey=1
         ntype=1
         bfnum=ntype+ntypes
         numsym=bfnum+ntypes*nbf
         bfsym=numsym+ntypes*nsym
         bfcode=bfsym+nbf
         if(symchk) then
            mosym=bfcode+nbf*nrefs
            molst=mosym+nsym  
            top=molst+nbf*nsym
            call mosin(z(bfsym),z(mosym),z(molst),nsym,nbf)
         else
            top=bfcode+nbf*nrefs
         endif
c
c     ----- read in the orbital codes -----
c
         na=0
         nb=0
         sym=0
c
         call scan(nbf,fzc,fzv,spe,uoc,cor,vir,alp,bet,occupd,
     #             speshl,virtul,frozen,rescor,resvir,opensh,
     #             valocc,valvir,multi,valence,multrf,words,
     #             norbs,ntypes,nsym,nrefs,inp,out,
     #             z(ntype),z(bfnum),z(numsym),z(bfsym),c(bfkey),
     #             z(bfcode),symchk)
c
c     ----- reorder the orbitals from the scf to the drt order -----
c
         iout=top
         orbsym=iout+nbf
         orbtbf=orbsym+norbs
         top=orbtbf+norbs
c
c     --------------------------------------------------------------
c           this just lets the code know that we are doing a
c           scattering calculation and puts the relevant data
c           in the correct palce for later useage.
c 
         if(scatyp(1:4).eq.'kohn') then
            write(out,1111)
1111        format(2x,' kohn-drt orbital order in use')
            scat=1
         else if(scatyp(1:8).eq.'r-matrix') then
            write(out,1112)
1112        format(2x,' r-matrix-drt orbital order in use')
            scat=1
         else
            scat=0
         end if
c
         call reordr(ntypes,nsym,norbs,nbf,speshl,opensh,multi,
     #               valocc,valvir,virtul,frozen,fzc,fzv,ncodes,
     #               z(ntype),z(bfnum),z(bfsym),z(bfcode),dela,delb,
     #               delele,na,nb,sym,levfrm,levocc,levopn,levmul,
     #               levval,orbfrm,z(iout),z(orbsym),z(orbtbf),scat)
c
c     ----- prepare the drt for the first (or only) reference -----
c
         nlevs=norbs+1
         orbcod=top
         start=orbcod+norbs*nrefs
c
         call movref(z(bfcode),nbf,nrefs,z(orbcod),norbs,z(iout))
c
c     ----- if doing groups, read in the a,b and symmetry of graph
c              head, and then the group occupancies, symmetries
c              and spins
c
         if (ngroup.ge.1) then
c            read (inp,nabs)
             read (inp,845) line
 845         format (a80)
             call locase(line,line)
             junk=na
             na=intkey(line,'na=',junk,' ')
             junk=nb
             nb=intkey(line,'nb=',junk,' ')
             junk=sym+1
             sym=intkey(line,'ns=',junk,' ')-1
         end if
c
c     ----- if doing p-and-q space , loop over p then q spaces -----
c
         if (pandq) then
            npq=2
         else
            npq=1
         end if
c
         pass=0
 2000    continue
c
c        ----- if on our second pass, then doing the intermediate
c              states. so, bump up excitation level, and do all
c              symmetries
c
         if (pass.eq.1) then
             exv=exv+1
             mxpass=nsym-1
             ninter=0
         else
             mxpass=0
         end if
         do 1500 tsym=0,mxpass
c
c           ----- if intermediate states, set symmetry to appropriate
c                 symmetry
c
            if (pass.eq.1) then
                ns=tsym
            else
                ns=sym
            end if
c
            do 1000 pq=1,npq
               noinp=.not.posinp(inpflg(1:lflg)//'-groups',cjunk)
               levpt1=start
               levnr1=levpt1+nlevs
               a1=levnr1+nlevs
               nrowmx=(maxcor-a1)/9
               b1=a1+nrowmx
               s1=b1+nrowmx
               x1=s1+nrowmx
               arc1=x1+nrowmx
               nlwks1=arc1+4*nrowmx
               botcor=nlwks1+nrowmx
               if (botcor.gt.maxcor) call lnkerr('core size')
c
               minref=2
               if (ngroup.ge.1) then
    3              continue
                   do 13 i=1,ngroup
                      numel(i)=-99999
                      symtp(i)=-1
                      spin(i)=-1
   13              continue
                   qspace=.false.
c                  read (inp,group)
                   read (inp,845) line
                   call locase(line,line)
                   call intarr(line,'numel',numel,ngroup,' ')
                   call intarr(line,'symtp',symtp,ngroup,' ')
                   call intarr(line,'spin',spin,ngroup,' ')
                   qspace=logkey(line,'qspace',.false.,' ')
c
c              ----- check for appropriate space ----
c
                   if (pandq) then
                       if (pq.eq.1.and.qspace) then
                           minref=minref+1
                           go to 3
                       else if (pq.eq.2.and..not.qspace) then
                           minref=minref+1
                           go to 3
                       end if
                   end if
c
                   do 14 i=2,ngroup
                      numel(i)=numel(i-1)+numel(i)
   14              continue
               end if
c
               call paldus(nlevs,na,nb,ns,nrowmx,norbs,ncodes,cor,vir,
     #                     levfrm,exv,nrows,nwks,out,
     #                     z(levpt1),z(levnr1),z(a1),z(b1),z(s1),z(x1),
     #                     z(orbcod),delele,z(orbsym),z(arc1),
     #                     z(nlwks1),numel,symtp,spin,ntypes)
               write(out,3000) nrows
 3000          format(/,25x,'[no. paldus rows = ', i5,']')          
c
               one=1
               write (out,4) ns+1,one,nwks,nwks,nwks,qspace
    4          format (/,t10,'symmetry ',i1,
     #         /,t10,'reference   configurations    unique      total'
     #               '    q',/,t10,'---------   --------------   '
     #                             '------      -----    -',/,t10,i5,
     #           i15,i14,i10,l8)
      
c
c     ----- collect the arrays together to free space -----
c
               levpt=start
               levnr=levpt+nlevs
               a=levnr+nlevs
               b=a+nrows
               s=b+nrows
               x=s+nrows
               arc=x+nrows
               nlwks=arc+4*nrows
               top=nlwks+nrows
c
               call movdrt(z(levpt),z(levpt1),z(levnr),z(levnr1),
     #                     z(a),z(a1),z(b),z(b1),z(s),z(s1),z(x),
     #                     z(x1),z(arc),z(arc1),z(nlwks),z(nlwks1),
     #                     nlevs,nrows)
c
c     ----- form drt's for other references -----
c
               do 10 ref=minref,nrefs
                  nwkold=nwks
                  levpt1=top
                  levnr1=levpt1+nlevs
                  a1=levnr1+nlevs
                  nrowmx=(maxcor-a1)/9
                  b1=a1+nrowmx
                  s1=b1+nrowmx
                  x1=s1+nrowmx
                  arc1=x1+nrowmx
                  nlwks1=arc1+4*nrowmx
                  top=nlwks1+nrowmx
                  if (top.gt.maxcor) call lnkerr('core size')
c
c     ----- if doing groups, read the group occupancies, symmetries
c              and spins
c
                  if (ngroup.ge.1) then
                      do 15 i=1,ngroup
                         numel(i)=-99999
                         symtp(i)=-1
                         spin(i)=-1
   15                 continue
                      qspace=.false.
c                     read (inp,group)
                      read (inp,845) line
                      call locase(line,line)
                      call intarr(line,'numel',numel,ngroup,' ')
                      call intarr(line,'symtp',symtp,ngroup,' ')
                      call intarr(line,'spin',spin,ngroup,' ')
                      qspace=logkey(line,'qspace',.false.,' ')
                      if (pandq) then
                          if (pq.eq.1.and.qspace) go to 9
                              if (pq.eq.2.and..not.qspace) go to 9
                      end if
                      do 16 i=2,ngroup
                         numel(i)=numel(i-1)+numel(i)
   16                 continue
                      orbcd=orbcod
                  else
                      orbcd=orbcod+(ref-1)*norbs
                  end if
c
                  call paldus(nlevs,na,nb,ns,nrowmx,norbs,ncodes,cor,
     #                        vir,levfrm,exv,nrows1,nwks1,out,
     #                        z(levpt1),z(levnr1),z(a1),z(b1),z(s1),
     #                        z(x1),z(orbcd),delele,z(orbsym),z(arc1),
     #                        z(nlwks1),numel,symtp,spin,ntypes)
c
c     ----- collect the arrays together to free space -----
c
                  nrows2=nrows1
                  levpt2=levpt1
                  levnr2=levpt2+nlevs
                  a2=levnr2+nlevs
                  b2=a2+nrows2
                  s2=b2+nrows2
                  x2=s2+nrows2
                  arc2=x2+nrows2
                  nlwks2=arc2+4*nrows2
                  top=nlwks2+nrows2
c
                  call movdrt(z(levpt2),z(levpt1),z(levnr2),
     #                        z(levnr1),z(a2),z(a1),z(b2),z(b1),
     #                        z(s2),z(s1),z(x2),z(x1),z(arc2),
     #                        z(arc1),z(nlwks2),z(nlwks1),nlevs,nrows2)
c
c        ----- merge the two drt's together -----
c
                  levpt1=top
                  levnr1=levpt1+nlevs
                  a1=levnr1+nlevs
                  nrowmx=(maxcor-a1)/10
                  b1=a1+nrowmx
                  s1=b1+nrowmx
                  x1=s1+nrowmx
                  arc1=x1+nrowmx
                  nlwks1=arc1+4*nrowmx
                  torow1=nlwks1
                  torow2=torow1+nrowmx
                  top=torow2+nrowmx
                  if (top.gt.maxcor) call lnkerr('core size')
c
                  call mrgdrt(out,nlevs,nrowmx,nrows1,nrows,nrows2,
     #                        z(a1),z(b1),z(s1),z(x1),z(arc1),
     #                        z(levpt1),z(levnr1),z(a),z(b),z(s),
     #                        z(x),z(arc),z(levpt),z(levnr),
     #                        z(a2),z(b2),z(s2),z(x2),z(arc2),
     #                        z(levpt2),z(levnr2),z(torow1),
     #                        z(torow2),.false.)
c
                  nrows2=nrows1
                  levpt2=start
                  levnr2=levpt2+nlevs
                  a2=levnr2+nlevs
                  b2=a2+nrows2
                  s2=b2+nrows2
                  x2=s2+nrows2
                  arc2=x2+nrows2
                  nlwks2=arc2+4*nrows2
                  top=nlwks2+nrows2
c
                  call movdrt(z(levpt2),z(levpt1),z(levnr2),z(levnr1),
     #                        z(a2),z(a1),z(b2),z(b1),z(s2),z(s1),
     #                        z(x2),z(x1),z(arc2),z(arc1),z(nlwks2),
     #                        z(nlwks1),nlevs,nrows2)
c
                  ntemp=100
                  temp=top
                  mxsim=7
                  levpt1=temp+ntemp
                  levnr1=levpt1+nlevs
                  a1=levnr1+nlevs
                  nrowmx=(maxcor-a1)/(9+mxsim)
                  b1=a1+nrowmx
                  s1=b1+nrowmx
                  x1=s1+nrowmx
                  ptsim=x1
                  arc1=x1+nrowmx
                  nlwks1=arc1+4*nrowmx
                  nsim=nlwks1
                  sim=nsim+nrowmx
                  top=sim+mxsim*nrowmx
                  if (top.gt.maxcor) call lnkerr('core size')
c
                  call squash(z(a1),z(b1),z(s1),z(arc1),z(nlwks1),
     #                        z(ptsim),z(nsim),z(sim),z(levpt1),
     #                        z(levnr1),z(a2),z(b2),z(s2),z(arc2),
     #                        z(levpt2),z(levnr2),z(temp),nrows,
     #                        nrows2,nrowmx,nlevs,nwks,out,
     #                        mxsim*nrowmx,ntemp)
c
c        ----- and finally, move the drt arrays down to open space -----
c
                  levpt=start
                  levnr=levpt+nlevs
                  a=levnr+nlevs
                  b=a+nrows
                  s=b+nrows
                  x=s+nrows
                  arc=x+nrows
                  nlwks=arc+4*nrows
                  top=nlwks+nrows
c
                  call movdrt(z(levpt),z(levpt1),z(levnr),z(levnr1),
     #                        z(a),z(a1),z(b),z(b1),z(s),z(s1),z(x),
     #                        z(x1),z(arc),z(arc1),z(nlwks),z(nlwks1),
     #                        nlevs,nrows)
c
                  write (out,5) ref,nwks1,nwks-nwkold,nwks,qspace
    5             format (t10,i5,i15,i14,i10,l8)
    9             continue
c
c         ----- reset top for those cases where we skipped a p or q
c               configuration list
c
                  top=nlwks+nrows
   10          continue
c
c     ----- if p-and-q space drt's and have done p, save it on a scratch
c           file
c
               if (pandq.and.pq.eq.1) then
                   call iosys('open pdrt as scratch',0,0,0,' ')
                   call iosys('write integer "nrows" to pdrt',1,
     $                         nrows,0,' ')
                   call iosys('write integer levpt to pdrt',nlevs,
     $                         z(levpt),0,' ')
                   call iosys('write integer levnr to pdrt',nlevs,
     $                         z(levnr),0,' ')
                   call iosys('write integer a to pdrt',nrows,
     $                         z(a),0,' ')
                   call iosys('write integer b to pdrt',nrows,
     $                         z(b),0,' ')
                   call iosys('write integer x to pdrt',nrows,
     $                         z(x),0,' ')
                   call iosys('write integer s to pdrt',nrows,
     $                         z(s),0,' ')
                   call iosys('write integer arc to pdrt',4*nrows,
     $                         z(arc),0,' ')
                   call iosys('write integer nlwks to pdrt',nrows,
     $                         z(nlwks),0,' ')
                   nwksp=nwks
                   nrowsp=nrows
c
c           ----- reposition input -----
c
                   rewind inp
                   noinp=.not.posinp(inpflg(1:lflg),cjunk)
               else
                   nwksq=nwks
               end if
 1000       continue
            if (pandq) then
c
c         ----- read p drt back in and merge the two drt's -----
c
                plevpt=top
                plevnr=plevpt+nlevs
                pa=plevnr+nlevs
                pb=pa+nrowsp
                ps=pb+nrowsp
                px=ps+nrowsp
                parc=px+nrowsp
                pnlwks=parc+4*nrowsp
c
c        ----- and space for the merged drt -----
c
                levpt1=pnlwks+4*nrowsp
                levnr1=levpt1+nlevs
                a1=levnr1+nlevs
                nrowsq=nrows
                nrows=nrowsp+nrowsq
                b1=a1+nrows
                s1=b1+nrows
                x1=s1+nrows
                arc1=x1+nrows
                nlwks1=arc1+4*nrows
                torow1=nlwks1+4*nrows
                torow2=torow1+nrows
                arcsv=torow2+nrows
                prowsv=arcsv+nlevs
                qrowsv=prowsv+nlevs
                top=qrowsv+nlevs
c
                if (top.gt.maxcor) call lnkerr('not enough core for pq')
c
                call iosys('read integer levpt from pdrt',nlevs,
     $                      z(plevpt),0,' ')
                call iosys('read integer levnr from pdrt',nlevs,
     $                      z(plevnr),0,' ')
                call iosys('read integer a from pdrt',nrowsp,
     $                      z(pa),0,' ')
                call iosys('read integer b from pdrt',nrowsp,
     $                      z(pb),0,' ')
                call iosys('read integer x from pdrt',nrowsp,
     $                      z(px),0,' ')
                call iosys('read integer s from pdrt',nrowsp,
     $                      z(ps),0,' ')
                call iosys('read integer arc from pdrt',4*nrowsp,
     $                      z(parc),0,' ')
                call iosys('read integer nlwks from pdrt',nrowsp,
     $                      z(pnlwks),0,' ')
c
c        ----- set up x arrays to differentiate between p and q -----
c
                do 23 i=1,nrowsp
                   z(px+i-1)=1
   23           continue
                do 24 i=1,nrowsq
                   z(x+i-1)=0
   24           continue
c
                call mrgdrt(out,nlevs,nrows,nrows1,nrowsq,nrowsp,
     #                      z(a1),z(b1),z(s1),z(x1),z(arc1),
     #                      z(levpt1),z(levnr1),z(a),z(b),z(s),z(x),
     #                      z(arc),z(levpt),z(levnr),z(pa),z(pb),
     #                      z(ps),z(px),z(parc),z(plevpt),z(plevnr),
     #                      z(torow1),z(torow2),.true.)
c
                call pqx(nrows1,nlevs,z(arc1),z(x1),z(prowsv),
     #                   z(qrowsv),z(arcsv),z(levpt1))
            else
                a1=a
                b1=b
                s1=s
                arc1=arc
                levpt1=levpt
                levnr1=levnr
                nrows1=nrows
                nlwks1=nlwks
            end if
c
c     ----- work out the number-of-lower-walks -----
c
            call lwrwks(z(a1),z(b1),z(s1),z(x1),z(nlwks1),z(arc1),
     #                  z(levpt1),z(levnr1),nlevs,nrows1,nwks,out)
c
c     ----------------------------------------------
c           begin to write the information to the needed file.
c    
            if (nwks.le.0) then
                if (pass.eq.0) call lnkerr('no walks !?!')
                symchr=itoc(ns)
                call iosys('write integer "'//symchr//' nwks" to '
     1                     //dsk,1,nwks,0,' ')
                call iosys('write integer "'//symchr//' nwks" to '//
     1                      dsk,1,nwks,0,' ')
                go to 1500
            end if
c
c     ----- and pack it down -----
c
            nrows=nrows1
            levpt=start
            levnr=levpt+nlevs
            a=levnr+nlevs
            b=a+nrows
            s=b+nrows
            x=s+nrows
            arc=x+nrows
            nlwks=arc+4*nrows
            top=nlwks+nrows
c
            call movdrt(z(levpt),z(levpt1),z(levnr),z(levnr1),z(a),
     #                  z(a1),z(b),z(b1),z(s),z(s1),z(x),z(x1),
     #                  z(arc),z(arc1),z(nlwks),z(nlwks1),nlevs,nrows)
c
c     ----- and sort the rows to a nice order -----
c
            levpt1=top
            levnr1=levpt1+nlevs
            a1=levnr1+nlevs
            b1=a1+nrows
            s1=b1+nrows
            x1=s1+nrows
            arc1=x1+nrows
            nlwks1=arc1+4*nrows
            top=nlwks1+nrows
c
            call srtdrt(z(levpt1),z(levpt),z(levnr1),z(levnr),
     #                  z(a1),z(a),z(b1),z(b),z(s1),z(s),z(x1),z(x),
     #                  z(arc1),z(arc),z(nlwks1),z(nlwks),nlevs,
     #                  nrows,nsym)
c
c
c     ----- and pack it down -----
c
            levpt=start
            levnr=levpt+nlevs
            a=levnr+nlevs
            b=a+nrows
            s=b+nrows
            x=s+nrows
            arc=x+nrows
            nlwks=arc+4*nrows
            top=nlwks+nrows
c
            call movdrt(z(levpt),z(levpt1),z(levnr),z(levnr1),z(a),
     #                  z(a1),z(b),z(b1),z(s),z(s1),z(x),z(x1),
     #                  z(arc),z(arc1),z(nlwks),z(nlwks1),nlevs,nrows)
c
c     ----- calculate the arc weights, i.e. the numbering scheme -----
c
            wt=top
            wtab=wt+4*nrows
            wtw=wtab+orbfrm
            wtx=wtw+orbfrm*nsym
            wty=wtx+orbfrm*nsym
            top=wty+orbfrm
c
            if (top.gt.maxcor) call lnkerr('drt: weight, not enough'//
     1                                     ' core')
c
            call weight(nlevs,nrows,levfrm,z(levpt),z(levnr),z(arc),
     #                  z(wt),z(x),z(nlwks))
c
            if (index(ops,'p-and-q').ne.0) then
                nwksq=z(nlwks)
                nwksp=z(nlwks+1)
            else
                nwksp=z(nlwks)
                nwksq=0
            end if
c
c    ------------ continue the writing ----------------------------
c
            call iosys('write integer nwksp to '//dsk,1,nwksp,0,' ')
            call iosys('write integer nwksq to '//dsk,1,nwksq,0,' ')
            call iosys('write integer nwksp to rwf',1,nwksp,0,' ')
            call iosys('write integer nwksq to rwf',1,nwksq,0,' ')
c
c
c
            call extern(orbfrm,nsym,z(a),z(b),z(s),z(levpt),z(levnr),
     #                  z(arc),z(wt),z(wtab),z(wtw),z(wtx),z(wty),
     #                  nrows,nlevs)
c
c     ----- the integral addressing scheme -----
c
            nnp=(norbs+1)*norbs/2
            nnpvir=(orbfrm+1)*orbfrm/2
            ijgrp=top
            ijadd=ijgrp+nnp
            kadd=ijadd+nnp
            ladd=kadd+norbs*nsym
            inint=ladd+norbs*nsym
            ijxx=inint+norbs
            ijww=ijxx+nnp
            klxx=ijww+nnp
            klww=klxx+nnpvir
            nklxx=klww+nnpvir
            nklww=nklxx+nsym*orbfrm
c           top=nklww+nsym*orbfrm
            inext=nklww+nsym*orbfrm
            jmnnxt=inext+norbs
            jmxnxt=jmnnxt+norbs
            ningrp=jmxnxt+norbs
            top=ningrp+nnp
c
            if (top.gt.maxcor) call lnkerr('drt: intgrl core problem')
            if (pass.eq.0) then
                call intgrl(z(bfsym),z(orbtbf),z(kadd),z(ladd),
     #                      z(ijadd),z(ijgrp),z(inint),z(inext),
     #                      z(jmnnxt),z(jmxnxt),z(ningrp),z(ijxx),
     #                      z(klxx),z(nklxx),z(ijww),z(klww),
     #                      z(nklww),nbf,norbs,nsym*norbs,nnp,
     #                      nsym,orbfrm,ioutq,ioutq,absmax,nmax,ngrp)
            end if
c
c     ----- print what is desired -----
c
            if (pass.eq.0.and.logkey(ops,'print=drt',.false.,' ')) then
                call prtdrt(z(a),z(b),z(s),z(x),z(nlwks),z(arc),
     #                      z(levpt),z(levnr),z(wt),nlevs,nrows,out,
     #                      z(orbsym),z(ijgrp),z(ijadd),z(kadd),
     #                      z(ladd),z(inint),z(ijxx),z(ijww),z(klxx),
     #                      z(klww),z(nklxx),z(nklww),nnpvir,orbfrm,
     #                      z(wtab),z(wtw),z(wtx),z(wty),nnp,nsym,norbs)
            end if
c
c     ----- put offsets into arc array -----
c
            call renumb(z(levpt),z(levnr),z(arc),nlevs,nrows)
c
c     ----- put information on the needed dsk file -----
c
            if (pass.eq.0) then
                call iosys('write integer "symmetries in ci" to '
     1                      //dsk,1,nsym,0,' ')
                call iosys('write integer "nrows" to '//dsk,1,
     1                      nrows,0,' ')
                call iosys('write integer "nlevs" to '//dsk,1,
     1                      nlevs,0,' ')
                call iosys('write integer "norbs" to '//dsk,1,
     1                      norbs,0,' ')
                call iosys('write integer "nrefs" to '//dsk,1,
     1                      nrefs,0,' ')
                call iosys('write integer "nwks" to '//dsk,1,
     1                      nwks,0,' ')
                call iosys('write integer "orbfrm" to '//dsk,1,
     1                      orbfrm,0,' ')
                call iosys('write integer "symorb" to '//dsk,1,
     1                      norbs*nsym,0,' ')
                call iosys('write integer "numij" to '//dsk,1,
     1                      nnp,0,' ')
                call iosys('write integer "ngroup" to '//dsk,1,
     1                      ngrp,0,' ')
                call iosys('write integer "nmax" to '//dsk,1,
     1                      nmax,0,' ')
                call iosys('write integer "maxb" to '//dsk,1,maxb,0,' ')
                call iosys('write integer "nijvir" to '//dsk,1,nnpvir,
     1                      0,' ')
c
                call iosys('write integer kadd to '//dsk,norbs*nsym,
     1                      z(kadd),0,' ')
                call iosys('write integer ladd to '//dsk,norbs*nsym,
     1                      z(ladd),0,' ')
                call iosys('write integer ijadd to '//dsk,nnp,
     1                      z(ijadd),0,' ')
                call iosys('write integer ijgrp to '//dsk,nnp,
     1                      z(ijgrp),0,' ')
                call iosys('write integer inint to '//dsk,norbs,
     1                      z(inint),0,' ')
                call iosys('write integer inext to '//dsk,norbs,
     1                      z(inext),0,' ')
                call iosys('write integer jmnnxt to '//dsk,norbs,
     1                      z(jmnnxt),0,' ')
                call iosys('write integer jmxnxt to '//dsk,norbs,
     1                      z(jmxnxt),0,' ')
                call iosys('write integer ningrp to '//dsk,ngrp,
     1                      z(ningrp),0,' ')
                call iosys('write integer iout to '//dsk,nbf,
     1                      z(iout),0,' ')
                call iosys('write integer orbsym to '//dsk,norbs,
     1                      z(orbsym),0,' ')
                call iosys('write integer bfsym to '//dsk,nbf,
     1                      z(bfsym),0,' ')
                call iosys('write integer bfcode to '//dsk,nrefs*nbf,
     1                      z(bfcode),0,' ')
                call iosys('write integer orbtbf to '//dsk,norbs,
     1                      z(orbtbf),0,' ')
                call iosys('write integer levpt to '//dsk,nlevs,
     1                      z(levpt),0,' ')
                call iosys('write integer levnr to '//dsk,nlevs,
     1                      z(levnr),0,' ')
                call iosys('write integer a to '//dsk,nrows,z(a),0,' ')
                call iosys('write integer b to '//dsk,nrows,z(b),0,' ')
                call iosys('write integer x to '//dsk,nrows,z(x),0,' ')
                call iosys('write integer s to '//dsk,nrows,z(s),0,' ')
                call iosys('write integer arc to '//dsk,4*nrows,
     1                      z(arc),0,' ')
                call iosys('write integer nlwks to '//dsk,nrows,
     1                      z(nlwks),0,' ')
                call iosys('write integer weight to '//dsk,4*nrows,
     1                      z(wt),0,' ')
                if (orbfrm.gt.0) then
                    call iosys('write integer wtab to '//dsk,orbfrm,
     1                          z(wtab),0,' ')
                    call iosys('write integer wtw to '//dsk,orbfrm*nsym,
     1                          z(wtw),0,' ')
                    call iosys('write integer wtx to '//dsk,orbfrm*nsym,
     1                          z(wtx),0,' ')
                    call iosys('write integer wty to '//dsk,orbfrm,
     1                          z(wty),0,' ')
                    call iosys('write integer ijxx to '//dsk,nnp,
     1                          z(ijxx),0,' ')
                    call iosys('write integer klxx to '//dsk,nnpvir,
     1                          z(klxx),0,' ')
                    call iosys('write integer nklxx to '//dsk,
     1                          nsym*orbfrm,z(nklxx),0,' ')
                    call iosys('write integer ijww to '//dsk,nnp,
     1                          z(ijww),0,' ')
                    call iosys('write integer klww to '//dsk,nnpvir,
     1                          z(klww),0,' ')
                    call iosys('write integer nklww to '//dsk,
     1                          nsym*orbfrm,z(nklww),0,' ')
                end if
            else
                symchr=itoc(ns)
                call iosys('write integer "'//symchr//' nwks" to '
     1                     //dsk,1,nwks,0,' ')
                call iosys('write integer "'//symchr//' nrows" to '
     1                     //dsk,1,nrows,0,' ')
                call iosys('write integer "'//symchr//' maxb" to '
     1                     //dsk,1,maxb,0,' ')
c
                call iosys('write integer "'//symchr//' levpt" to '
     1                     //dsk,nlevs,z(levpt),0,' ')
                call iosys('write integer "'//symchr//' levnr" to '
     1                     //dsk,nlevs,z(levnr),0,' ')
                call iosys('write integer "'//symchr//' a" to '//dsk,
     1                      nrows,z(a),0,' ')
                call iosys('write integer "'//symchr//' b" to '//dsk,
     1                      nrows,z(b),0,' ')
                call iosys('write integer "'//symchr//' x" to '//dsk,
     1                      nrows,z(x),0,' ')
                call iosys('write integer "'//symchr//' s" to '//dsk,
     1                      nrows,z(s),0,' ')
                call iosys('write integer "'//symchr//' arc" to '//dsk,
     1                      4*nrows,z(arc),0,' ')
                call iosys('write integer "'//symchr//' nlwks" to '
     1                     //dsk,nrows,z(nlwks),0,' ')
                call iosys('write integer "'//symchr//' weight" to '
     1                     //dsk,4*nrows,z(wt),0,' ')
                ninter=ninter+nwks
            end if
c
 1500    continue
         pass=pass+1
         if (pass.eq.1.and.(logkey(ops,'ci=full',.false.,' ').or.
     $       logkey(ops,'mcscf=cas',.false.,' '))) go to 2000
c
         if (pass.eq.2) then
             write (out,400) ninter
 400         format (t10,'number of intermediate walks:',i10)
         end if
         if(pandq) then
            write(ioutq,1600) nwksp, nwksq
            call iosys('write integer nwksp to chk',1,nwksp,0,' ')
            call iosys('write integer nwksq to chk',1,nwksq,0,' ')
         endif
c
c      ----------------------------------------------------------
c               close the dsk file and exit
c
         call iosys('close '//dsk,0,0,0,' ')
 5000 continue   
      call iosys('close chk',0,0,0,' ')
 1600 format(/,5x,'writing nwksp = ',i8,' and nwksq = ',i8,' to chk')
c
c     ----- and exit with grace -----
c
c
      call getmem(-maxcor,p,idum,'m806',idum)
      call chainx(0)
c
c
      stop
      end
