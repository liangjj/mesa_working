*deck @(#)mn903.f	5.1  11/6/94
      subroutine mn903(a,z,maxcor,oper,cvec,svec,tunit,tdunit,ttype,
     #                 mcroot,vn)
c***begin prologue     mn903
c***date written       850101  (yymmdd) 
c***revision date      910607  (yymmdd)     
c
c   7 june    1991    rlm at lanl
c     adding nroots-at-a-time option
c   5 june    1991    rlm at lanl
c     passing down cvec,svec to hc and tdm.
c  10 may     1990    rlm at lanl
c     dynamical core allocation/ 32bit-64bit conversions
c  30 march   1988    bhl at llnl
c     logkey test of print of h.times.c commented out
c  28 january 1988    bhl at brl
c     mcscf=ci=   options inserted for cycles and convergence
c
c***keywords           
c***author             saxe, paul 
c***source             @(#)mn903.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       mn903
c
      implicit integer (a-z)
      character*2 itoc
      character*(*) oper
      character*2 vn
      character*4096 ops
      character*8 prtflg
      character*(*) ttype
      character*(*) tunit
      character*(*) tdunit
      character*16 unit
      character*16 dunit
c
      integer a(maxcor)
      integer drtpt(8)
      integer numrow(8)
      integer numsym(8)
      integer offsym(8)
c
      logical logkey
c
      real*8 cnverg
      real*8 cvec(*)
      real*8 fpkey
      real*8 svec(*)
      real*8 thresh
      real*8 fpjunk
      real*8 z(*)
c
      common /io/ inp,iout
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
      unit=tunit
      dunit=tdunit
c
c     ----- decipher this links options -----
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      nroots=max(mcroot,1)
c
      if(oper.eq.'ci') then
c
         junk=nroots
         nroots=intkey(ops,'ci=nroots',junk,' ')
         nattim=intkey(ops,'ci=nroots-at-a-time',nroots,' ')
c
c        ----- for the mcscf, tighten up the default convergence -----
c
         if (logkey(ops,'mcscf',.false.,' ')) then
            fpjunk=1.0d-07
            cnverg=fpkey(ops,'mcscf=ci=convergence',fpjunk,' ')
            thresh=fpkey(ops,'mcscf=ci=threshold',cnverg,' ')
            junk=intkey(ops,'mcscf=ci=cycles',20,' ')
            mxiter=intkey(ops,'mcscf=ci=iterations',junk,' ')
         else
            fpjunk=1.0d-05
            cnverg=fpkey(ops,'ci=convergence',fpjunk,' ')
            thresh=fpkey(ops,'ci=threshold',cnverg,' ')
            junk=intkey(ops,'ci=cycles',20,' ')
            mxiter=intkey(ops,'ci=iterations',junk,' ')
         end if
c
         nguess=intkey(ops,'ci=nguess',nroots,' ')
         maxwks=intkey(ops,'ci=intermediates',-1,' ')
         mxiter=mxiter*nattim+nguess
c
      end if
c
c     ----- get sizes for dividing core ------
c
      call iosys('read integer norbs from rwf',1,norbs,0,' ')
      call iosys('read integer nlevs from rwf',1,nlevs,0,' ')
      call iosys('read integer nwks  from rwf',1,nwks, 0,' ')
      call iosys('read integer nrows from rwf',1,nrows,0,' ')
      call iosys('read integer "symmetries in ci" from rwf',
     #            1,nsym,0,' ')
c
c     ----- read in dimensions of intermediate drt's for each symmetry
c
      inrows=0
      do 1 sym=1,nsym
         drtpt(sym)=inrows+1
         call iosys('read integer "'//itoc(sym-1)//'nwks" from rwf',
     $               1,nwksym,0,' ')
         if (nwksym.eq.0) then
            numrow(sym)=0
         else
            call iosys('read integer "'//itoc(sym-1)//
     $                 'nrows" from rwf',1,numrow(sym),0,' ')
         end if
         inrows=inrows+numrow(sym)
    1 continue
      nnp=(norbs+1)*norbs/2
c
c     ----- integral addressing array and size of integral list ----
c
      ijpt=1
      orbsym=ijpt+nnp
      top=orbsym+norbs
c
      call getscm(top,a,maxcor,'m903',0)
c
      call iosys('read integer orbsym from rwf',norbs,a(orbsym),
     $            0,' ')
c
      call intad(a(orbsym),norbs,nsym,numsym,offsym,a(ijpt),ngsym,nnp)
      nnp0=numsym(1)
c
c     ----- allocate core -----
c
      dnarc=top
      uparc=dnarc+4*nrows
      dnnwks=uparc+4*nrows
      upnwks=dnnwks+nrows
      dnwt=upnwks+nrows
      upwt=dnwt+4*nrows
      levpt=upwt+4*nrows
      levnr=levpt+nlevs
      rowoff=levnr+nlevs
      aval=rowoff+nrows
      bval=aval+nrows
      sval=bval+nrows
      idnarc=sval+nrows
      iuparc=idnarc+4*inrows
      idnnwk=iuparc+4*inrows
      iupnwk=idnnwk+inrows
      idnwt=iupnwk+inrows
      iupwt=idnwt+4*inrows
      ilevpt=iupwt+4*inrows
      ilevnr=ilevpt+nlevs*nsym
      irowof=ilevnr+nlevs*nsym
      iaval=irowof+inrows
      ibval=iaval+inrows
      isval=ibval+inrows
      irowsv=isval+inrows
      jrowsv=irowsv+nlevs
      iwtsv=jrowsv+nlevs
      jwtsv=iwtsv+nlevs
      pagesv=jwtsv+nlevs
      segsv=pagesv+nlevs
      refwlk=segsv+nlevs
c
      acoef=iadtwp(refwlk+nguess)
      bcoef=acoef+nlevs
      h=bcoef+nlevs
      g=h+nnp0
      if (oper.eq.'h times c') then
         root=g+ngsym
         dvec=root
         dmat=dvec
         top1=wpadti(dvec)
         h1=iadtwp(top1)
         g1=h1+nnp
         top2=wpadti(g1+nnp**2)
      else if (oper.eq.'ci') then
         mxvc=min(64,nroots)
         if(prtflg.ne.'minimum') then
            write(iout,700) mxvc,nroots
  700       format(5x,'davidson vectors in use:      ',i8,
     $            /5x,'davidson vectors possible:    ',i8)
         end if
         root=g+ngsym
         dvec=root+mxiter
         dmat=dvec+mxiter**2
         c=dmat+((mxiter+1)*mxiter/2)
         s=c+mxvc*nwks
         top1=wpadti(s+mxvc*nwks)
         h1=iadtwp(top1)
         g1=h1+nnp
         top2=wpadti(g1+nnp**2)
      else if (oper.eq.'density') then
         c=g+ngsym
         h1=c+nwks
         g1=h1+nnp
         top1=wpadti(g1+nnp**2)
         top2=wpadti(g1+nnp**2)
      else if (oper.eq.'transition density') then
         h1=g+ngsym
         g1=h1+nnp
         top1=wpadti(g1+nnp**2)
         top2=wpadti(g1+nnp**2)
      end if
c
c
c..bhl
      top=max(top1,top2)+250000
c     ----- get core needed -----
      call getscm(top,z,maxcor,'need for m903',0)
c
      if (oper.eq.'ci'.or.oper.eq.'h times c') then
c
c        ----- open integrals unit and get them -----
c
         if (ttype.eq.'ci') then
            call iosys('read real "mo one-electron integrals" from '//
     $                  unit,nnp,z(h1),0,' ')
            call iosys('read real "mo two-electron integrals" from '//
     $                  unit,nnp**2,z(g1),0,' ')
         else if (ttype.eq.'mcscf') then
            call iosys('read real "h mcscf" from '//unit,
     $                 nnp,z(h1),0,' ')
            call iosys('read real "g mcscf" from '//unit,
     $                 nnp**2,z(g1),0,' ')
         else
            call lnkerr('bad type passed to mn903')
         end if
      end if
c
c     ----- see if a dividing level was entered as an option -----
c
      if(oper.eq.'ci') then
         levdiv=intkey(ops,'ci=dividing-level',(nlevs+1)/2,' ')
      else
         levdiv=(nlevs+1)/2
      end if
c
c
      if (oper.eq.'h times c') then
c
c     ----- multiply h.c
c
c        trial vector passed down in cvec; product returned in svec.
         call hc(a(uparc),a(dnarc),a(upwt),a(dnwt),a(upnwks),a(dnnwks),
     #        a(aval),a(bval),a(levpt),a(levnr),a(rowoff),
     #        a(iuparc),a(idnarc),a(iupwt),a(idnwt),a(iupnwk),
     #        a(idnnwk),a(iaval),a(ibval),a(ilevpt),a(ilevnr),
     #        inrows,nrows,nlevs,norbs,
     #        z(acoef),a(irowsv),a(jrowsv),a(iwtsv),a(jwtsv),
     #        a(pagesv),a(segsv),
     #        levdiv,nnp,cvec,nwks,svec,z(h),z(g),z(bcoef),prtflg,
     #        nsym,numrow,drtpt,a(top1),a(sval),a(isval),
     #        maxwks,a(ijpt),nnp0,ngsym,numsym,offsym,a(orbsym),z(h1),
     $        z(g1),ops)
      else if (oper.eq.'ci') then
c
c     ----- do the ci calculation -----
c
         froot=max(1,mcroot)
         call ci(a(uparc),a(dnarc),a(upwt),a(dnwt),a(upnwks),a(dnnwks),
     #        a(aval),a(bval),a(levpt),a(levnr),a(rowoff),
     #        a(iuparc),a(idnarc),a(iupwt),a(idnwt),a(iupnwk),
     #        a(idnnwk),a(iaval),a(ibval),a(ilevpt),a(ilevnr),
     #        inrows,nrows,nlevs,norbs,
     #        z(acoef),a(irowsv),a(jrowsv),a(iwtsv),a(jwtsv),
     #        a(pagesv),a(segsv),
     #        levdiv,nnp,z(c),nwks,z(s),z(h),z(g),
     #        z(bcoef),z(root),z(dvec),z(dmat),nroots,mxiter,
     #        cnverg,thresh,nguess,prtflg,a(refwlk),
     #        nsym,numrow,drtpt,a(top1),a(sval),a(isval),
     #        maxwks,a(ijpt),nnp0,ngsym,numsym,offsym,a(orbsym),z(h1),
     $        z(g1),ops,nattim,mxvc,froot,vn)
      else if (oper.eq.'density') then
c
c     ----- form the density matrix -----
c
c        get the appropriate vector.
         froot=max(1,mcroot)
         call iosys('read real "'//vn//' root '//itoc(froot)
     $               //'" from rwf',nwks,z(c),0,' ')
         call dm(a(uparc),a(dnarc),a(upwt),a(dnwt),a(upnwks),a(dnnwks),
     #        a(aval),a(bval),a(levpt),a(levnr),a(rowoff),
     #        a(iuparc),a(idnarc),a(iupwt),a(idnwt),a(iupnwk),
     #        a(idnnwk),a(iaval),a(ibval),a(ilevpt),a(ilevnr),
     #        inrows,nrows,nlevs,norbs,
     #        z(acoef),a(irowsv),a(jrowsv),a(iwtsv),a(jwtsv),
     #        a(pagesv),a(segsv),
     #        levdiv,nnp,z(c),nwks,z(h),z(g),
     #        z(bcoef),prtflg,
     #        nsym,numrow,drtpt,a(top1),a(sval),a(isval),
     #        maxwks,a(ijpt),nnp0,ngsym,numsym,offsym,a(orbsym),z(h1),
     $        z(g1),ops,dunit)
      else if (oper.eq.'transition density') then
c
c     ----- form the transition density matrix -----
c
c        vectors should come down in cvec,svec
         call tdm(a(uparc),a(dnarc),a(upwt),a(dnwt),a(upnwks),
     $        a(dnnwks),
     #        a(aval),a(bval),a(levpt),a(levnr),a(rowoff),
     #        a(iuparc),a(idnarc),a(iupwt),a(idnwt),a(iupnwk),
     #        a(idnnwk),a(iaval),a(ibval),a(ilevpt),a(ilevnr),
     #        inrows,nrows,nlevs,norbs,
     #        z(acoef),a(irowsv),a(jrowsv),a(iwtsv),a(jwtsv),
     #        a(pagesv),a(segsv),
     #        levdiv,nnp,cvec,nwks,svec,z(h),z(g),
     #        z(bcoef),prtflg,
     #        nsym,numrow,drtpt,a(top1),a(sval),a(isval),
     #        maxwks,a(ijpt),nnp0,ngsym,numsym,offsym,a(orbsym),z(h1),
     $        z(g1),ops,dunit)
      end if
c
c
      return
      end
