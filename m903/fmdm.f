*deck @(#)fmdm.f	5.1  11/6/94
      subroutine fmdm(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,a,b,sval,
     #                  levpt,levnr,rowoff,
     #                  iuparc,idnarc,iupwt,idnwt,iupnwk,idnnwk,
     #                  ia,ib,isval,ilevpt,ilevnr,inrows,
     #                  nrows,nlevs,norbs,
     #                  acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,segsv,
     #                  levdiv,nnp,c,nwks,h,g,bcoef,coeffs,
     #                  maxwks,nsym,numrow,drtpt,z,iz,maxcor,
     #                  nnp0,ngsym,numsym,offsym,ijpt)
c
c***begin prologue    fmdm
c***date written      860911   (yymmdd)
c***revision date     yymmdd   (yymmdd)
c***keywords          ci, coupling coefficients, guga
c***author            saxe, paul,    (lanl)
c***purpose           to form the product of a vector with the hamiltonian
c                     matrix.
c***description
c
c          fmdm forms the product vector s of tha hamiltonian matrix
c      and a vector c. the hamiltonian matrix is formed on the 'fly'
c      from the one- and two-electron integrals in h and g, using the
c      drt information in uparc, dnarc, upwt, dnwt, upnwks, a, b,
c      levpt, and levnr. nlevs is the number of levels in the drt;
c      nrows, the number of rows.
c          the two electron coupling coefficients are formed from a
c      sum over all intermediate states k of the form:
c           e(ij,kl) = sum(k)   e(ik,ij)*e(kj,kl)
c
c***references
c
c***routines called
c***end prologue      fmdm
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c
c
c     ----- external unmodified arrays -----
c
      real*8 c(nwks),h(nnp0),g(ngsym),coeffs(420)
      integer ijpt(nnp),numsym(nsym),offsym(nsym)
      integer uparc(4,nrows),dnarc(4,nrows),upwt(4,nrows),dnwt(4,nrows)
      integer upnwks(nrows),dnnwks(nrows),a(nrows),b(nrows)
      integer levpt(nlevs),levnr(nlevs),rowoff(nrows)
      integer iuparc(4,inrows),idnarc(4,inrows),iupwt(4,inrows)
      integer idnwt(4,inrows)
      integer iupnwk(inrows),idnnwk(inrows),ia(inrows),ib(inrows)
      integer ilevpt(nlevs,nsym),ilevnr(nlevs,nsym)
      integer sval(nrows),isval(inrows)
      integer numrow(nsym),drtpt(nsym)
c
c     ----- external scratch arrays -----
c
      real*8 acoef(nlevs),bcoef(nlevs)
      real*8 z(*)
      integer iz(maxcor)
      integer irowsv(nlevs),jrowsv(nlevs),iwtsv(nlevs),jwtsv(nlevs)
      integer pagesv(nlevs),segsv(nlevs)
c
c     ----- common blocks -----
c
      common /io/ inp,iout
c
c     ----- zero the storage for the density matrices -----
c
      call rzero(h,nnp0)
      call rzero(g,ngsym)
c
c     ----- one-electron section -----
c
      jmin=levpt(levdiv)+1
      jmax=levpt(levdiv)+levnr(levdiv)
      nj=jmax-jmin+1
      botnum=1
      topnum=botnum+nj
      dnnum=topnum+nj
      upnum=dnnum+nj
      top=upnum+nj
      left=maxcor-top+1
c
      call onedm(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,a,b,
     #           levpt,levnr,rowoff,
     #           nrows,nlevs,norbs,
     #           acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,segsv,
     #           levdiv,nnp,c,nwks,h,bcoef,coeffs,
     #           maxwks,z(iadtwp(top)),iz(top),left,iz(botnum),
     #           iz(topnum),iz(dnnum),iz(upnum),jmin,jmax,nnp0,ijpt)
c
c
c     ----- two-electron section -----
c
      do 1 sym=1,nsym
         pt=drtpt(sym)
         jmin=levpt(levdiv)+1
         jmax=levpt(levdiv)+levnr(levdiv)
         nj=jmax-jmin+1
         botnum=1
         topnum=botnum+nj
         dnnum=topnum+nj
         upnum=dnnum+nj
         top=upnum+nj
         left=maxcor-top+1
c
         nnps=numsym(sym)
         ptsym=offsym(sym)+1
c
         call twodm(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,a,b,
     #              levpt,levnr,rowoff,
     #              iuparc(1,pt),idnarc(1,pt),iupwt(1,pt),
     #              idnwt(1,pt),iupnwk(pt),idnnwk(pt),
     #              ia(pt),ib(pt),ilevpt(1,sym),ilevnr(1,sym),
     #              numrow(sym),nrows,nlevs,norbs,
     #              acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,segsv,
     #              levdiv,nnps,c,nwks,h,g(ptsym),bcoef,coeffs,
     #              maxwks,sym-1,sval,isval(pt),z(iadtwp(top)),
     #              iz(top),left,iz(botnum),iz(topnum),iz(dnnum),
     #              iz(upnum),jmin,jmax,nnp,ijpt)
    1 continue
c
      return
      end
