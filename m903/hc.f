*deck @(#)hc.f	5.1  11/6/94
      subroutine hc(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,a,b,
     #              levpt,levnr,rowoff,
     #              iuparc,idnarc,iupwt,idnwt,iupnwk,idnnwk,ia,ib,
     #              ilevpt,ilevnr,inrows,
     #              nrows,nlevs,norbs,
     #              acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,segsv,
     #              levdiv,nnp,c,nwks,s,h,g,bcoef,
     #              prtflg,nsym,numrow,drtpt,
     #              z,sval,isval,maxwks,ijpt,nnp0,ngsym,numsym,
     #              offsym,orbsym,h1,g1,ops)
c
      implicit integer (a-z)
c
c     ----- external unmodified arrays -----
c
      real*8 h1(nnp),g1(nnp,nnp)
      real*8 h(nnp0),g(ngsym)
      integer numsym(nsym),offsym(nsym),ijpt(nnp),orbsym(norbs)
      integer uparc(4,nrows),dnarc(4,nrows),upwt(4,nrows),dnwt(4,nrows)
      integer upnwks(nrows),dnnwks(nrows),a(nrows),b(nrows)
      integer levpt(nlevs),levnr(nlevs),rowoff(nrows)
      integer iuparc(4,inrows),idnarc(4,inrows),iupwt(4,inrows)
      integer idnwt(4,inrows)
      integer iupnwk(inrows),idnnwk(inrows),ia(inrows),ib(inrows)
      integer ilevpt(nlevs,nsym),ilevnr(nlevs,nsym)
      integer numrow(nsym),drtpt(nsym)
      integer sval(nrows),isval(inrows)
c
c     ----- external scratch arrays -----
c
      real*8 c(nwks),s(nwks)
      real*8 acoef(nlevs),bcoef(nlevs)
      real*8 z(*)
      integer irowsv(nlevs),jrowsv(nlevs),iwtsv(nlevs),jwtsv(nlevs)
      integer pagesv(nlevs),segsv(nlevs)
c
c     ----- external unmodified scalars -----
c
      character*(*) ops
      character*8 prtflg
      integer nsym
c
c     ----- local arrays -----
c
      real*8 coeffs(20,21)
c
c     ----- local scalars -----
c
      real*8 t
      logical debug
c
c     ----- external functions -----
c
c
c     ----- local data statements -----
c
      parameter (debug=.false.)
c
c     ----- common blocks -----
c
      common /io/ inp,iout
c
c     ----- build the coefficient table -----
c
      call rzero(coeffs,20*21)
c
      do 700 i=3,20
         t = float(i-2)
         coeffs(i,1) = sqrt(t/(t+1.0d+00))
         coeffs(i,2) = -coeffs(i,1)
         coeffs(i,3) = coeffs(i,1)/sqrt(2.0d+00)
         coeffs(i,4) = -coeffs(i,3)
         coeffs(i,5) = sqrt((t+1.0d+00)/t)
         coeffs(i,6) = -coeffs(i,5)
         coeffs(i,7) = coeffs(i,5)/sqrt(2.0d+00)
         coeffs(i,8) = -coeffs(i,7)
         coeffs(i,9) = sqrt((t+2.0d+00)/(t*2.0d+00))
         coeffs(i,10) = -coeffs(i,9)
         coeffs(i,11) = sqrt(t/(2.0d+00*(t+2.0d+00)))
         coeffs(i,12) = -coeffs(i,11)
         coeffs(i,13) = sqrt(2.0d+00/(t*(t+1.0d+00)))
         coeffs(i,14) = sqrt(t*(t+2.0d+00))/(t+1.0d+00)
         coeffs(i,15) = -sqrt(t*(t+2.0d+00))/(t+1.0d+00)
         coeffs(i,16) = sqrt((t-1.0d+00)*(t+2.0d+00)/(t*(t+1.0d+00)))
         coeffs(i,17) = -coeffs(i,16)
         coeffs(i,18)=-sqrt(2.0d+00/(t*(t+2.0d+00)))
         coeffs(i,19) = 1.0d+00/t
         coeffs(i,20) = -1.0d+00/t
         coeffs(i,21) = -sqrt(2.0d+00)/t
  700 continue
c
c     ----- get the main and intermediate drt's -----
c
      call gdrt(nrows,nlevs,levpt,levnr,dnarc,uparc,dnnwks,
     #            upnwks,dnwt,upwt,a,b,sval,
     #            inrows,numrow,drtpt,nsym,ilevpt,ilevnr,idnarc,
     #            iuparc,idnnwk,iupnwk,idnwt,iupwt,ia,ib,isval)
c
c     ----- create the offset pointers for the rows on the dividing level
c
      n=0
      do 1 row=levpt(levdiv)+1,levpt(levdiv)+levnr(levdiv)
         rowoff(row)=n
         n=n+upnwks(row)*dnnwks(row)
    1 continue
      if (n.ne.dnnwks(1)) then
         call lnkerr('did not get nwks while creating offsets')
      end if
c
      n=0
      do 3 sym=1,nsym
         do 2 row=ilevpt(levdiv,sym)+1,
     #                ilevpt(levdiv,sym)+ilevnr(levdiv,sym)
            pt=drtpt(sym)-1+row
            n=max(n,iupnwk(pt)*idnnwk(pt))
    2   continue
    3 continue
c
c
c     ----- finish core allocation -----
c
c
      if (maxwks.le.0) then
         maxwks=n
      else
         maxwks=min(n,maxwks)
      end if
c
      if (prtflg.ne.'minimum') then
         write (iout,920) maxwks,n
 920     format (/,' m903:',/,5x,'attempting to use ',i7,
     #        ' intermediates at a time',/,
     #        5x,'of a possible     ',i7)
      end if
c
      call getscm(0,z,canget,'m903 available',0)
      top=min(canget,wpadti(2*nnp*maxwks+10000))
      call getscm(top,z,maxcor,'m903 ec and gec',0)
c
      if(debug) then
          write (iout,921) (rowoff(i),i=levpt(levdiv)+1,levpt(levdiv)+
     #                                                  levnr(levdiv))
  921     format (/,' row offsets:',/,(1x,10i8))
      end if
c
c     ----- form combinations of integrals -----
c
      call fmints(h1,g1,norbs,nnp)
c
c     ----- symmetry block the integrals -----
c
      call rfmint(orbsym,norbs,nsym,numsym,offsym,ijpt,nnp0,
     #            nnp,ngsym,h,h1,g,g1)
c
c     ----- multiply the c vector by the hamiltonian -----
c
      call rzero(s,nwks)
c
      call multhc(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,
     #     a,b,sval,levpt,levnr,rowoff,
     #     iuparc,idnarc,iupwt,idnwt,iupnwk,idnnwk,
     #     ia,ib,isval,ilevpt,ilevnr,inrows,
     #     nrows,nlevs,norbs,
     #     acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,
     #     segsv,levdiv,nnp,c,nwks,s,h,
     #     g,bcoef,coeffs,maxwks,
     #     nsym,numrow,drtpt,z,z,maxcor,
     #     nnp0,ngsym,numsym,offsym,ijpt)
c
      if(debug) then
         write (iout,443) (c(i),i=1,nwks)
  443    format (/,5x,' trial:',/,(10x,10f10.6))
         write (iout,444) (s(i),i=1,nwks)
  444    format (5x,' product:',/,(10x,10f10.6))
      endif
c
c
      return
      end
