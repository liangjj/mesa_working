*deck @(#)ci.f	5.1  11/6/94
      subroutine ci(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,a,b,
     #              levpt,levnr,rowoff,
     #              iuparc,idnarc,iupwt,idnwt,iupnwk,idnnwk,ia,ib,
     #              ilevpt,ilevnr,inrows,
     #              nrows,nlevs,norbs,
     #              acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,segsv,
     #              levdiv,nnp,c,nwks,s,h,g,bcoef,
     #              root,dvec,dmat,nroots,mxiter,cnverg,thresh,
     #              nguess,prtflg,refwlk,nsym,numrow,drtpt,
     #              z,sval,isval,maxwks,ijpt,nnp0,ngsym,numsym,
     #              offsym,orbsym,h1,g1,ops,nattim,mxvc,froot,vn)
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
      real*8 c(nwks,mxvc),s(nwks,mxvc),root(nroots),dvec(mxiter,mxiter)
      real*8 dmat((mxiter+1)*mxiter/2)
      real*8 acoef(nlevs),bcoef(nlevs)
      real*8 z(*)
      integer irowsv(nlevs),jrowsv(nlevs),iwtsv(nlevs),jwtsv(nlevs)
      integer pagesv(nlevs),segsv(nlevs),refwlk(nroots)
c
c     ----- external unmodified scalars -----
c
      character*(*) ops
      character*2 vn
      real*8 cnverg,thresh
      character*8 prtflg
      integer nroots,mxiter,nguess,nsym
c
c     ----- local arrays -----
c
      real*8 coeffs(20,21)
c
c     ----- local scalars -----
c
      real*8 t,rep,fzcore,dmin,czero,eguess,eci
      character*16 status
      character*3 answer
      logical debug
c
c     ----- external functions -----
c
      character*4 itoc
      logical logkey
c
c     ----- common blocks -----
c
      common /io/ inp,iout
c
c     ----- local data statements -----
c
      parameter (debug=.false.)
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
c     ----- get the nuclear repulsion and frozen core energies -----
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,rep,0,' ')
      call iosys('does "frozen core energy" exist on rwf',0,0,0,answer)
      if (answer.eq.'no') then
         fzcore=0.0d+00
      else
         call iosys('read real "frozen core energy" from rwf',1,fzcore,
     #               0,' ')
      end if
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
         write (iout,921) (rowoff(i),i=levpt(levdiv)+1,
     #                                 levpt(levdiv)+levnr(levdiv))
  921    format (/,' row offsets:',/,(1x,10i8))
      end if
c
c     ----- form the diagonals in lexical order -----
c
      call rzero(c,nwks)
      call dndiag(dnarc,dnwt,a,b,dnnwks,nrows,nlevs,coeffs,acoef,bcoef,
     #            irowsv,iwtsv,pagesv,segsv,c,nwks,h1,g1,nnp,jrowsv)
c
c     ----- convert from lexical to divided order -----
c
      call reord(dnarc,upwt,dnwt,dnnwks,rowoff,nrows,nlevs,
     #            irowsv,iwtsv,jwtsv,segsv,levdiv,c,s,nwks)
c
      if (logkey(ops,'print=ci=diagonals',.false.,' ')) then
         write (iout,230) (s(i,1),i=1,nwks)
 230     format (' diagonals:',/,(1x,7f10.6))
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
c     ----- save diagonals on rwf, and then pass to bliu -----
c
      call iosys('write real "h diagonals" to rwf',nwks,s,0,' ')
c
      if (prtflg.eq.'minimum') then
         status='noprint'
      else
         status=' '
      end if
c
      call bliu('initialize with diagonals',status,s,thresh,nwks,
     #           mxiter,nroots,iout,nattim,rep+fzcore,cnverg,dvec)
c
c
c     ----- find the lowest diagonals for the guess -----
c
      do 20 guess=1,nguess
         dmin=1.0d+03
         do 15 i=1,nwks
            if (s(i,1).lt.dmin) then
               do 10 j=1,guess-1
                  if (i.eq.refwlk(j)) go to 14
   10          continue
               refwlk(guess)=i
               dmin=s(i,1)
   14          continue
            end if
   15    continue
   20 continue
c
c     ----- multiply the guess vectors by the hamiltonian -----
c
      do 30 guess=1,nguess
         call rzero(c,nwks)
         call rzero(s,nwks)
         c(refwlk(guess),1)=1.0d+00
c
         call multhc(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,
     #               a,b,sval,levpt,levnr,rowoff,
     #               iuparc,idnarc,iupwt,idnwt,iupnwk,idnnwk,
     #               ia,ib,isval,ilevpt,ilevnr,inrows,
     #               nrows,nlevs,norbs,
     #               acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,
     #               segsv,levdiv,nnp,c,nwks,s,h,
     #               g,bcoef,coeffs,maxwks,
     #               nsym,numrow,drtpt,z,z,maxcor,
     #               nnp0,ngsym,numsym,offsym,ijpt)
c
         if (logkey(ops,'print=ci=trials',.false.,' ')) then
            write (iout,443) guess,(c(i,1),i=1,nwks)
 443        format (/,5x,' trial:',i3,/,(10x,10f10.6))
            write (iout,444) (s(i,1),i=1,nwks)
 444        format (5x,' product:',/,(10x,10f10.6))
         end if
c
         call bliu('with vectors',0,c,s,nwks,mxiter,0,0,0,
     $              dmat,root,dvec)
c
   30 continue
c
c     ----- now commence the iterations -----
c
c.rlm     mxvc=1
   40 continue
         call bliu('solve',status,c,s,nwks,mxiter,0,0,0,dmat,root,dvec)
         if (status.eq.'converged') go to 100
c
         nvc=0
   50    continue
            nvc=nvc+1
            call bliu('new trial',status,c(1,nvc),s(1,nvc),nwks,mxiter,
     #                0,0,0,0,root,dvec)
            if (status.ne.'done'.and.nvc.lt.mxvc) go to 50
            if (status.eq.'done'.and.nvc.eq.1) go to 40
            if (status.eq.'done') nvc=nvc-1
c
            call rzero(s,nwks*nvc)
c
            do 60 ivc=1,nvc
               call multhc(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,
     #                  a,b,sval,levpt,levnr,rowoff,
     #                  iuparc,idnarc,iupwt,idnwt,iupnwk,idnnwk,
     #                  ia,ib,isval,ilevpt,ilevnr,inrows,
     #                  nrows,nlevs,norbs,
     #                  acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,
     #                  segsv,levdiv,nnp,c(1,ivc),nwks,s(1,ivc),h,
     #                  g,bcoef,coeffs,maxwks,
     #                  nsym,numrow,drtpt,z,z,maxcor,
     #                  nnp0,ngsym,numsym,offsym,ijpt)
c
               if(debug) then
                  write(iout,83) ivc,(c(iq,ivc),iq=1,nwks)
   83             format(1x,'vector ivc=',i5,/,(1x,10f12.6))
                  write(iout,84) ivc,(s(iq,ivc),iq=1,nwks)
   84             format(1x,'sigma  ivc=',i5,/,(1x,10f12.6))
               end if
               call bliu('with vectors',0,c(1,ivc),s(1,ivc),
     #                    nwks,mxiter,0,0,0,0,0,0)
   60       continue
c
         if (status.eq.'done') go to 40
         nvc=0
      go to 50
c
c
  100 continue
c
c     ----- recover each root's ci vector, find the most important
c           configuration, and write vector to rwf
c
      call iosys('read real diagonals from bliu',-1,s,0,' ')
c
      if (prtflg.ne.'minimum') then
         write (iout,90)
   90    format (///,' root   reference  guess energy    ci energy  ',
     #               '   c(0)')
      end if
c
      do 200 iroot=1,nroots
         call bliu('get vector',status,c,0,nwks,0,iroot,0,0,0,0,0)
         if (status.ne.'ok') go to 200
         call iosys('write real "'//vn//' root '//itoc(iroot)
     #               //'" to rwf',nwks,c,0,' ')
c
         czero=0.0d+00
         ref=0
         do 110 i=1,nwks
            if (abs(c(i,1)).gt.abs(czero)) then
               czero=c(i,1)
               ref=i
            end if
  110    continue
c
         eguess=s(ref,1)+rep+fzcore
         eci=root(iroot)+rep+fzcore
cps      edav=eci+(eci-eguess)*(1.0d+00-czero**2)
c
         if (iroot.eq.froot) then
            call iosys('write real energy to rwf',1,eci,0,' ')
         end if
c
         if (prtflg.ne.'minimum') then
            write (iout,120) iroot,ref,eguess,eci,czero
  120       format (1x,i3,i10,2g18.9,f8.4)
         end if
c
  200 continue
c
      call bliu('finish',0,0,0,0,0, 0,0,0,0,0,0)
c
      return
      end
