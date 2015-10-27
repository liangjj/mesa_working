*deck @(#)twotdm.f	5.1  11/6/94
      subroutine twotdm(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,a,b,
     #                  levpt,levnr,rowoff,
     #                  iuparc,idnarc,iupwt,idnwt,iupnwk,idnnwk,
     #                  ia,ib,ilevpt,ilevnr,inrows,
     #                  nrows,nlevs,norbs,
     #                  acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,segsv,
     #                  levdiv,nnp,c,nwks,s,h,g,bcoef,coeffs,
     #                  maxwks,sym,sval,isval,z,iz,maxcor,botnum,
     #                  topnum,dnnum,upnum,jmin,jmax,nij,ijpt)
c
c***begin prologue     twotdm
c***date written       860916  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords
c***author             saxe, paul (lanl)
c***source
c***purpose            two-body portion of h.c using summation of states
c                      and one-body terms.
c***description
c
c***references
c***routines called
c***end prologue       twotdm
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c
c
c     ----- external unmodified arrays -----
c
      real*8 s(nwks)
      real*8 c(nwks),h(nnp),g(nnp,nnp),coeffs(420)
      integer ijpt(nij)
      integer uparc(4,nrows),dnarc(4,nrows),upwt(4,nrows),dnwt(4,nrows)
      integer upnwks(nrows),dnnwks(nrows),a(nrows),b(nrows)
      integer levpt(nlevs),levnr(nlevs),rowoff(nrows)
      integer iuparc(4,inrows),idnarc(4,inrows),iupwt(4,inrows)
      integer idnwt(4,inrows)
      integer iupnwk(inrows),idnnwk(inrows),ia(inrows),ib(inrows)
      integer ilevpt(nlevs),ilevnr(nlevs)
      integer sval(nrows),isval(inrows)
c
c     ----- external scratch arrays -----
c
      real*8 z(*)
      real*8 acoef(nlevs),bcoef(nlevs)
      integer iz(maxcor)
      integer irowsv(nlevs),jrowsv(nlevs),iwtsv(nlevs),jwtsv(nlevs)
      integer pagesv(nlevs),segsv(nlevs)
      integer botnum(jmin:jmax),topnum(jmin:jmax)
      integer dnnum(jmin:jmax),upnum(jmin:jmax)
c
c     ----- local arrays -----
c
      integer nc(10)
c
c     ----- common blocks -----
c
      common /io/ inp,iout
c
      maxwp=iadtwp(maxcor)
      call izero(nc,10)
      minrow=ilevpt(levdiv)+1
      maxrow=ilevpt(levdiv)+ilevnr(levdiv)
      do 4000 irow=minrow,maxrow
c
c        ----- determine how much we can do at a time -----
c
         kdnnwk=idnnwk(irow)
         kupnwk=iupnwk(irow)
c
c        ----- form the partial loops and store results -----
c
         pt=1
         call izero(botnum,jmax-jmin+1)
         call izero(topnum,jmax-jmin+1)
         call izero(dnnum,jmax-jmin+1)
         call izero(upnum,jmax-jmin+1)
         maxov=0
         do 200 jrow=jmin,jmax
            deltaa=a(jrow)-ia(irow)
            deltab=b(jrow)-ib(irow)
c
            if (deltaa.eq.0.and.(deltab.eq.-1.or.deltab.eq.1).or.
     #          deltaa.eq.-1.and.deltab.eq.1.or.
     #          deltaa.eq.1.and.deltab.eq.-1) then
c
c
                  call top(iuparc,iupwt,ia,ib,iupnwk,inrows,
     #                     uparc,upwt,a,b,upnwks,nrows,
     #                     acoef,irowsv,jrowsv,iwtsv,jwtsv,
     #                     pagesv,segsv,
     #                     irow,jrow,levdiv,nlevs,coeffs,
     #                     z(pt),z(pt),wptoin(4),(maxwp-pt)/4,ntops)
                  topnum(jrow)=ntops
                  pt=pt+ntops*4
c
                  call bottom(idnarc,idnwt,ia,ib,idnnwk,inrows,
     #                        dnarc,dnwt,a,b,dnnwks,nrows,
     #                        acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,
     #                        segsv,irow,jrow,levdiv,nlevs,coeffs,
     #                        z(pt),z(pt),wptoin(4),(maxwp-pt)/4,nbots)
c
               maxov=maxov+nbots
                  botnum(jrow)=nbots
                  pt=pt+nbots*4
c
            end if
  200    continue
c
c        ----- loops entirely in the top or bottom half of drt -----
c
         do 300 jrow=jmin,jmax
            deltaa=a(jrow)-ia(irow)
            deltab=b(jrow)-ib(irow)
            deltau=xor(xor(sval(jrow),sval(1)),xor(isval(irow),sym))
            deltal=xor(sval(jrow),isval(irow))
c
            if (deltaa.ne.0.or.deltab.ne.0) go to 300
c
            jupnwk=upnwks(jrow)
            jdnnwk=dnnwks(jrow)
c
            if (deltau.eq.0) then
               call dnloop(idnarc,idnwt,ia,ib,idnnwk,inrows,
     #                     dnarc,dnwt,a,b,dnnwks,nrows,
     #                     acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,
     #                     segsv,irow,jrow,levdiv,nlevs,coeffs,
     #                     z(pt),z(pt),wptoin(4),(maxwp-pt)/4,ndnlps,
     #                     nij,ijpt)
c
               dnnum(jrow)=ndnlps
               pt=pt+ndnlps*4
c
               call upwalk(iuparc,iupwt,ia,ib,inrows,
     #                     uparc,a,b,nrows,
     #                     irowsv,jrowsv,iwtsv,segsv,
     #                     irow,jrow,levdiv,nlevs,z(pt),jupnwk)
c
               pt=pt+jupnwk
            end if
            if (deltal.eq.0) then
c
               call uploop(iuparc,iupwt,ia,ib,iupnwk,inrows,
     #                     uparc,upwt,a,b,upnwks,nrows,
     #                     acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,
     #                     segsv,irow,jrow,levdiv,nlevs,coeffs,
     #                     z(pt),z(pt),wptoin(4),(maxwp-pt)/4,nuplps,
     #                     nij,ijpt)
               upnum(jrow)=nuplps
               pt=pt+nuplps*4
c
               call dnwalk(idnarc,idnwt,ia,ib,inrows,
     #                     dnarc,a,b,nrows,
     #                     irowsv,jrowsv,iwtsv,segsv,
     #                     irow,jrow,levdiv,nlevs,z(pt),jdnnwk)
               pt=pt+jdnnwk
            end if
  300    continue
c
c        ----- allocate scratch space -----
c
         asort=pt
         jsort=asort+kdnnwk
         jnsort=jsort+kdnnwk
         aov=jnsort+kdnnwk
         jov=aov+maxov
         kov=jov+maxov
         jnov=kov+maxov
         ec=jnov+maxov
c
c        ----- determine how much we can do at a time -----
c
         left=maxwp-ec
         ninter=left/nnp/2
         nup=min(kupnwk,ninter/kdnnwk)
         if (nup.lt.1) then
         write(iout,*)' maxwp ec ',maxwp,ec
         write(iout,*)' kupnwk ninter kdnnwk ',kupnwk,ninter,
     1   kdnnwk
          call lnkerr('not enough core')
         end if
         nkwks=nup*kdnnwk
c
         es=ec+nnp*nkwks
         topcor=es+nnp*nkwks
c
         if (topcor.gt.maxwp) then
            call lnkerr('core allocation for ec & es')
         end if
c
         maxk=0
 2100    continue
            mink=maxk+1
            maxk=min(kupnwk,maxk+nup)
            nupks=maxk-mink+1
            nkwks=nupks*kdnnwk
c
            call rzero(z(ec),nnp*nkwks)
            call rzero(z(es),nnp*nkwks)
c
            pt=1
c
            do 3000 jrow=jmin,jmax
               deltaa=a(jrow)-ia(irow)
               deltab=b(jrow)-ib(irow)
c
               if (deltaa.eq.0.and.(deltab.eq.-1.or.deltab.eq.1).or.
     #             deltaa.eq.-1.and.deltab.eq.1.or.
     #             deltaa.eq.1.and.deltab.eq.-1) then
c
                  ntops=topnum(jrow)
                  pt1=pt+ntops*4
                  nbots=botnum(jrow)
                  if (ntops.eq.0) go to 2900
                  if (nbots.eq.0) go to 2900
c
                  call ectb(z(ec),nnp,nkwks,z(pt),z(pt),wptoin(4),
     #                      ntops,z(pt1),z(pt1),nbots,dnnwks(jrow),
     #                      idnnwk(irow),rowoff(jrow),c,nwks,
     #                      z(asort),z(jsort),z(jnsort),z(aov),
     #                      z(jov),z(kov),z(jnov),maxov,mink,maxk,
     #                      nij,ijpt)
                  call ectb(z(es),nnp,nkwks,z(pt),z(pt),wptoin(4),
     #                      ntops,z(pt1),z(pt1),nbots,dnnwks(jrow),
     #                      idnnwk(irow),rowoff(jrow),s,nwks,
     #                      z(asort),z(jsort),z(jnsort),z(aov),
     #                      z(jov),z(kov),z(jnov),maxov,mink,maxk,
     #                      nij,ijpt)
                  nc(1)=nc(1)+nbots
                  nc(2)=nc(2)+ntops
 2900             continue
                  pt=pt+4*(ntops+nbots)
               end if
 3000       continue
c
c
            do 3100 jrow=jmin,jmax
               deltaa=a(jrow)-ia(irow)
               deltab=b(jrow)-ib(irow)
               deltau=xor(xor(sval(jrow),sval(1)),xor(isval(irow),sym))
               deltal=xor(sval(jrow),isval(irow))
c
               if (deltaa.ne.0.or.deltab.ne.0) go to 3100
c
               jupnwk=upnwks(jrow)
               jdnnwk=dnnwks(jrow)
c
               if (deltau.eq.0) then
                  ndnlps=dnnum(jrow)
                  pt1=pt+4*ndnlps
                  nc(3)=nc(3)+ndnlps
                  nc(4)=nc(4)+jupnwk
c
                  call ecdn(z(ec),nnp,nkwks,z(pt),z(pt),wptoin(4),
     #                      ndnlps,kdnnwk,jdnnwk,jupnwk,rowoff(jrow),
     #                      c,nwks,z(pt1),mink,maxk)
                  call ecdn(z(es),nnp,nkwks,z(pt),z(pt),wptoin(4),
     #                      ndnlps,kdnnwk,jdnnwk,jupnwk,rowoff(jrow),
     #                      s,nwks,z(pt1),mink,maxk)
                  pt=pt1+jupnwk
               end if
c
               if (deltal.eq.0) then
c
                  nuplps=upnum(jrow)
                  pt1=pt+nuplps*4
                  nc(5)=nc(5)+nuplps
                  nc(6)=nc(6)+jdnnwk
c
                  call ecup(z(ec),nnp,nkwks,z(pt),z(pt),wptoin(4),
     #                      nuplps,kdnnwk,jdnnwk,rowoff(jrow),
     #                      c,nwks,z(pt1),mink,maxk)
                  call ecup(z(es),nnp,nkwks,z(pt),z(pt),wptoin(4),
     #                      nuplps,kdnnwk,jdnnwk,rowoff(jrow),
     #                   s,nwks,z(pt1),mink,maxk)
c
                  pt=pt1+jdnnwk
               end if
 3100       continue
c
ci            call ebc(z(gec),g,z(ec),nnp,nnp,nkwks)
            call apbct(g,z(es),z(ec),nnp,nkwks,nnp)
c
         if (maxk.lt.kupnwk) go to 2100
c
 4000 continue
cps      write (iout,4001) sym,nc
 4001 format (1x,i2,10i6)
c
c
c
c
      return
      end
