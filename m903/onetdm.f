*deck @(#)onetdm.f	5.1  11/6/94
      subroutine onetdm(uparc,dnarc,upwt,dnwt,upnwks,dnnwks,a,b,
     #                  levpt,levnr,rowoff,
     #                  nrows,nlevs,norbs,
     #                  acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,segsv,
     #                  levdiv,nnp,c,nwks,s,h,bcoef,coeffs,
     #                  maxwks,z,iz,maxcor,botnum,topnum,dnnum,upnum,
     #                  jmin,jmax,nnp0,ijpt)
c
c***begin prologue     onetdm
c***date written       860916  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords
c***author             saxe, paul (lanl)
c***source
c***purpose            one-body portion of h.c using summation of states
c                      and one-body terms.
c***description
c
c***references
c***routines called
c***end prologue       onetdm
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c
c
c     ----- external unmodified arrays -----
c
      real*8 s(nwks)
      real*8 c(nwks),h(nnp0),coeffs(420)
      integer ijpt(nnp)
      integer uparc(4,nrows),dnarc(4,nrows),upwt(4,nrows),dnwt(4,nrows)
      integer upnwks(nrows),dnnwks(nrows),a(nrows),b(nrows)
      integer levpt(nlevs),levnr(nlevs),rowoff(nrows)
c
c     ----- external scratch arrays -----
c
      real*8 acoef(nlevs),bcoef(nlevs),z(*)
      integer iz(maxcor)
      integer irowsv(nlevs),jrowsv(nlevs),iwtsv(nlevs),jwtsv(nlevs)
      integer pagesv(nlevs),segsv(nlevs)
      integer topnum(jmin:jmax),botnum(jmin:jmax)
      integer upnum(jmin:jmax),dnnum(jmin:jmax)
c
c     ----- common blocks -----
c
      common /io/ inp,iout
c
c     ----- loop through rows, forming one-body matrix elements -----
c
      maxwp=iadtwp(maxcor)
      do 2000 irow=jmin,jmax
         kdnnwk=dnnwks(irow)
         kupnwk=upnwks(irow)
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
            deltaa=a(jrow)-a(irow)
            deltab=b(jrow)-b(irow)
c
            if (deltaa.eq.0.and.(deltab.eq.-1.or.deltab.eq.1).or.
     #          deltaa.eq.-1.and.deltab.eq.1.or.
     #          deltaa.eq.1.and.deltab.eq.-1) then
c
c
               call top(uparc,upwt,a,b,upnwks,nrows,
     #                  uparc,upwt,a,b,upnwks,nrows,
     #                  acoef,irowsv,jrowsv,iwtsv,jwtsv,
     #                  pagesv,segsv,
     #                  irow,jrow,levdiv,nlevs,coeffs,
     #                  z(pt),z(pt),wptoin(4),(maxwp-pt)/4,ntops)
               topnum(jrow)=ntops
               pt=pt+ntops*4
c
               call bottom(dnarc,dnwt,a,b,dnnwks,nrows,
     #                     dnarc,dnwt,a,b,dnnwks,nrows,
     #                     acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,
     #                     segsv,irow,jrow,levdiv,nlevs,coeffs,
     #                     z(pt),z(pt),wptoin(4),(maxwp-pt)/4,nbots)
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
         call dnloop(dnarc,dnwt,a,b,dnnwks,nrows,
     #               dnarc,dnwt,a,b,dnnwks,nrows,
     #               acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,
     #               segsv,irow,irow,levdiv,nlevs,coeffs,
     #               z(pt),z(pt),wptoin(4),(maxwp-pt)/4,ndnlps,nnp,ijpt)
         dnnum(irow)=ndnlps
         pt=pt+ndnlps*4
c
         call upwalk(uparc,upwt,a,b,nrows,
     #               uparc,a,b,nrows,
     #               irowsv,jrowsv,iwtsv,segsv,
     #               irow,irow,levdiv,nlevs,z(pt),kupnwk)
         pt=pt+kupnwk
c
         call uploop(uparc,upwt,a,b,upnwks,nrows,
     #               uparc,upwt,a,b,upnwks,nrows,
     #               acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,
     #               segsv,irow,irow,levdiv,nlevs,coeffs,
     #               z(pt),z(pt),wptoin(4),(maxwp-pt)/4,nuplps,nnp,ijpt)
         upnum(irow)=nuplps
         pt=pt+nuplps*4
c
         call dnwalk(dnarc,dnwt,a,b,nrows,
     #               dnarc,a,b,nrows,
     #               irowsv,jrowsv,iwtsv,segsv,
     #               irow,irow,levdiv,nlevs,z(pt),kdnnwk)
         pt=pt+kdnnwk
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
         ninter=left/nnp0
         nup=min(kupnwk,ninter/kdnnwk)
         if (nup.lt.1) then
            write(iout,*)' maxwp ec nnp0 ',maxwp,ec,nnp0
            write(iout,*)' kupnwk ninter kdnnwk ',kupnwk,ninter,
     #                     kdnnwk
            call lnkerr('not enough core step 2 ')
         end if
c
         maxk=0
  100    continue
            mink=maxk+1
            maxk=min(kupnwk,maxk+nup)
            nupks=maxk-mink+1
            nkwks=nupks*kdnnwk
c
            if (ec+nnp0*nkwks-1.gt.maxwp) then
               call lnkerr('core allocation error')
            end if
c
            call rzero(z(ec),nnp0*nkwks)
c
c
            pt=1
            do 1000 jrow=jmin,jmax
               deltaa=a(jrow)-a(irow)
               deltab=b(jrow)-b(irow)
c
               if (deltaa.eq.0.and.(deltab.eq.-1.or.deltab.eq.1).or.
     #             deltaa.eq.-1.and.deltab.eq.1.or.
     #             deltaa.eq.1.and.deltab.eq.-1) then
c
                  ntops=topnum(jrow)
                  pt1=pt+topnum(jrow)*4
                  nbots=botnum(jrow)
c
                  call ectb(z(ec),nnp0,nkwks,z(pt),z(pt),wptoin(4),
     #                      ntops,z(pt1),z(pt1),nbots,dnnwks(jrow),
     #                      dnnwks(irow),rowoff(jrow),c,nwks,
     #                      z(asort),z(jsort),z(jnsort),z(aov),
     #                      z(jov),z(kov),z(jnov),maxov,mink,maxk,
     #                      nnp,ijpt)
c
                  pt=pt+4*(ntops+nbots)
c
               end if
 1000       continue
c
c
            ndnlps=dnnum(irow)
            pt1=pt+4*ndnlps
c
            call ecdn(z(ec),nnp0,nkwks,z(pt),z(pt),wptoin(4),ndnlps,
     #                kdnnwk,kdnnwk,kupnwk,rowoff(irow),
     #                c,nwks,z(pt1),mink,maxk)
c
            pt=pt1+kupnwk
            nuplps=upnum(irow)
            pt1=pt+nuplps*4
c
            call ecup(z(ec),nnp0,nkwks,z(pt),z(pt),wptoin(4),nuplps,
     #                kdnnwk,kdnnwk,rowoff(irow),
     #                c,nwks,z(pt1),mink,maxk)
c
ci            call apbc(s(rowoff(irow)+1+(mink-1)*kdnnwk),
ci     #                h,z(ec),1,nnp0,nkwks)
            call apbc(h,z(ec),s(rowoff(irow)+1+(mink-1)*kdnnwk),
     $           nnp0,nkwks,1)
c
c
         if (maxk.lt.kupnwk) go to 100
c
 2000 continue
c
      return
      end
