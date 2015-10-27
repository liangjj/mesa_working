*deck @(#)gdrt.f	5.1  11/6/94
      subroutine gdrt(nrows,nlevs,levpt,levnr,dnarc,uparc,dnnwks,
     #                  upnwks,dnwt,upwt,aval,bval,sval,
     #                  inrows,numrow,drtpt,nsym,ilevpt,ilevnr,idnarc,
     #                  iuparc,idnnwk,iupnwk,idnwt,iupwt,iaval,
     #                  ibval,isval)
c
c***begin prologue     gdrt
c***date written       860916  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords
c***author             saxe, paul (lanl)
c***source
c***purpose            to retrieve the main and intermediate drt's.
c***description
c
c***references
c***routines called
c***end prologue       gdrt
c
      implicit integer (a-z)
c
      integer levpt(nlevs),levnr(nlevs),dnarc(4,nrows),uparc(4,nrows)
      integer dnnwks(nrows),upnwks(nrows),dnwt(4,nrows),upwt(4,nrows)
      integer aval(nrows),bval(nrows),sval(nrows)
      integer ilevpt(nlevs,nsym),ilevnr(nlevs,nsym)
      integer idnarc(4,inrows),iuparc(4,inrows)
      integer idnnwk(inrows),iupnwk(inrows),idnwt(4,inrows)
      integer iupwt(4,inrows),iaval(inrows),ibval(inrows),isval(inrows)
      integer numrow(nsym),drtpt(nsym)
      character*2 itoc,symchr
c
c     ----- read in the main drt information -----
c
      call iosys('read integer levpt from rwf',nlevs,levpt,0,' ')
      call iosys('read integer levnr from rwf',nlevs,levnr,0,' ')
      call iosys('read integer arc   from rwf',4*nrows,dnarc,0,' ')
      call iosys('read integer nlwks from rwf',nrows,dnnwks,0,' ')
      call iosys('read integer weight from rwf',4*nrows,dnwt,0,' ')
      call iosys('read integer a from rwf',nrows,aval,0,' ')
      call iosys('read integer b from rwf',nrows,bval,0,' ')
      call iosys('read integer s from rwf',nrows,sval,0,' ')
c
c     ----- create the remaining arrays -----
c
      call updrt(dnarc,uparc,dnnwks,upnwks,dnwt,upwt,levpt,levnr,
     #           nlevs,nrows)
c
c     ----- and the intermediate drt's for all symmetries -----
c
      do 1 sym=1,nsym
         symchr=itoc(sym-1)
         pt=drtpt(sym)
         n=numrow(sym)
c
         if (n.eq.0) go to 1
c
         call iosys('read integer "'//symchr//'levpt" from rwf',nlevs,
     #                         ilevpt(1,sym),0,' ')
         call iosys('read integer "'//symchr//'levnr" from rwf',nlevs,
     #                         ilevnr(1,sym),0,' ')
         call iosys('read integer "'//symchr//'arc"   from rwf',4*n,
     #                        idnarc(1,pt),0,' ')
         call iosys('read integer "'//symchr//'nlwks" from rwf',n,
     #                        idnnwk(pt),0,' ')
         call iosys('read integer "'//symchr//'weight" from rwf',4*n,
     #                        idnwt(1,pt),0,' ')
         call iosys('read integer "'//symchr//'a" from rwf',n,
     #                         iaval(pt),0,' ')
         call iosys('read integer "'//symchr//'b" from rwf',n,
     #                         ibval(pt),0,' ')
         call iosys('read integer "'//symchr//'s" from rwf',n,
     #                         isval(pt),0,' ')
c
c        ----- create the remaining arrays -----
c
         call updrt(idnarc(1,pt),iuparc(1,pt),idnnwk(pt),iupnwk(pt),
     #              idnwt(1,pt),iupwt(1,pt),ilevpt(1,sym),ilevnr(1,sym),
     #              nlevs,n)
    1 continue
c
c
      return
      end
