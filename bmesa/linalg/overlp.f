      subroutine overlp (smat,ipvt,soln,flago,wts,ncst,nlag,jind,nsts
     1 ,nbk,nbkpls,brk,c,sc,spln,nflag,ncmax,ntchn,ldim,msize,nptmx,nsol
     2 ,nds1,nds2,nds3,nds4,ien,ops)
      implicit integer(a-z)
      logical logkey
      common /io/ inp, iout
      character *(*) ops
      character *80 wrtspn
      character *3 itoc
      real *8smat, soln, flago, wts, brk, c, sc, spln
      dimension smat(nflag,nsol), flago(nptmx,ldim,nflag)
      dimension soln(msize,nsol), wts(msize), ipvt(nflag), ncst(nsts)
      dimension jind(ncmax,nsts), nlag(nsts)
      dimension brk(nds1), c(nds2), sc(nds3), spln(nds4)
      wrtspn(n1)='write real "spline vecs'//itoc(n1)//'" to wavefn witho
     1ut rewinding'
      ilen=nsol*ntchn*1203
      bndcnt=0
      point=0
      do 50 is=1,nsts
      nlg=nlag(is)
      nsch=ncst(is)
      if (nlg.eq.0) go to 50
      do 40 i=1,nlg
      bndcnt=bndcnt+1
      do 30 j=1,nsol
      smat(bndcnt,j)=0.d+00
      pcnt=point
      do 20 nc=1,nsch
      lp=jind(nc,is)
      do 10 k=1,nptmx
      pcnt=pcnt+1
   10 smat(bndcnt,j)=smat(bndcnt,j)+flago(k,lp,bndcnt)*wts(pcnt)*soln
     1 (pcnt,j)
   20 continue
   30 continue
   40 continue
   50 point=point+nsch*nptmx
      do 60 it=1,ntchn
      do 60 jt=1,nflag
   60 smat(jt,it)=-smat(jt,it)
      if (logkey(ops,'print=lam=overlap',.false.,' ')) then
      write (iout,120)
      call matprt (smat,nflag,nsol,nflag,nsol,0,0,0,0,0,0,0)
      endif
      call sgefa (smat(1,ntchn+1),nflag,nflag,ipvt,info)
      do 70 it=1,ntchn
   70 call sgesl (smat(1,ntchn+1),nflag,nflag,ipvt,smat(1,it),0)
      if (logkey(ops,'print=lam=overlap',.false.,' ')) then
      call matprt (smat,nflag,nsol,nflag,nsol,0,0,0,0,0,0,0)
      endif
*
*          now calculate the physical solutions by summing
*          over the lagrange multipliers times the particular
*                                    solutions
*
      uper=ntchn+1
      do 90 i=1,ntchn
      do 90 j=1,msize
      icnt=0
      do 80 k=uper,nsol
      icnt=icnt+1
   80 soln(j,i)=soln(j,i)+smat(icnt,i)*soln(j,k)
   90 continue
      if (logkey(ops,'lam=output=wavefunctions',.false.,' ')) then
      call iosys ('create real "spline vecs'//itoc(ien)//'" on wavefn'
     1 ,ilen,0,0,0)
      do 110 isol=1,nsol
      icnt=0
      do 100 ii=1,nsts
      nsch=ncst(ii)
      do 100 jj=1,nsch
      call bsplin (nptmx,rgs,soln(icnt+1,isol),3,nbk,brk,c,iflag,sc)
      call iosys (wrtspn(ien),1203,spln,0,0)
  100 icnt=icnt+nptmx
  110 continue
      endif
      return
c
  120 format (/,5x,'overlap matrix')
      end
