*deck plmtof.f
c***begin prologue     plmtof
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           plm, y(l,m), projections
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            reconstruct a three dimensional function from
c***                   its radial flm(r) projection and p(l,m)*phifn(m)
c***description        the input projections and the angular functions
c***                   are summed to recreate the original three 
c***                   dimensional function. this can be used to make
c***                   a pointwise least squares comparison with the
c***                   exact original function.
c***references         none
c
c***routines called
c***end prologue       plmtof
      subroutine plmtof (fns,fs,plm,phifn,ylm,flm,scr,pt,l,m,nr,nth,nph,
     1                   nang,nonsep,prnt)
      implicit integer (a-z)
      real*8 fns, fs, plm, phifn, ylm, flm, scr, pt
      logical prnt, nonsep
      character*80 title
      dimension fns(nr,nang), fs(nr,nth,nph), plm(nth,*)
      dimension phifn(nph,*), ylm(nang,*), flm(nr,*), scr(*), pt(nr)
      common /io/ inp, iout
      if (nonsep) then
          dimang=nang
          call rzero(fns,nr*nang)
          cntlm=1
          do 10 mu=0,m
             ndim=l-mu+1
             no=2
             if (mu.eq.0) then
                 no=1
             endif
             do 20 count=1,no
                call apbct(fns,flm(1,cntlm),ylm(1,cntlm),nr,
     1                     ndim,nang)
                cntlm=cntlm+ndim
   20        continue
   10     continue
          if (prnt) then
              title='fitted function'
              call prntfm(title,fns,nr,dimang,nr,dimang,iout)
          endif
      else
          dimang=nth*nph
          call rzero(fs,nr*dimang)
          flmcnt=1
          plmcnt=1
          phicnt=1
          do 30 mu=0,m
             ndim=l-mu+1
             no=2
             if (mu.eq.0) then
                 no=1
             endif
             cntscr=1
             do 40 count=1,no
                call ebct(scr(cntscr),flm(1,flmcnt),plm(1,plmcnt),
     1                    nr,ndim,nth)
                flmcnt=flmcnt+ndim
                cntscr=cntscr+nr*nth
   40        continue
             call apbct(fs,scr,phifn(1,phicnt),nr*nth,no,nph)
             plmcnt=plmcnt+ndim
             phicnt=phicnt+no
   30     continue
          if (prnt) then
              title='fitted function'
              call prntfm(title,fs,nr,dimang,nr,dimang,iout)
          endif
      endif
      return
      end

