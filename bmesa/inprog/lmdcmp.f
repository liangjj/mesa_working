*deck lmdcmp.f
c***begin prologue     lmdcmp
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           lmdcmp, link 2702, y(l,m), projections
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            decompose a function into its legendre projections
c***description        a three dimensional function is decomposed into
c***                   its projection onto p(l,m) phi(m) angular
c***                   functions. the angular functions in the phi
c***                   variable are the normalized sin(m*phi) and
c***                   cos(m*phi). it is assumed that the function is
c***                   stored with no blank spaces in any column direction.
c***references         none
c
c***routines called
c***end prologue       lmdcmp
      subroutine lmdcmp (f,plm,phifn,r,flm,scr,mval,nl,nr,nth,nph,
     1                   nm,prnt)
      implicit integer (a-z)
      real *8 f, plm, phifn, flm, scr, r
      logical prnt
      character*2 itoc, tmp
      character*80 title
      dimension f(nr,nth,nph), plm(nth,*), phifn(nph,*)
      dimension scr(nr,nth), flm(nr,*), r(nr)
      dimension mval(nm), nl(nm)
      common /io/ inp, iout
      plmcnt=1
      flmcnt=1
      do 10 mu=1,nm
         call ebc(scr,f,phifn(1,mu),nr*nth,nph,1)
         if (prnt) then
             write(iout,1) mval(mu)
             tmp=itoc(mval(mu))
             len=length(tmp)
         endif
         call ebc(flm(1,flmcnt),scr,plm(1,plmcnt),nr,nth,nl(mu))
         call vmmul(r,flm(1,flmcnt),flm(1,flmcnt),nr,
     1              nl(mu))
         if (prnt) then
             title='p(l,'//tmp(1:len)//') decomposition'
             call prntfm(title,flm(1,flmcnt),nr,nl(mu),nr,
     1                   nl(mu),iout)
         endif                    
         flmcnt=flmcnt+nl(mu)     
         plmcnt=plmcnt+nl(mu)
   10 continue
    1 format(/,5x,'the legendre decomposition matrix for m = ',i3)
      return
      end
