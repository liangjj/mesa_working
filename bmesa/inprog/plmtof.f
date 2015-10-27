c $Header: plmtof.f,v 1.2 92/12/31 14:44:00 bis Exp $
*deck plmtof.f
c***begin prologue     plmtof
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           plm, link 2702, y(l,m), projections
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            reconstruct a three dimensional function from
c***                   its radial y(l,m) projection and p(l,m)*phifn(m)
c***description        the input projections and the angular functions
c***                   are summed to recreate the original three 
c***                   dimensional function. this can be used to make
c***                   a pointwise least squares comparison with the
c***                   exact original function.
c***references         none
c
c***routines called
c***end prologue       plmtof
      subroutine plmtof (f,plm,phifn,flm,scr,lmax,m,nr,nth,nph,
     1                   num,prnt)
      implicit integer (a-z)
      real*8 f, plm, phifn, flm, scr
      logical prnt
      character*2 itoc
      character*80 title
      dimension f(nr,nth,nph,2), plm(nth,m:lmax), phifn(nph,num)
      dimension flm(nr,m:lmax,num), scr(nr,nth,num)
      common /io/ inp, iout
      dim=lmax-m+1
      do 10 i=1,num
         call ebct(scr(1,1,i),flm(1,m,i),plm(1,m),nr,dim,nth)
   10 continue
      call apbct(f(1,1,1,2),scr,phifn,nr*nth,num,nph)
      if (prnt) then
          title='exact and fitted function for m = '//itoc(m)
          prd=nr*nth*nph
          call prntfm(title,f,prd,2,prd,2,iout)
      endif
      return
      end
