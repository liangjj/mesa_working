c $Header: decomp.f,v 1.2 92/12/24 15:29:55 bis Exp $
*deck decomp.f
c***begin prologue     decomp
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           decomp, link 2702, y(l,m), projections
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
c***end prologue       decomp
      subroutine decomp (f,plm,phifn,flm,scr,lmax,m,nr,nth,nph,num,prnt)
      implicit integer (a-z)
      real *8 f, plm, phifn, flm, scr
      logical prnt
      character*2 itoc, tmp
      character*80 title
      dimension f(nr,nth,nph), plm(nth,m:lmax), phifn(nph,num)
      dimension flm(nr,m:lmax,num), scr(nr,nth,num)
      common /io/ inp, iout
c     the packing of the nr,nth labels allows us to consider f
c     a matrix having (nr*nth) rows. that is the column is tightly
c     packed.
      dim=lmax-m+1
      call ebc(scr,f,phifn,nr*nth,nph,num)
      do 10 i=1,num
         call ebc(flm(1,m,i),scr(1,1,i),plm(1,m),nr,nth,dim)
   10 continue
      if (prnt) then
          do 20 i=1,num
             tmp=itoc(m)
             len=length(itoc(m))
             title='p(l,'//tmp(1:len)//')-'//itoc(i)//' decomposition'
             call prntfm(title,flm(1,m,i),nr,dim,nr,dim,iout)
   20     continue
      endif
      return
      end
