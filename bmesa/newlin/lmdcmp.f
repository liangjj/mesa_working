c $Header: lmdcmp.f,v 1.2 92/12/24 15:29:55 bis Exp $
*deck lmdcmp.f
c***begin prologue     lmdcmp
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           lmdcmp, link 2702, y(l,m), projections
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            decompose a three dimensional function 
c***                   into its y(l,m) projections.
c***description        a function labelled space co-ordinates
c***                   is decomposed into its projection onto 
c***                   p(l,m) phi(m) angular functions. the angular
c***                   functions in the phi variable are the 
c***                   normalized sin(m*phi) or cos(m*phi). it is
c***                   assumed that the function is stored with no
c***                   blank spaces in any column direction.
c***references         none
c
c***routines called
c***end prologue       lmdcmp
      subroutine lmdcmp (f,plm,phifn,flm,scr,lmax,nr,nth,nph,prnt)
      implicit integer (a-z)
      real *8 f, plm, phifn, flm, scr
      logical prnt
      character*80 title
      dimension f(nr,nth,nph), plm(nth,lmax), phifn(nph)
      dimension flm(nr,lmax), scr(nr,nth)
      common /io/ inp, iout
c     the packing of the nr,nth labels allows us to consider f
c     a matrix having (nr*nth) rows. that is the column is tightly
c     packed.
c     project onto this specific phifn
      call ebc(scr,f,phifn,nr*nth,nph,1)
c     follow with projection onto the appropriate p(l,m)
      call ebc(flm,scr,plm,nr,nth,lmax)
      if (prnt) then
          title='p(l,m) decomposition'
          call prntfm(title,flm,nr,lmax,nr,lmax,iout)
      endif
      return
      end

