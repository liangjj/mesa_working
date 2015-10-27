*deck @(#)dxcint.f	5.3 4/17/95
      subroutine dxcint(values,d,nbf,nnp,dkay,nexch,grdwts,
     $     ndmat,nat,mxgrd,exc,slater,becke,lyp,vwn,calce,calc,
     $     dmcut,dencut,eexch,prteexch,
     $     mxgblk,mxgbsiz,valuesi,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,xyzgrid,charge,bigl,
     $     ops,nderiv,grade,gradwts,ian,ptrad,radshls,vwts,rnuc,amu,
     $     pwtx,rr,radii,akl,rmax,lmax,nomega,nradial,grdtyp,adjust,
     $     minesz)
c***begin prologue     dxcint.f
c***date written       930521   (yymmdd)  
c***revision date      4/17/95
c   may 13, 1994       rlm at lanl
c      modifying kmatrix(m511) to do derivatives.
c
c***keywords           k matrix, coulomb matrix 
c***author             RUSSO, thomas (lanl)
c***source             @(#)dxcint.f	5.3 4/17/95
c***purpose            to form the k matrix given basis and gradients
c                      on grid, and density matx
c***description
c     
c    
c
c***references         Johnson, Gill and Pople, J. Chem. Phys.,98,5612(1993)
c
c***routines called
c
c***end prologue       dxcint.f
      implicit none
c     --- input variables -----
      integer nbf,nnp,nexch,ndmat,nat,mxgrd,mxgbsiz,mxgblk
      integer nprim,ntypes,nbtype,ncont,mxcont,bigl
      integer nderiv
      integer rmax,lmax,nomega,nradial,minesz
      character*(*) ops
      character*8 calc
      real*8 dencut,dmcut
c     --- input arrays (unmodified) ---
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
      integer ian(nat)
      real*8 d(nnp,ndmat)
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 xyzgrid(mxgrd,3,nat),gradwts(mxgrd,3,nat,nat)
      logical adjust
      character*(*) grdtyp(nat)
c     --- input arrays (scratch) ---
      integer valuesi(*)
      integer ptrad(rmax,nat),radshls(nat)
      real*8 values(*),grdwts(*)
      real*8 vwts(mxgrd),rnuc(nat,nat),amu(nat,nat),pwtx(nat),
     $       rr(nat),radii(nat),akl(nat,nat)
c     --- output arrays ---
      real*8 dkay(nbf,nbf,3,nexch)
      real*8 charge(nat,ndmat)
      real*8 grade(3,nat)
c     --- output variables ---
      real*8 exc,eexch
c     --- scratch arrays ---
c     --- local variables ---
c     --- local variables ---
      integer inp,iout,dengrida,dengridb,dengrada,dengradb,ga,gb,gab,
     $        fout,therest,phi,grad,hess,scr,scr2,queue,tea,kay,
     $        ktmp,dtmp,tmpgwt,nzptrs,maxl,accum,drhoaa
      integer i,coord
      logical slater,becke,lyp,vwn,calce,prteexch
      logical debug
      integer wpadti,iadtwp
c
      parameter (debug=.false.)
c
      common /io/ inp,iout
c
 1000 format(5x,'the derivative k-matrix')
 1010 format(5x,'coordinate:',i3)
 1020 format(5x,'exchange-correlation energy:',1x,f15.6)
c
c     --- initialize array which will hold derivatives of the K matrix.
      call rzero(dkay,nbf*nbf*3*nexch)
      call rzero(grade,3*nat)
c
c     --- initialize variables
      exc=0.0d+00
      eexch=0.0d+00
c
c     --- evaluate derviatives of the K matrix.
      if (calc .eq. 'open' .or. calc .eq. 'general' ) then 
         dtmp=1
         dengrida=dtmp+nnp*2
         dengridb=dengrida+mxgbsiz
         if (becke .or. lyp) then
            dengrada=dengridb+mxgbsiz
            dengradb=dengrada+3*mxgbsiz
            ga=dengradb+3*mxgbsiz
            gb=ga+mxgbsiz
            if (lyp) then
               gab=gb+mxgbsiz
            else
               gab=gb
            endif
            queue=gab+mxgbsiz
            fout=queue+mxgbsiz*3
         else
            dengrada=1
            dengradb=1
            ga=1
            gb=1
            gab=1
            queue=1
            fout=dengridb+mxgbsiz
         endif
      else if (calc .eq. 'closed') then
         dtmp=1
         dengrida=dtmp+nnp
         if (becke .or. lyp) then
            dengrada=dengrida+mxgbsiz
            ga=dengrada+3*mxgbsiz
c           queue isn't used anymore
            queue=ga
            fout=queue+mxgbsiz*3
         else
            dengrada=1
            ga=1
            queue=1
            fout=dengrida+mxgbsiz
         endif
      endif
      phi=fout+mxgbsiz*5
      if (lyp) phi=phi+mxgbsiz
      grad=phi+mxgbsiz*nbf
      hess=grad+mxgbsiz*3*nbf
      scr=hess+mxgbsiz*nbf*6
c     use scr for gridden scratch too.
      tea=scr+max(mxgbsiz,minesz*nbf)
      tmpgwt=tea+mxgbsiz*nbf
      accum=tmpgwt+mxgbsiz
      nzptrs=wpadti(accum+nbf*3)
      maxl=nzptrs+mxgbsiz
      therest=iadtwp(maxl+nat)
      if(calc.eq.'closed') then
c         call vmove(values(dtmp),d,nnp)
c
         call dxclos(values(dengrida),
     $     values(dengrada),values(ga),
     $     values(fout),values(therest),values(phi),
     $     values(grad),values(hess),values(scr),slater,becke,lyp,vwn,
     $     d,nbf,nnp,dkay,grdwts,gradwts,ndmat,nat,mxgrd,
     $     calce,exc,dmcut,dencut,calc,valuesi(nzptrs),
     $     values(tmpgwt),eexch,prteexch,values(queue),values(tea),
     $     mxgbsiz,mxgblk,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,xyzgrid,valuesi(maxl),bigl,charge,grade,
     $     values(accum),ian,ptrad,radshls,vwts,rnuc,amu,pwtx,rr,radii,
     $     akl,rmax,lmax,nomega,nradial,grdtyp,adjust,minesz)
      else 
         call vmove(values(dtmp),d,nnp*2)
c
c$$$         call dxopen(values(dengrida),values(dengridb),
c$$$     $     values(dengrada),values(dengradb),values(ga),values(gb),
c$$$     $     values(gab),values(fout),values(therest),values(phi),
c$$$     $     values(grad),values(scr),values(scr2),slater,
c$$$     $     becke,lyp,vwn,values(dtmp),nbf,nnp,dkay,nexch,grdwts,
c$$$     $     ndmat,nat,mxgrd,
c$$$     $     calce,exc,ngrid,dmcut,dencut,calc,valuesi(nzptrs),
c$$$     $     values(tmpgwt),eexch,prteexch,values(queue),values(tea),
c$$$     $     values(kay),values(ktmp),ngb,gblksiz,mxgbsiz,mxgblk,
c$$$     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
c$$$     $     mxcont,nprim,ntypes,nbtype,ncont,
c$$$     $     start,nocart,nobf,maxmom,minmom,mintyp,
c$$$     $     nx,ny,nz,xyzgrid,charge,valuesi(maxl),bigl)
         call lnkerr('Gack, no open shells yet') 
      endif
      if (debug) then
         write(iout,1020) exc
      endif
c
c
      return
      end
