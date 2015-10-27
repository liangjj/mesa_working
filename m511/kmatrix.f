*deck @(#)kmatrix.f	5.3 4/18/95
      subroutine kmatrix(values,d,dlast,nbf,nnp,kmat,nexch,grdwts,
     $     ndmat,nat,mxgrd,exc,slater,becke,lyp,vwn,calce,calc,
     $     dmcut,dencut,eexch,prteexch,
     $     mxgbsiz,valuesi,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,ian,xyzgrid,charge,maxl,dograd,bigl,
     $     vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,
     $     rmax,lmax,nomega,nradial,grdtyp,adjust,minesz)
c
c***begin prologue     kmatrix.f
c***date written       930521   (yymmdd)  
c***revision date      4/18/95      
c
c***keywords           k matrix, coulomb matrix 
c***author             RUSSO, thomas (lanl)
c***source             @(#)kmatrix.f	5.3   4/18/95
c***purpose            to form the k matrix given basis and gradients
c                      on grid, and density matx
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       kmatrix.f
      implicit none
c     --- input variables ---
      integer nbf,nnp,nexch,ndmat,nat,mxgrd,mxgbsiz
      integer nprim,ntypes,nbtype,ncont,mxcont
      integer rmax,lmax,nomega,nradial
      integer minesz
      character*8 calc
      real*8 dencut,dmcut
      logical dograd,adjust
c     --- input arrays (unmodified) ---
      integer ian(nat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 d(nnp,ndmat),dlast(nnp,ndmat)
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 xyzgrid(mxgrd,3)
      character*(*) grdtyp(nat)
c     --- input arrays (scratch) ---
      real*8 values(*),grdwts(mxgrd)
      real*8 vwts(mxgrd),rnuc(nat,nat),amu(nat,nat)
      real*8 pwtx(nat),rr(nat),radii(nat),akl(nat,nat)
      integer ptrad(rmax)
      integer valuesi(*)
      integer maxl(nat)
c     --- output arrays ---
      real*8 kmat(nnp,nexch)
      real*8 charge(nat,ndmat)
c     --- output variables ---
      real*8 exc,eexch
c     --- local variables ---
      integer inp,iout,dengrida,dengridb,dengrada,dengradb,ga,gb,gab,
     $        fout,therest,phi,grad,scr,scr2,i,queue,tea,
     $        dtmp,tmpgwt,nzptrs,bigl,ngbtot
      integer phibar,gradbar
      logical slater,becke,lyp,vwn,calce,prteexch,dorhoup,called
      integer wpadti,iadtwp
      character*3 answer
c
c
      common /io/ inp,iout
      save called
      data called/.false./
c
c     --- initialize the output array
      call rzero(kmat,nnp*nexch)
      exc=0.0d0
      eexch=0.0d0
c
c     --- allocate some space
      if (calc .eq. 'open') then 
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
            fout=queue+minesz*3
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
            queue=ga+mxgbsiz
            fout=queue+minesz*3
         else
            dengrada=1
            ga=1
            queue=1
            fout=dengrida+mxgbsiz
         endif
      endif
      phi=fout+mxgbsiz*5
      if (lyp) phi=phi+mxgbsiz
      phibar=phi+mxgbsiz*nbf
      grad=phibar+nbf
      gradbar=grad+mxgbsiz*3*nbf
      scr=gradbar+nbf
      if(calc.eq.'open') then
         scr=gradbar+3*nbf
      endif
c     use scr space for strip mining density formation too. (takes minesz*nbf)
      scr2=scr+max(mxgbsiz,minesz*nbf)
      tea=scr2+minesz
      tmpgwt=tea+minesz*nbf
      nzptrs=wpadti(tmpgwt+mxgbsiz)
      therest=iadtwp(nzptrs+mxgbsiz)
c
c     --- create disk files for the density
      if(calc.eq.'closed') then
         if(.not.called) then
c           first time through work with full density.
            call vmove(values(dtmp),d,nnp)
            dorhoup=.false.
            call iosys('does "closed density" exist on rwf',
     $                  0,0,0,answer)
            ngbtot=1+(mxgrd/mxgbsiz)
            if(answer.eq.'no') then
               call iosys('create real "closed density" on rwf',
     $                     mxgbsiz*ngbtot*nat,0,0,' ')
               if(becke.or.lyp) then
                  call iosys('create real "closed density gradient" '
     $                     //' on rwf',mxgbsiz*ngbtot*nat*3,0,0,' ')
               endif
            endif
         else
c           subsequent entries work with density difference.
            call vsub(values(dtmp),d,dlast,nnp)
            dorhoup=.true.
         endif
c
         call kmtxguts(values(dengrida),
     $     values(dengrada),values(ga),
     $     values(fout),values(therest),values(phi),
     $     values(grad),values(scr),slater,becke,lyp,vwn,
     $     values(dtmp),nbf,nnp,kmat,nexch,grdwts,ndmat,nat,mxgrd,calce,
     $     exc,dmcut,dencut,calc,valuesi(nzptrs),values(tmpgwt),
     $     eexch,prteexch, values(queue),values(tea),
     $     dorhoup,
     $     mxgbsiz,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,ian,xyzgrid,charge,maxl,dograd,bigl,
     $     vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,
     $     rmax,lmax,nomega,nradial,grdtyp,adjust,minesz,values(phibar),
     $     values(gradbar))
      else 
         if (.not.called) then
            call vmove(values(dtmp),d,nnp*2)
            dorhoup=.false.
            call iosys('does "open density" exist on rwf',
     $                  0,0,0,answer)
            if(answer.eq.'no') then
               ngbtot=1+(mxgrd/mxgbsiz)
               call iosys('create real "open density" on rwf',
     $                     2*mxgbsiz*ngbtot*nat,0,0,' ')
               if(becke.or.lyp) then
                  call iosys('create real "open density gradient" '
     $                     //' on rwf',2*mxgbsiz*ngbtot*nat*3,0,0,' ')
               endif
            endif
         else
            call vsub(values(dtmp),d,dlast,nnp*2)
            dorhoup=.true.
         endif
         call kmopen(values(dengrida),values(dengridb),
     $     values(dengrada),values(dengradb),values(ga),values(gb),
     $     values(gab),values(fout),values(therest),values(phi),
     $     values(grad),values(scr),values(scr2),slater,
     $     becke,lyp,vwn,values(dtmp),nbf,nnp,kmat,nexch,grdwts,ndmat,
     $     nat,mxgrd,
     $     calce,exc,dmcut,dencut,calc,valuesi(nzptrs),
     $     values(tmpgwt),eexch,prteexch,values(queue),values(tea),
     $     dorhoup,mxgbsiz,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,ian,xyzgrid,charge,maxl,dograd,bigl,
     $     vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,
     $     rmax,lmax,nomega,nradial,grdtyp,adjust,minesz,values(phibar),
     $     values(gradbar))
      endif
c      
c
      do 10 i=1,nexch
         call vneg(kmat(1,i),kmat(1,i),nnp)
 10   continue 
      called=.true.
c
c
      return
      end
