*deck @(#)basgrd.f	5.1 11/6/94
      subroutine basgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     nat,nprim,ntypes,nbtype,nnp,ncont,
     $     start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,ops,xyzgrid,mxgrd,xyzpow,scr,rsq,
     $     s,r,ngrid,phi,grad,dograd,mxcont,wts,ngb,gblksiz,
     $     mxgbsiz,mxgblk,maxl)
c***beginprologue     basgrd.f
c***date written       900513   (yymmdd)
c***revision date      11/6/94
c
c***keywords           basis, grid, amplitude
c***author             martin, richard (lanl), RUSSO, thomas (lanl)
c***source             @(#)basgrd.f	5.1  11/6/94
c***description
c     evaluate all of the contracted basis functions on a grid given by
c     xyzgrid.
c
c***references
c
c***routines called
c
c***end prologue       basgrd.f
c
      implicit none
c
c --- input variables ---
      character*(*) ops
      integer nat,nprim,nbtype,nnp,ncont,ntypes,mxcont,mxgrd,mxgbsiz,
     $     mxgblk,nbf
      logical dograd
c --- input arrays (unmodified) ---
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 xyzgrid(mxgrd,3,nat),wts(mxgrd,nat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
      integer ngrid(nat),ngb(nat),gblksiz(mxgblk,nat)
c --- local variables ---
      integer ngbtot,i,iatom,ioff,iblk,ng
      character*3 answer
c --- output arrays ---
      real*8 phi(mxgbsiz,nbf),grad(mxgbsiz,3,nbf)
c --- scratch arrays ---
      real*8 xyzpow(mxgbsiz,3,*),scr(mxgbsiz),s(mxgbsiz,mxcont),
     $     r(mxgbsiz,mxcont),rsq(mxgbsiz)
      integer maxl(nat)
c
      integer inp,iout
      common/io/inp,iout
c
c
      ngbtot=0
      do 1 i=1,nat
         ngbtot=ngbtot+ngb(nat)
 1    continue 
c
c     ----- the main loop -----
c     the plan is: for each atom, determine the amplitude of all basis 
c     functions on the grid for that atom. 
      call iosys('does "basis set on grid" exist on rwf',
     $           0,0,0,answer)
      if(answer.eq.'no') then
         call iosys(
     $        'create real "basis set on grid" on rwf',
     $        mxgbsiz*ngbtot*nbf*nat,0,0,' ')
      endif
      if(dograd) then
         call iosys('does "basis set gradients on grid"'
     $              //' exist on rwf',0,0,0,answer)
         if(answer.eq.'no') then
            call iosys(
     $           'create real "basis set gradients on grid"'
     $           //' on rwf',mxgbsiz*ngbtot*nbf*nat*3,0,0,' ')
         endif
      endif
c
      do 100 iatom=1,nat
         ioff=0
         do 101 iblk=1,ngb(iatom)
            ng=gblksiz(iblk,iatom)
            call bfgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $                 nat,nprim,ntypes,nbtype,nnp,ncont,
     $                 start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $                 nx,ny,nz,xyzgrid(1,1,iatom),mxgrd,xyzpow,
     $                 scr,rsq,s,r,phi,grad,dograd,mxcont,mxgbsiz,
     $                 ng,ioff,maxl)
c           write out this atom's grid stuff
c
            call iosys(
     $           'write real "basis set on grid" '//
     $           'on rwf without rewinding',mxgbsiz*nbf,phi,0,' ')
            if(dograd) then
               call iosys(
     $              'write real "basis set gradients on grid" '//
     $              'on rwf without rewinding',mxgbsiz*nbf*3,grad,0,' ')
            endif
            ioff=ioff+ng
 101     continue 
 100  continue
c
c
      return
      end
