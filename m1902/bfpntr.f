*deck @(#)bfpntr.f	5.1  11/6/94
      subroutine bfpntr(natoms,nprim,ncont,ntypes,nbtype,ptprim,
     $                  noprim,ptcont,nocont,start,nocart,nobf,minmom,
     $                  maxmom,mintyp,nx,ny,nz,cont,ex,need)
c***begin prologue     bfpntr.f
c***date written       930519  (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           basis set 
c***author             martin, richard (lanl) 
c***source             @(#)basis.f	5.1   11/6/94
c***purpose            returns information from the rwf file
c                      specifying the basis set. 
c***description
c                      the arrays specifying the basis set are
c                      loaded into z(implicitly equivalenced to a)
c                      and returned along with pointers which
c                      access them.
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       basis.f
      implicit none
c     --- input variables -----
      integer natoms
      integer need
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
      integer mxcont,maxl
      integer nprim,ncont,ntypes,nbtype,lenxyz
      integer ptprim,noprim,ptcont,nocont,start
      integer nocart,nobf,minmom,maxmom,mintyp
      integer nx,ny,nz
      integer cont,ex
c     --- scratch arrays ---
c     --- local variables ---
      integer iadtwp,wpadti
      integer inp,iout
c
      common/io/inp,iout
c
c     --- get some basic information ---
      call iosys('length of exponents on rwf',nprim,0,0,' ')
      call iosys('length of "contraction coefficients" on rwf',
     $            ncont,0,0,' ')
c
c     ntypes is the total number of basis function types defined
c     in m102.  this total is composed of two sets: the first nbtype
c     are used for actual basis functions( e.g. s,p,d,f,g,h,sp,...).
c     the remainder are used to describe effective core potential
c     projector types (e.g. us,up,ud,uf,ug,s-u,...).
      call iosys('length of "number of pure functions" on rwf',
     $           ntypes,0,0,' ')
      call iosys('read "number basis types" from rwf',
     $           1,nbtype,0,' ')
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
      call iosys('read integer maxcont from rwf',1,mxcont,0,' ')
      call iosys('read integer maxl from rwf',1,maxl,0,' ')
c
c     --- assign pointers for basis set information ---
      ptprim=1
      noprim=ptprim+ntypes*natoms
      ptcont=noprim+ntypes*natoms
      nocont=ptcont+ntypes*natoms
      start=nocont+ntypes*natoms
      nocart=start+ntypes*natoms
      nobf=nocart+ntypes
      minmom=nobf+ntypes
      maxmom=minmom+ntypes
      mintyp=maxmom+ntypes
      nx=mintyp+ntypes
      ny=nx+lenxyz
      nz=ny+lenxyz
      cont=iadtwp(nz+lenxyz)
      ex=cont+ncont
      need=wpadti(ex+nprim)
c
c
      return
      end
