*deck @(#)basis.f	5.1  11/6/94
      subroutine basis(natoms,nbf,nprim,ncont,ntypes,nbtype,lenxyz,
     $                 mxcont,maxl,
     $                 ptprim,noprim,ptcont,nocont,start,
     $                 nocart,nobf,minmom,maxmom,mintyp,
     $                 nx,ny,nz,cont,ex,top,z,a)
c***begin prologue     basis.f
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
      integer natoms,nbf
      integer top
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer a(1)
      real*8 z(1)
c     --- output variables ---
      integer mxcont,maxl
      integer nprim,ncont,ntypes,nbtype,lenxyz
      integer ptprim,noprim,ptcont,nocont,start
      integer nocart,nobf,minmom,maxmom,mintyp
      integer nx,ny,nz
      integer cont,ex
c     --- scratch arrays ---
c     --- local variables ---
      integer iadtwp,wpadti,ngot
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
      ptprim=top
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
      top=wpadti(ex+nprim)
c
c     --- have we enough room ---
      call getscm(0,z(1),ngot,'basis',0)
      if (top.ge.ngot) then
         call lnkerr('not enough core in basis')
      endif
c
c     read in basis set information from readwrite file.
      call iosys('read integer "pointer to primitives" from rwf',
     $           -1,a(ptprim),0,' ')
      call iosys('read integer "number of primitives" from rwf',
     $           -1,a(noprim),0,' ')
      call iosys('read integer "pointer to contraction coefficients"'//
     $           'from rwf',-1,a(ptcont),0,' ')
      call iosys('read integer "number of contraction coefficients"'//
     $           'from rwf',-1,a(nocont),0,' ')
      call iosys('read integer "pointer to first function" from rwf',
     $           -1,a(start),0,' ')
      call iosys('read integer "number of cartesians" from rwf',
     $           -1,a(nocart),0,' ')
      call iosys('read integer "number of pure functions" from rwf',
     $           -1,a(nobf),0,' ')
      call iosys('read integer "minimum momentum" from rwf',-1,
     $           a(minmom),0,' ')
      call iosys('read integer "maximum momentum" from rwf',-1,
     $           a(maxmom),0,' ')
      call iosys('read integer "pointer to cartesians" from rwf',
     $           -1,a(mintyp),0,' ')
      call iosys('read integer "power of x" from rwf',-1,a(nx),0,' ')
      call iosys('read integer "power of y" from rwf',-1,a(ny),0,' ')
      call iosys('read integer "power of z" from rwf',-1,a(nz),0,' ')
      call iosys('read real "contraction coefficients" from rwf',
     $           -1,z(cont),0,' ')
      call iosys('read real exponents from rwf',-1,z(ex),0,' ')
c
c
      return
      end
