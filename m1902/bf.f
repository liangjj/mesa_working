*deck @(#)bf.f	5.1  11/6/94
      subroutine bf(ptprim,noprim,ptcont,nocont,start,
     $              nocart,nobf,minmom,maxmom,mintyp,
     $              nx,ny,nz,cont,ex,z,a)
c***begin prologue     bf.f
c***date written       930519  (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           basis set 
c***author             martin, richard (lanl) 
c***source             @(#)bf.f	5.1   11/6/94
c***purpose            returns information from the rwf file
c                      specifying the basis set. 
c***description
c                      the arrays specifying the basis set are
c                      loaded into z(implicitly equivalenced to a)
c    
c
c***references
c
c***routines called
c
c***end prologue       basis.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer a(1)
      real*8 z(1)
c     --- output variables ---
      integer ptprim,noprim,ptcont,nocont,start
      integer nocart,nobf,minmom,maxmom,mintyp
      integer nx,ny,nz
      integer cont,ex
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
c
      common/io/inp,iout
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
