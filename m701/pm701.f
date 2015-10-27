*deck @(#)pm701.f	5.1  11/6/94
      subroutine pm701(z,a)
c***begin prologue     pm701.f
c***date written       870629   (yymmdd)
c***revision date      11/6/94
c
c   22 december 1987   bhl at brl
c      scaled mcscf ao core density
c      c..bhl..scale..density
c
c   21 november 1987   pws at lanl
c      adding capability for mcscf gradients.
c
c***keywords           hf density matrices, primitive basis density
c***author             saxe, paul (lanl)
c***source             @(#)pm701.f	5.1   11/6/94
c
c***purpose            to transform the hartree-fock one-particle density
c     matrices to the primitive basis and organize them effectiveley.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       pm701.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a(*)
      real*8 z(*)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer nat,nbf,npf,ntypes,nbtype
      integer nnp,nnprim,prtoao,dpr,t1,dao,t2,top
      integer ngot
      integer wpadti,iadtwp
      character*4096 ops
      character*8 prtflg
      logical prnt
c
      data prnt/.true./
      save prnt
c
      common /io/     inp,iout
c
 1000 format(' m701: transform one-electron density matrices ')
c
c     --- recover the options string ---
      call iosys('read character options from rwf',-1,0,0,ops)
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
      if(prnt) then
         write(iout,1000)
      endif
c
c     --- get lengths, dimensions, etc. needed for core allocation 
      call iosys('read integer "number of atoms" from rwf',1,nat,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbf,0,' ')
      call iosys('read integer "number of primitive functions" '//
     $           'from rwf',1,npf,0,' ')
      call iosys('length of "number of pure functions" on rwf',
     $           ntypes,0,0,' ')
      call iosys('read integer "number basis types" from rwf',
     $           1,nbtype,0,' ')
c
      nnp=nbf*(nbf+1)/2
      nnprim=(npf+1)*npf/2
c
c     --- divide up core ---
      prtoao=iadtwp(1)
      dpr=prtoao+npf*nbf
      t1=dpr
      dao=dpr+npf**2
      t2=dao
c
      top=wpadti(t2+max(nbf*npf,nnprim))
c
      call getscm(top,a,ngot,'m701 transformation',0)
c
c     --- run through the density matrices, transforming each ---
      call tr1dm(z(t1),z(t2),z(dao),z(prtoao),z(dpr),nbf,nnp,npf,nnprim,
     $           ops)
c
c     --- and exit gracefully ---
c
      call chainx(0)
c
c
      stop
      end
