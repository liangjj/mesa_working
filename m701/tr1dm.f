*deck @(#)tr1dm.f	5.1  11/6/94
      subroutine tr1dm(t1,t2,dao,prtoao,dpr,nbf,nnp,npf,nnprim,ops)
c***begin prologue     tr1dm.f
c***date written       860820  (yymmdd)
c***revision date      11/6/94
c
c 21 november 1987     pws at lanl
c    adding sections for mcscf gradients.
c
c 29 june 1986   pws at lanl
c    moving from m712 to m701, and changing to handle one density
c    matrix at a time to save on core space.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)tr1dm.f	5.1   11/6/94
c***purpose            to transform the hf density matrices from the ao
c     to the primitive basis
c***description
c
c n.b.  this routine is set up so that 't1' and 'dpr' may be implicitly
c       equivalenced in the call, as may 't2' and 'dao'
c
c***references
c***routines called
c***end prologue       tr1dm.f
      implicit none
c     --- input variables -----
      integer nbf,nnp,npf,nnprim
      character*(*) ops
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 prtoao(npf,nbf)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 dao(nnp)
      real*8 dpr(npf,npf)
      real*8 t1(nbf,nbf)
c     --- t2 also holds a upper triangle 'nnprim' long
      real*8 t2(nbf,npf)
c     --- local variables ---
      integer inp,iout
      integer i,ndmat,dmat
      character*3 answer
      logical logkey
c
      common /io/     inp,iout
c
 1000 format(5x,'the primitive-to-ao transformation matrix')
 1010 format (5x,'the ao hf density matrix for shell ',i2)
 1020 format (5x,'the primitive hf density matrix for shell ',
     $           i2)
 1030 format (5x,'the ao mcscf core density matrix')
 1040 format (5x,'the primitive mcscf core density matrix')
 2030 format (5x,'the ao mcscf active density matrix')
 2040 format (5x,'the primitive mcscf active density matrix')
c
c     --- get the primitive-to-ao transformation matrix ---
      call iosys('read real t(prim,cont) from rwf',
     $            npf*nbf,prtoao,0,' ')
      if (logkey(ops,'m701=print=primitive-to-ao-matrix',.false.,' '))
     $     then
         write (iout,1000)
         call matout(prtoao,npf,nbf,npf,nbf,iout)
      end if
c
c     --- transform all the one-particle density matrices to the 
c         primitive basis
      if (logkey(ops,'hf',.false.,' ').or.
     $    logkey(ops,'dft',.false.,' ')) then
         call iosys('read integer "number of hf density matrices" '//
     $              'from rwf',1,ndmat,0,' ')
         call iosys('does "hf primitive density" exist on rwf',
     $               0,0,0,answer)
         if (answer.eq.'no') then
            call iosys('create real "hf primitive density" on rwf',
     $                  nnprim*ndmat,0,0,' ')
         end if
         do 1 dmat=1,ndmat
            call iosys('read real "hf density matrix" from rwf',nnp,
     $                  dao,(dmat-1)*nnp,' ')
            if (logkey(ops,'m701=print=ao-density-matrix',.false.,' '))
     $           then
               write (iout,1010) dmat
               call print(dao,nnp,nbf,iout)
            end if
c
            call trtosq(t1,dao,nbf,nnp)
            call ebct(t2,t1,prtoao,nbf,nbf,npf)
            call ebc(dpr,prtoao,t2,npf,nbf,npf)
            call sqtotr(t2,dpr,npf,nnprim)
            call iosys('write real "hf primitive density" to rwf',
     $                  nnprim,t2,(dmat-1)*nnprim,' ')
            if (logkey(ops,'m701=print=primitive-density-matrix',
     $           .false.,' ')) then
               write (iout,1020) dmat
               call print(t2,nnprim,npf,iout)
            end if
 1       continue
      end if
c
c     --- check for mcscf density matrices ---
      if (logkey(ops,'mcscf',.false.,' ')) then
c
         call iosys('read real "mcscf ao core density" from rwf',nnp,
     $               dao,0,' ')
c
c        --- scale density
         do 1029 i=1,nnp
            dao(i)=dao(i)*.5d0
 1029    continue
         if (logkey(ops,'m701=print=ao-density-matrix=mcscf',
     $        .false.,' ')) then
            write (iout,1030)
            call print(dao,nnp,nbf,iout)
         end if
c
         call trtosq(t1,dao,nbf,nnp)
         call ebct(t2,t1,prtoao,nbf,nbf,npf)
         call ebc(dpr,prtoao,t2,npf,nbf,npf)
         call sqtotr(t2,dpr,npf,nnprim)
         call iosys('write real "mcscf primitive core density" to rwf',
     $               nnprim,t2,0,' ')
         if (logkey(ops,'m701=print=primitive-density-matrix=mcscf',
     $        .false.,' ')) then
            write (iout,1040)
            call print(t2,nnprim,npf,iout)
         end if
c
         call iosys('read real "mcscf ao active density" from rwf',nnp,
     $               dao,0,' ')
         if (logkey(ops,'m701=print=ao-density-matrix=mcscf',
     $        .false.,' ')) then
            write (iout,2030)
            call print(dao,nnp,nbf,iout)
         end if
c
         call trtosq(t1,dao,nbf,nnp)
         call ebct(t2,t1,prtoao,nbf,nbf,npf)
         call ebc(dpr,prtoao,t2,npf,nbf,npf)
         call sqtotr(t2,dpr,npf,nnprim)
         call iosys('write real "mcscf primitive active density" '//
     $              'to rwf',nnprim,t2,0,' ')
         if (logkey(ops,'m701=print=primitive-density-matrix=mcscf',
     $        .false.,' ')) then
            write (iout,2040)
            call print(t2,nnprim,npf,iout)
         end if
      end if
c
c
      return
      end
