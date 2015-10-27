*deck %W%  %G%
      subroutine formjk(xs,x2,nbfx,nnpx,
     $                  d,nnp,nbf,jmat,kmat,ncoul,
     $                  nexch,ndmat,dsq,
     $                  t1,t2,t3,t4,t5,t6)
c
c***begin prologue     formjk
c***date written       870521   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           j and k matrices, coulomb matrices,
c                      exchange matrices
c***author             martin, richard (lanl)
c***source             %W%   %G%
c
c***purpose            to form j and k matrices given auxiliary integrals and
c                      density matrices.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       formjk
c
      implicit integer (a-z)
c
c     ----- input arguments,unmodified -----
      integer nbfx,nnpx,nbf,nnp
      integer ncoul,nexch
      real*8 xs(nnp,nbfx),x2(nnpx),d(nnp,ndmat)
c
c     ----- output arrays -----
      real*8 jmat(nnp,ncoul),kmat(nnp,nexch)
c
c     ----- scratch arrays -----
      real*8 dsq(nbf,nbf,nexch)
      real*8 t1(nbf*nbf,nbfx)
      real*8 t2(nbfx,nbfx)
      real*8 t3(nbf*nbf,nbfx)
      real*8 t4(nbf*nbf,nbfx)
      real*8 t5(nbf,nbf,ncoul)
      real*8 t6(nbf,nbf,nexch)
c
      real*8 zero,one
      logical called,debug,debug2
      real*8 jnk
      parameter (zero=0.0d+0,one=1.0d+0,debug=.false.,debug2=.true.)
c
      data called /.false./
      save called
c
      common /io/ inp,iout
c
c
      foo=0
      if(.not.called) then
         foo=1
         called=.true.
      endif
c
c
      call rzero(jmat,nnp*ncoul)
      call rzero(kmat,nnp*nexch)
      call rzero(t5,nbf*nbf*ncoul)
      call rzero(t6,nbf*nbf*nexch)
c
c     ----- square up the density matrices -----
c
      nbfsq=nbf*nbf
      do 60 i=1,nexch
         call trtosq(dsq(1,1,i),d(1,i),nbf,nnp)
         if(debug) then
           write(iout,*) 'density matrix'
           call print(d(1,i),nnp,nbf,iout)
         endif
   60 continue
c
c      ----- square up the expansion overlap matrices -----
      do 70 fx=1,nbfx
         call trtosq(t1(1,fx),xs(1,fx),nbf,nnp)
         if(debug) then
            write(iout,*) 'expansion overlap:,function ',fx
            call print(xs(1,fx),nnp,nbf,iout)
         endif
   70 continue
c
c     ----- square up the auxiliary two-electron integrals -----
      call trtosq(t2,x2,nbfx,nnpx)
      if(debug) then
         write(iout,*) 'auxiliary 2e'
         call print(x2,nnpx,nbfx,iout)
      endif
c
c --- INCREDIBLE KLUDGY DEBUGGING TRASH----
c Try forming the "two-electron integrals" using the auxiliary expansion
c but do it without having to store nbf**4 matrices, i.e. do it in slow
c trashy fortran without blas calls and only if debug2 is set.
c
      if (debug2 .and. foo.eq.1 ) then
         write(iout,*)"m505: Fitted 2-e integrals"
         do 6901 i69=1,nbf
            do 6902 j69=1,nbf
               boff=i69+(j69-1)*nbf
               do 6903 k69=1,nbf
                  do 6904 l69=1,nbf
                     jnk=0.0
                     koff=k69+(l69-1)*nbf
                     do 6905 mp=1,nbfx
                        do 6906 m=1,nbfx
                           jnk=jnk+t1(boff,m)*t2(m,mp)*t1(koff,mp)
 6906                   continue 
 6905                continue 
                     if ( abs(jnk) .gt. 1d-16 ) then
                        write(iout,*)i69," ",j69," ",k69," ",l69,
     $                       " ",jnk
                     endif
 6904             continue 
 6903          continue 
 6902       continue 
 6901    continue 
      endif
c
c     ----- form the coulomb matrices -----
c     contract density with xs; i.e. xs(ls,x)*p(ls)
c     work oin this so that i just use the triangle stuff.
c     is density halved on diagonal, then multiply by two at end.
      do 80 shell=1,ncoul
         call sgemv('t',nbfsq,nbfx,one,t1(1,1),nbfsq,dsq(1,1,shell),1,
     $               zero,t3(1,1),1)
         if(debug) then
            write(iout,*) 't3',(t3(i,1),i=1,nbfx)
         endif
         call sgemv('n',nbfx,nbfx,one,t2(1,1),nbfx,t3(1,1),1,
     $               zero,t4(1,1),1)
         if(debug) then
            write(iout,*) 't4',(t4(i,1),i=1,nbfx)
         endif
         call sgemv('n',nbfsq,nbfx,one,t1(1,1),nbfsq,t4(1,1),1,
     $               zero,t5(1,1,shell),1)
         if(debug) then
            write(iout,*) 't5',((t5(i,j,shell),i=1,nbf),j=1,nbf)
         endif
   80 continue
c
c     ----- form exchange matrices -----
      do 130 shell=1,nexch
         do 110 fx=1,nbfx
            call sgemm('n','n',nbf,nbf,nbf,one,t1(1,fx),nbf,
     $                 dsq(1,1,shell),nbf,zero,t3(1,fx),nbf)        
            if(debug) then
               write(iout,*) 't3,fx',fx,(t3(i,fx),i=1,nbf*nbf)
            endif
  110    continue
         call sgemm('n','n',nbfsq,nbfx,nbfx,one,t3(1,1),nbfsq,
     $              t2(1,1),nbfx,zero,t4(1,1),nbfsq)        
         do 111 fx=1,nbfx
            if(debug) then
               write(iout,*) 't4,fx',fx,(t4(i,fx),i=1,nbf*nbf)
            endif
  111    continue
         do 120 fx=1,nbfx
            call sgemm('n','t',nbf,nbf,nbf,one,t4(1,fx),nbf,
     #                 t1(1,fx),nbf,one,t6(1,1,shell),nbf)
            if(debug) then
               write(iout,*) 't6,fx',fx,
     $                      ((t6(i,j,shell),i=1,nbf),j=1,nbf)
            endif
  120    continue
  130 continue
c
c
c     ----- return jmat,kmat as lower triangles -----
      do 140 shell=1,ncoul
         call sqtotr(jmat(1,shell),t5(1,1,shell),nbf,nnp)
  140 continue
      do 150 shell=1,nexch
         call sqtotr(kmat(1,shell),t6(1,1,shell),nbf,nnp)
  150 continue
c
c
      return
      end
