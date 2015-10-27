*deck @(#)mccg1.f	5.1  11/6/94
      subroutine mccg1(nk,nob,cg,
     $     nf35,r,g,buf,lbufso,sg,rabcx,incor)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mccg1.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cc
cc
      dimension g(2)
      dimension cg(nob,nk)
      real*8 buf(lbufso),sg(*),rabcx(*)
      dimension r(nob,nk)
      common /pcpack/ ipkt, nhext, ipkes
c
      common /io/ inp,iout
c
c-----------------------------------------------------------------------
c
c --- description     this routine makes contributions from a block
c                     of 1-symmetry integrals and the corresponding
c                     density matrix to c*g.
c
c --- input
c
c     nk              no of k active orbitals.
c     nob             no of k orbitals, including core and virtual.
c     g(--)           density matrix elements.
c     nf35            fortran no for dataset containing transformed
c                     integrals.
c
c --- working storage
c
c     r(nob,nk)
c
c --- output
c
c     cg(nob,nk)      c*g
c
c-----------------------------------------------------------------------
c
c
c52
c
c
      nnk=nk*(nk+1)/2
      nknk=nk*nk
      nijkl=nnk*(nnk+1)/2
      mk=nk*nnk
c
c  square and scale the density
c
      ix=0
      do 5 i=1,nnk
         ix=ix+i
         g(ix)=g(ix)*2.0
 5    continue
c
      call trtosq(sg,g,nnk,nijkl)
c
      ix=0
      do 10 k=1,nk
         ix=ix+k
         jx=ix
         do 20 l=1,nnk
            sg(jx)=sg(jx)*2.0
            jx=jx+nnk
 20      continue
 10   continue
c
      ix=1
      jx=1
      do 30 i=1,nnk
         call trtosq(g(jx),sg(ix),nk,nnk)
         jx=jx+nknk
         ix=ix+nnk
 30   continue
c
c
      if(incor.eq.0) then
c
c  rabcx is read in buffers
c
         call iosys('rewind abcx on rwf',0,0,0,' ')
c
         nr=nk*nob
         lpass=lbufso/nr
         npass=(nnk-1)/lpass+1
         ix=1
         nni=nnk
c
         do 35 i=1,npass
c
            call iosys('read real abcx from rwf without rewinding',
     $           lbufso,buf,0,' ')
            nm=min(lpass,nni)
            nni=nni-nm
            mk=nm*nk
c
            call apbct(cg,buf,g(ix),nob,mk,nk)
            ix=ix+nk*mk
c
 35      continue
c
      else
c
c  rabcx is held in core
c
         call apbct(cg,rabcx,g,nob,mk,nk)
c
      end if
c
c
c     do 900 i = 1, nk
c 900 write (iout,9000) (cg(j,i),j=1,nob)
c9000 format(' *mccg1 cg '//4(1x,f16.8))
c
c
      return
 6000 write (iout,9001)
 9001 format(//'0****** mccg1',6x,' error reading i-vectors')
      call lnkerr(' ')
      end
