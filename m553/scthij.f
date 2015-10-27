*deck @(#)scthij.f	5.1  11/6/94
      subroutine scthij(mixinv,hess,ci,cj,nbfi,nbfj,leni,lenj,mixi,
     $     mixj,loci,locj,ipt,jpt,isym,iorb,jsym,jorb,tv,htran,
     $     shmo,lshmo,lbhmo,asort,ncors,nmix)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)scthij.f	5.1   11/6/94
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
      dimension hess(nbfi,nbfj),htran(leni,lenj),hmo(2),tv(2),
     $     ci(nbfi,2),cj(nbfj,2),loci(2),locj(2),mixi(2),
     $     mixj(2),mixinv(nbfi,2)
cc
      dimension shmo(*),lshmo(*),lbhmo(*),asort(*)
cc

      common / number / zero,pt5,one,two,four,eight
c
c     this program transforms the hessian from the ao
c     basis to the mo basis
c
c     hess(i,j) hessian in the ao basis
c     ci        mo's for symmetry i
c     cj        mo's for symmetry j
c     mixi      mo pointer vector for orbital i
c     mixj      mo pointer vector for orbital j
c     loci      hessian pointer vector for i-a mixings
c     locj      hessian pointer vector for j-b mixings
c     leni      the number of i-a mixings
c     lenj      the number of j-b mixings
c     nbfi      the number of basis functions in symmetry i
c     nbfj      the number of basis functions in symmetry j
c     hmo       the hessian in the mo basis
c     htran     temporary storage for the half-transformed hessian
c     tv        temporary vector storage used in the transformation
c     iorb      hessian label of orbital i
c     jorb      hessian label of orbital j
c     ipair     pointer array into a triangular matrix
c
cc--------------------
c    transform section
cc--------------------
c
      ipair(i)=(i*(i-1))/2
c
      do 70 i=1,leni
         im=mixi(i)
         if(im.gt.iorb)go to 21
c-----------------
c   positive phase
c-----------------
         do 20 j=1,nbfj
            xx=zero
            do 10 k=1,nbfi
               xx=xx+ci(k,im)*hess(k,j)
 10         continue
            tv(j)=xx
 20      continue
         go to 41
 21      continue
c-----------------
c   negative phase
c-----------------
         do 40 j=1,nbfj
            xx=zero
            do 30 k=1,nbfi
               xx=xx+ci(k,im)*hess(k,j)
 30         continue
            tv(j)=-xx
 40      continue
 41      continue
c------------------------------------------------------
c  finish off the transform for this row of the hessian
c------------------------------------------------------
         do 60 j=1,lenj
            jm=mixj(j)
            xx=zero
            do 50 k=1,nbfj
               xx=xx+tv(k)*cj(k,jm)
 50         continue
c---------------------------------------------
c   adjust the phase for the j orbital mixings
c---------------------------------------------
            if(jm.gt.jorb)go to 55
            htran(i,j)=xx
            go to 60
 55         continue
            htran(i,j)=-xx
 60      continue
c
 70   continue
c
c--------------------------------------------
c   put away the transformed hessian elements
c--------------------------------------------
c
      if(isym.ne.jsym) go to 79
      ilen=mixinv(jorb,iorb)
      jlen=mixinv(iorb,jorb)
      if(ilen.eq.0.or.jlen.eq.0) go to 79
      htran(ilen,jlen)=two*htran(ilen,jlen)
 79   continue
c
      iix=0
c
      do 90 i=1,leni
         ii=(loci(i)-1)*nmix
         do 80 j=1,lenj
c
            iix=iix+1
            shmo(iix)=htran(i,j)
            lshmo(iix)=ii+locj(j)
c
            if(locj(j).ne.loci(i)) then
               iix=iix+1
               shmo(iix)=htran(i,j)
               lshmo(iix)=(locj(j)-1)*nmix+loci(i)
            endif
cc
ccc      hmo(ix)=hmo(ix)+htran(i,j)
cc
 80      continue
 90   continue
c
      call sorter('with bin',asort,asort,0,iix,lshmo,
     $     lbhmo,shmo,0,0,0,.false.)
c
      return
      end