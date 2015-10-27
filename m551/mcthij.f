*deck @(#)mcthij.f	1.1  11/30/90
      subroutine mcthij(mixinv,hess,ci,cj,nbfi,nbfj,leni,lenj,mixi,
     $     mixj,loci,locj,ipt,jpt,isym,iorb,jsym,jorb,tv,htran,hmo,
     $     itflag,nmix)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcthij.f	1.1   11/30/90
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
cmp   extended dummy hess,htran,hmo,tv,ci,cj,loci,locj,
cmp    1 mixi,mixj,mixinv
cc
      dimension hess(nbfi,nbfj),htran(leni,lenj),hmo(2),tv(2),
     $     ci(nbfi,2),cj(nbfj,2),loci(2),locj(2),mixi(2),
     $     mixj(2),mixinv(nbfi,2)
c
      common /io/ inp,iout
c
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
c     write(iout,7799) ilen,jlen,htran(ilen,jlen)
c7799 format('  jlen ilen htran ',2x,2i5,2x,f20.12)
      htran(ilen,jlen)=two*htran(ilen,jlen)
cc    htran(jlen,ilen)=two*htran(jlen,ilen)
 79   continue
c
c     if(iocore.eq.0) go to 200
c
c      write(iout,*)'  nmix   hmo(2) ',nmix,hmo(2)
c
      do 90 i=1,leni
         in=(loci(i)-1)*nmix
         do 80 j=1,lenj
c
            hmo(in+locj(j))=hmo(in+locj(j))+htran(i,j)
c
            if(loci(i).ne.locj(j)) then
               jpt=(locj(j)-1)*nmix+loci(i)
               hmo(jpt)=hmo(jpt)+htran(i,j)
            endif
c
 80      continue
 90   continue
c
      return
      end
