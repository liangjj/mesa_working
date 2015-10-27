*deck @(#)mcthii.f	1.1  11/30/90
      subroutine mcthii(hess,ci,nbfi,leni,mixi,
     $     loci,ipt,isym,iorb,tv,htran,hmo,itflag,nmix)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcthii.f	1.1   11/30/90
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
cmp   extended dummy hess,htran,hmo,tv,ci,loci,mixi
cc
      dimension hess(2),htran(2),hmo(2),tv(2),
     $     ci(nbfi,2),loci(2),mixi(2)
c
      common / number / zero,pt5,one,two,four,eight
c     this program transforms a diagonal block of
c     the hessian from the ao basis to the mo basis
c
c     hess(i,i) hessian in the ao basis
c     ci        mo's for symmetry i
c     mixi      mo pointer vector for orbital i
c     loci      hessian pointer vector for i-a mixings
c     leni      the number of i-a mixings
c     nbfi      the number of basis functions in symmetry i
c     hmo       the hessian in the mo basis
c     htran     temporary storage for the half-transformed hessian
c     tv        temporary vector storage used in the transformation
c     isym      symmetry of orbital i
c     iorb      hessian label of orbital i
c
c     note --   the diagonal block is folded before it is transformed
c
c---------------------------c
c   fold the diagonal block
c---------------------------c
c
c     call printm(ci,nbfi,3)
c     write(iout,9901)
c     call printm(hess,nbfi,3)
c9901 format('  mcthii  hess ')
      ix=0
      do 1 i=1,nbfi
         jx=i
         js=ix+1
         je=ix+i
         do 2 j=js,je
            xx=hess(j)+hess(jx)
            hess(j)=xx
            hess(jx)=xx
            jx=jx+nbfi
 2       continue
         ix=ix+nbfi
    1 continue
c     write(iout,9902)
c9902 format('  hess  mcthii folded ')
c     call printm(hess,nbfi,3)
c
cc--------------------
c    transform section
cc--------------------
c
c
      iix=0
      do 70 i=1,leni
         im=mixi(i)
         if(im.gt.iorb)go to 21
c-----------------
c   positive phase
c-----------------
         ix=0
         do 20 j=1,nbfi
            xx=zero
            do 10 k=1,nbfi
               ix=ix+1
               xx=xx+ci(k,im)*hess(ix)
 10         continue
            tv(j)=xx
 20      continue
         go to 41
 21      continue
c-----------------
c   negative phase
c-----------------
         ix=0
         do 40 j=1,nbfi
            xx=zero
            do 30 k=1,nbfi
               ix=ix+1
               xx=xx+ci(k,im)*hess(ix)
 30         continue
            tv(j)=-xx
 40      continue
c
 41      continue
c------------------------------------------------------
c  finish off the transform for this row of the hessian
c------------------------------------------------------
         do 60 ii=1,i
            iim=mixi(ii)
            xx=zero
            do 50 k=1,nbfi
               xx=xx+tv(k)*ci(k,iim)
 50         continue
c---------------------------------------------
c   adjust the phase for the ii orbital mixings
c---------------------------------------------
            iix=iix+1
            if(iim.gt.iorb)go to 55
            htran(iix)=xx
            go to 60
 55         continue
            htran(iix)=-xx
 60      continue
c
 70   continue
c
c--------------------------------------------
c   put away the transformed hessian elements
c--------------------------------------------
c
c     if(iocore.eq.0) go to 100
c
      iix=0
      ii=0
      do 90 i=1,leni
         ix=(loci(i)-1)*nmix
         ii=ii+i
         htran(ii)=htran(ii)*.5
         do 80 j=1,i
            iix=iix+1
            hmo(ix+loci(j))=hmo(ix+loci(j))+htran(iix)
            hmo((loci(j)-1)*nmix+loci(i))=hmo((loci(j)-1)*nmix+loci(i))+
     $           htran(iix)
 80      continue
 90   continue
c
      return
c
c 100 continue
c
c     call swrit(nfhm,htran,intowp(iix))
c
cc    write(nfhm) ipt,ipt,leni,leni,iix,(htran(i),i=1,iix)
c
c     return
      end
