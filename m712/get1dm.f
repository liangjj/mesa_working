*deck @(#)get1dm.f	5.1  11/6/94
      subroutine get1dm(d,nnp,ndmat,dij,ni,nj,nfi,nfj,is,js)
c
c***begin prologue     get1dm
c***date written       870702   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           hf density matrices
c***author             saxe, paul (lanl)
c***source             @(#)get1dm.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       get1dm
c
      implicit integer (a-z)
c
      real*8 d(nnp,ndmat)
      real*8 dij(ni,nj,nfi,nfj,ndmat)
      integer ni,nj
      integer nfi,nfj
      integer is,js
c
c     ----- extract just the local ij portion of the full density -----
c
      do 50 dmat=1,ndmat
         do 40 jf=1,nfj
            jpos=js+jf-nfj
            do 30 if=1,nfi
               ipos=is+if-nfi
               do 20 j=1,nj
                  jj=jpos+j*nfj
                  do 10 i=1,ni
                     ii=ipos+i*nfi
                     ij=max(ii,jj)*(max(ii,jj)-1)/2+min(ii,jj)
                     dij(i,j,if,jf,dmat)=d(ij,dmat)
 10               continue
 20            continue
 30         continue
 40      continue
 50   continue
c
c
      return
      end
