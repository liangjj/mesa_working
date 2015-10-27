      subroutine elim(nbf,lnew,mnew,nnew,eta,nfirst,nlast,
     x allrho,nstate,maxp,maxnbf,maxprim)
c
c  eliminate basis functions for which there are no nonzero
c  density matrix elements
c
      implicit real*8 (a-h,o-z)
      dimension lnew(maxprim), mnew(maxprim), nnew(maxprim)
      dimension eta(maxprim,5)
      dimension nfirst(maxnbf), nlast(maxnbf)
      dimension allrho((maxnbf*(maxnbf+1))/2,maxp)
      nbfold=nbf
      ia = 1
c
c test and eliminate function ia if corresponding density
c matrix elements are all = 0.
c
  1   test = 0.0
      do 10 ib=1,nbf
      if(ib.gt.ia) then
      indx = ib*(ib-1)/2 + ia
      else
      indx = ia*(ia-1)/2 + ib
      endif
      do 5 is=1,nstate
  5   test=test+ abs(allrho(indx,is))
  10  continue
c
      if(test.lt.1.e-8) then
c
c eliminate function ia from density matrices
c
      nnp1 = nbf*(nbf+1)/2
      do 20 ib=1,nbf
      if(ib.gt.ia) then
      indx = ib*(ib-1)/2 + ia
      else
      indx = ia*(ia-1)/2 + ib
      endif
      ntop = nnp1-ib
      indx = indx -ib +1
      do 15 jndx=indx,ntop
      do 15 js=1,nstate
  15  allrho(jndx,js) = allrho(jndx+1,js)
  20  continue
c
c eliminate function ia from basis function description vectors
c
      ilast = nlast(nbf)
      ilow = nfirst(ia)
      ihi = nlast(ia)
      nbfm1 = nbf-1
      ndel=nfirst(ia+1)-nfirst(ia)
      do 25 i=ia,nbfm1
      nfirst(i) = nfirst(i+1) - ndel
  25  nlast(i) = nlast(i+1) - ndel
      ihp1 = ihi+1
      do 30 i=ihp1,ilast
      isub = i - (ihi-ilow) - 1
      lnew(isub) = lnew(i)
      mnew(isub) = mnew(i)
      nnew(isub) = nnew(i)
      do 27 j=1,5
  27  eta(isub,j) = eta(i,j)
  30  continue
c
c  now this function has been eliminated both from the basis
c  function vectors and from the density matrices,  so reset
c  nbf and leave the pointer (ia) alone because it now points
c  to a new function.
c
      nbf=nbf-1
      else
c
c  we did not eliminate ia so move the pointer
c
      ia = ia + 1
      endif
c
c  if ia points to nbf+1, we are finished.  otherwise go back
c  and test function ia
c
c
      if(ia.le.nbf) go to 1
c
      write(66,100) nbfold,nbf
  100 format(///,' numbers of contracted functions before and after ',/,
     x ' elimination of noncontributors are',2i4)
      write(66,200)
200   format(///,' basis functions indexed for potential generation')
      do 500 i=1,nbf
      istrt=nfirst(i)
      il=nlast(i)
      write(66,101) i, istrt,il
101   format(//,' basis fcn ',i3,' nfirst and nlast are ',2i4)
      do 400 ii=istrt,il
      write(66,102) lnew(ii),mnew(ii),nnew(ii),(eta(ii,jj),jj=1,5)
102   format(1x,3i3,5f10.5)
400   continue
500   continue
      return
      end
