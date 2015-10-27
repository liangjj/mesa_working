      subroutine mcijki(dj,dkii,dkki,dkkj,ni,nj,nk,nl,nij,nki,nii,
     $     nkjt,nkj,den)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
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
c     implicit real*8(a-h,o-p,r-z),integer*2(q)
cc
cmp   extended dummy dj,dkii,dkki,dkkj,den
cc
      common /number/ zero,pt5,one,two,four,eight
      dimension dj(2),dkii(nii,nkjt),dkki(nki,nij),dkkj(nkj,2)
      dimension den(2)
c----------------------------c
c     mk > mi > mj
c----------------------------c
c
      ijki = 0
      ki=0
      ixx=0
c
      do 60 i = 1, ni
         do 50 k = 1, nk
            ki=ki+1
            kj=k
            ij=i
            do 40 j = 1, nj
               itot=i
               ix=ixx
               kii=k
               do 30 ii = 1, i
                  ix=ix+1
                  ijki = ijki + 1
                  dkii(ix,kj) = den(ijki)
                  dkkj(kj,itot)=den(ijki)
                  dkki(kii,ij) = den(ijki)
                  dj(ijki) = den(ijki)
                  kii=kii+nk
                  itot=itot+ni
 30            continue
               if(i.eq.ni) go to 25
               i1 = i + 1
               iix=ix+i
               do 20 ii = i1,ni
                  ijki = ijki + 1
                  dkii(iix,nkj+kj) = den(ijki)
                  dkkj(kj,itot)=den(ijki)
                  dkki(kii,ij) = den(ijki)
                  dj(ijki) = den(ijki)
                  iix=iix+ii
                  kii=kii+nk
                  itot=itot+ni
 20            continue
 25            continue
               ij=ij+ni
               kj=kj+nk
 40         continue
 50      continue
         itt=itt+ni
         ixx=ixx+i
 60   continue
c
      ix=0
      do 80 i=1,ni
         ix=ix+i
         do 70 kj=1,nkj
            dkii(ix,nkj+kj)=zero
 70      continue
 80   continue
c
c     reorder dkii  for yoshimine convention   ??!!?!"*&
c
      if(ni.eq.1) return
      ix=1
      do 95 i=2,ni
         i1=i-1
         do 94 j=1,i1
            ix=ix+1
            do 90 kj=1,nkj
               xx=dkii(ix,nkj+kj)
               dkii(ix,nkj+kj)=dkii(ix,kj)
               dkii(ix,kj)=xx
 90         continue
 94      continue
         ix=ix+1
 95   continue
c
c
      return
      end
