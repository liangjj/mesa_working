      subroutine mcijkl(dj,dkki,dkkj,ni,nj,nk,nl,nij,nkl,nki,njl,nil,
     $     nkj,den,mi,mj,mk,ml)
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
c     implicit real*8(a-h,o-p,r-z),             integer*2(q)
cc
cmp   extended dummy dj,dkki,dkkj,den
cc
      dimension dj(2),dkki(nki,njl),dkkj(nkj,nil)
      dimension den(2)
c------------------------------------------------------c
c     mi  mj  mk  ml   (mk>mi)   (mk>ml)   (mi>mj)
c------------------------------------------------------c
      ijkl = 0
c
      if(ml.lt.mi.or.ml.lt.mj) go to 100
c
c---------------------c
c    ml > mi,mj
c---------------------c
c
      do 60 l = 1, nl
         do 50 k = 1, nk
            jxj=0
            lxj=l
            do 40 j = 1, nj
               ixi=0
               lxi=l
               do 30 i = 1, ni
                  ijkl = ijkl + 1
                  dkki(ixi+k,lxj) = den(ijkl)
                  dkkj(jxj+k,lxi) = den(ijkl)
                  lxi=lxi+nl
                  ixi=ixi+nk
                  dj(ijkl) = den(ijkl)
 30            continue
               jxj=jxj+nk
               lxj=lxj+nl
 40         continue
 50      continue
 60   continue
      go to 400
c
 100  continue
c
      if(ml.gt.mj) go to 200
c
c---------------------c
c    mi , mj > ml
c---------------------c
c
      ixl=0
      jxl=0
      do 160 l = 1, nl
         do 150 k = 1, nk
            jxj=0
            do 140 j = 1, nj
               ixi=0
               jl=jxl+j
               do 130 i = 1, ni
                  ijkl = ijkl + 1
                  dkki(ixi+k,jl) = den(ijkl)
                  dkkj(jxj+k,ixl+i) = den(ijkl)
                  ixi=ixi+nk
                  dj(ijkl) = den(ijkl)
 130           continue
               jxj=jxj+nk
 140        continue
 150     continue
         ixl=ixl+ni
         jxl=jxl+nj
 160  continue
      go to 400
c
 200  continue
c
c---------------------c
c    mi > ml > mj
c---------------------c
c
      lxl=0
      do 260 l = 1, nl
         do 250 k = 1, nk
            jxj=0
            lxj=l
            do 240 j = 1, nj
               ixi=0
               do 230 i = 1, ni
                  ijkl = ijkl + 1
                  dkki(ixi+k,lxj) = den(ijkl)
                  dkkj(jxj+k,lxl+i) = den(ijkl)
                  ixi=ixi+nk
                  dj(ijkl) = den(ijkl)
 230           continue
               jxj=jxj+nk
               lxj=lxj+nl
 240        continue
 250     continue
         lxl=lxl+ni
 260  continue
c
 400  continue
c
      return
      end
