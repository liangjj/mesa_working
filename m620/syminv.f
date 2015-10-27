*deck @(#)syminv.f	5.1  11/6/94
      subroutine syminv(am,nv,ng)
      implicit real*8 (a-h,o-z)
      dimension am(1)
c     diagonalize
      do 24 i=1,nv
      i0=i*nv-(i*(i-1))/2
      ii=i0+i
      li=(i-1)*nv+i-((i-1)*(i-2))/2
      am(li)=float(i)
      bigajj=0.0d0
      do 2 j=i,nv
      jj=j*nv+j-(j*(j-1))/2
      if (abs(am(jj))-bigajj)2,2,3
    3 bigajj=abs(am(jj))
      am(li)=float(j)
    2 continue
      call aintch(am,nv,i,i0,ii,li)
c     resume diagonalize
   20 if(abs(am(ii))-10.0d0**(-30)) 22,22,21
   22 ng=1
      go to 62
   21 am(ii)=1.0d0/am(ii)
      iq0=(i-1)*nv-((i-1)*(i-2))/2
      i1=i+1
      if (i1-nv)28,28,29
   28 do 23 k=i1,nv
      ik=i0+k
      ikq=iq0+k
   23 am(ikq)=-am(ii)*am(ik)
      do 24 j=i1,nv
      j0=j*nv-(j*(j-1))/2
      ij=i0+j
      do 24 k=j,nv
      jk=j0+k
      ikq=iq0+k
   24 am(jk)=am(jk)+am(ikq)*am(ij)
c     restore
   29 if (nv-1) 61,61,70
   70 i=nv-1
   55 i0=i*nv-(i*(i-1))/2
      ii=i0+i
      li=(i-1)*nv+i-((i-1)*(i-2))/2
      i1=i+1
      do 50 j=i1,nv
      ij=i0+j
   50 am(ij)=0.0d0
      iq0=(i-1)*nv-((i-1)*(i-2))/2
      do 52 j=i1,nv
      j0=j*nv-(j*(j-1))/2
      ij=i0+j
      ijq=iq0+j
      do 30 k=j,nv
      ik=i0+k
      jk=j0+k
      am(ik)=am(ik)+am(ijq)*am(jk)
      if (j-k)53,30,53
   53 ikq=iq0+k
      am(ij)=am(ij)+am(ikq)*am(jk)
   30 continue
   52 am(ii)=am(ii)+am(ijq)*am(ij)
c     interchange
      call aintch(am,nv,i,i0,ii,li)
      i=i-1
      if (i)61,61,55
   61 ng=0
   62 return
      end
