*deck @(#)lpd1.f	5.1  11/6/94
      subroutine lpd1(nprim,ex,igbegn,igend,jgbegn,jgend,maxp2,
     $     sp,mini,maxi,minj,maxj,
     $     scr,dmini,dmaxi,dminj,dmaxj,nx,ny,nz,nderiv,
     $     iatom,jatom,katom)
      implicit integer(a-z)
      real*8 ex(nprim)
c last dim of sp could be 7 or 28, depending on nderiv=1 or 2
      real*8 sp(maxp2,minj:maxj,mini:maxi,*)
      real*8 scr(maxp2,dminj:dmaxj,dmini:dmaxi)
      dimension nx(*),ny(*),nz(*)
c     this hard-wires this code into handling l=5(h functions) at the most.
c     i.e., derivatives of g functions
      parameter (highl=5)
      integer ipntr(0:highl,0:highl,0:highl)
      integer jpntr(0:highl,0:highl,0:highl)
      real*8 two,twoi,twoj
      integer putmwhr(21)
      common/io/inp,iout
c
c I'll be generating the second derivatives in the following order:
c xx yy zz xy xz yz x'x' y'y' z'z' x'y' x'z' y'z' xx' xy' xz' yx' yy' yz'
c   zx' zy' zz'
c but they need to be stuck into sp in a different place.  Here's the 
c mapping
c
      data putmwhr/8,10,13,9,11,12,17,22,28,21,26,27,14,18,23,15,19,
     $     24,16,20,25/
      parameter (two=2.0d+00)
c
c     prepare pointer table.
      call izero(ipntr,125)
      call izero(jpntr,125)
      do 10 i=dmini,dmaxi
         ipntr(nx(i),ny(i),nz(i))=i
   10 continue
      do 20 j=dminj,dmaxj
         jpntr(nx(j),ny(j),nz(j))=j
   20 continue
      do 100 i=mini,maxi
         i0=ipntr(nx(i),ny(i),nz(i))
         ixp=ipntr(nx(i)+1,ny(i),nz(i))
         iyp=ipntr(nx(i),ny(i)+1,nz(i))
         izp=ipntr(nx(i),ny(i),nz(i)+1)
         if(nx(i).ne.0) then
            ixm=ipntr(nx(i)-1,ny(i),nz(i))
         endif
         if(ny(i).ne.0) then
            iym=ipntr(nx(i),ny(i)-1,nz(i))
         endif
         if(nz(i).ne.0) then
            izm=ipntr(nx(i),ny(i),nz(i)-1)
         endif
         do 90 j=minj,maxj
            j0=jpntr(nx(j),ny(j),nz(j))
            jxp=jpntr(nx(j)+1,ny(j),nz(j))
            jyp=jpntr(nx(j),ny(j)+1,nz(j))
            jzp=jpntr(nx(j),ny(j),nz(j)+1)
            if(nx(j).ne.0) then
               jxm=jpntr(nx(j)-1,ny(j),nz(j))
            endif
            if(ny(j).ne.0) then
               jym=jpntr(nx(j),ny(j)-1,nz(j))
            endif
            if(nz(j).ne.0) then
               jzm=jpntr(nx(j),ny(j),nz(j)-1)
            endif
c
c
            iprim=0
            do 80 jgauss=jgbegn,jgend
               twoj=two*ex(jgauss)
               do 70 igauss=igbegn,igend
                  twoi=two*ex(igauss)
                  iprim=iprim+1
                  sp(iprim,j,i,1)=sp(iprim,j,i,1)+scr(iprim,j0,i0)
c                 center i. 
                  sp(iprim,j,i,2)=sp(iprim,j,i,2)+twoi*scr(iprim,j0,ixp)
                  if(nx(i).ne.0) then
                     sp(iprim,j,i,2)=sp(iprim,j,i,2)
     $                               -float(nx(i))*scr(iprim,j0,ixm)
                  endif
                  sp(iprim,j,i,3)=sp(iprim,j,i,3)+twoi*scr(iprim,j0,iyp)
                  if(ny(i).ne.0) then
                     sp(iprim,j,i,3)=sp(iprim,j,i,3)
     $                               -float(ny(i))*scr(iprim,j0,iym)
                  endif
                  sp(iprim,j,i,4)=sp(iprim,j,i,4)+twoi*scr(iprim,j0,izp)
                  if(nz(i).ne.0) then
                     sp(iprim,j,i,4)=sp(iprim,j,i,4)
     $                               -float(nz(i))*scr(iprim,j0,izm)
                  endif
c                 center j.
                  sp(iprim,j,i,5)=sp(iprim,j,i,5)+twoj*scr(iprim,jxp,i0)
                  if(nx(j).ne.0) then
                     sp(iprim,j,i,5)=sp(iprim,j,i,5)
     $                               -float(nx(j))*scr(iprim,jxm,i0)
                  endif
                  sp(iprim,j,i,6)=sp(iprim,j,i,6)+twoj*scr(iprim,jyp,i0)
                  if(ny(j).ne.0) then
                     sp(iprim,j,i,6)=sp(iprim,j,i,6)
     $                               -float(ny(j))*scr(iprim,jym,i0)
                  endif
                  sp(iprim,j,i,7)=sp(iprim,j,i,7)+twoj*scr(iprim,jzp,i0)
                  if(nz(j).ne.0) then
                     sp(iprim,j,i,7)=sp(iprim,j,i,7)
     $                               -float(nz(j))*scr(iprim,jzm,i0)
                  endif
 70            continue
 80         continue
 90      continue
 100  continue
c
c
      return
      end
