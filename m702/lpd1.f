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
c     this hardwires this code to handle l=5 (h functions) at the most
      parameter (highl=5)
      integer ipntr(0:highl,0:highl,0:highl)
      integer jpntr(0:highl,0:highl,0:highl)
      real*8 two,twoi,twoj
      integer putmwhr(21)
      integer nfoo1(3),nfoo2(3),nfoo3(3)
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
      call izero(ipntr,(highl+1)**3)
      call izero(jpntr,(highl+1)**3)
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
         nfoo1(1)=nx(i)
         nfoo1(2)=ny(i)
         nfoo1(3)=nz(i)
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
            nfoo2(1)=nx(j)
            nfoo2(2)=ny(j)
            nfoo2(3)=nz(j)
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
                  if (nderiv .eq. 2) then
c
c do the AA second derivs first.  First do the xx,yy,zz.
c                     
                     theint=1
                     do 6661 coord=1,3
c
c $%$)ing ugly g!*-d!*( hack to work around movd2e and save flops.  
c Ickypoo.
c
                        if (iatom .ne. katom) then
                           nfoo3(1)=nfoo1(1)
                           nfoo3(2)=nfoo1(2)
                           nfoo3(3)=nfoo1(3)
                           nfoo3(coord)=nfoo3(coord)+2
                           ip2=ipntr(nfoo3(1),nfoo3(2),nfoo3(3))
                           if (nfoo1(coord).gt.1) then
                              nfoo3(coord)=nfoo1(coord)-2
                              im2=ipntr(nfoo3(1),nfoo3(2),nfoo3(3))
                           endif
                           sp(iprim,j,i,putmwhr(theint))=
     $                          sp(iprim,j,i,putmwhr(theint))+
     $                          twoi*twoi*scr(iprim,j0,ip2)-
     $                          twoi*float(2*nfoo1(coord)+1)*
     $                          scr(iprim,j0,i0)
                           if (nfoo1(coord).gt.1) then
                              sp(iprim,j,i,putmwhr(theint))=
     $                             sp(iprim,j,i,putmwhr(theint))+
     $                             float(nfoo1(coord)*
     $                             (nfoo1(coord)-1))*scr(iprim,j0,im2)
                           endif
                        endif
                        theint=theint+1
 6661                continue 
c
c next xy,xz,yz
c
                     do 6662 coord1=1,2
                        do 6663 coord2=coord1+1,3
                           if (iatom .ne. katom) then
                              nfoo3(1)=nfoo1(1)
                              nfoo3(2)=nfoo1(2)
                              nfoo3(3)=nfoo1(3)
                              nfoo3(coord1)=nfoo3(coord1)+1
                              nfoo3(coord2)=nfoo3(coord2)+1
                              ipp=ipntr(nfoo3(1),nfoo3(2),nfoo3(3))
                              if (nfoo1(coord2).gt.0) then
                                 nfoo3(coord1)=nfoo1(coord1)+1
                                 nfoo3(coord2)=nfoo1(coord2)-1
                                 ipm=ipntr(nfoo3(1),nfoo3(2),nfoo3(3))
                              endif
                              if (nfoo1(coord1).gt.0) then
                                 nfoo3(coord1)=nfoo1(coord1)-1
                                 nfoo3(coord2)=nfoo1(coord2)+1
                                 imp=ipntr(nfoo3(1),nfoo3(2),nfoo3(3))
                                 if (nfoo1(coord2).gt.0) then
                                    nfoo3(coord1)=nfoo1(coord1)-1
                                    nfoo3(coord2)=nfoo1(coord2)-1
                                    imm=ipntr(nfoo3(1),nfoo3(2),
     $                                   nfoo3(3))
                                 endif
                              endif
                              sp(iprim,j,i,putmwhr(theint))=
     $                             sp(iprim,j,i,putmwhr(theint))+
     $                             twoi*twoi*scr(iprim,j0,ipp)
                              if (nfoo1(coord2).gt.0) then
                                 sp(iprim,j,i,putmwhr(theint))=
     $                                sp(iprim,j,i,putmwhr(theint))-
     $                                twoi*float(nfoo1(coord2))*
     $                                scr(iprim,j0,ipm)
                                 if (nfoo1(coord1).gt.0) then
                                    sp(iprim,j,i,putmwhr(theint))=
     $                                   sp(iprim,j,i,putmwhr(theint))+
     $                                   float(nfoo1(coord1)*
     $                                   nfoo1(coord2))*
     $                                   scr(iprim,j0,imm)
                                 endif
                              endif
                              if (nfoo1(coord1).gt.0) then
                                 sp(iprim,j,i,putmwhr(theint))=
     $                                sp(iprim,j,i,putmwhr(theint))-
     $                                twoi*float(nfoo1(coord1))*
     $                                scr(iprim,j0,imp)
                              endif
                           endif
                           theint=theint+1
 6663                   continue 
 6662                continue 
c
c now do the bb guys, x'x',y'y',z'z':
c                     
                     do 6664 coord=1,3
                        if (jatom .ne. katom) then
                           nfoo3(1)=nfoo2(1)
                           nfoo3(2)=nfoo2(2)
                           nfoo3(3)=nfoo2(3)
                           nfoo3(coord)=nfoo3(coord)+2
                           jp2=jpntr(nfoo3(1),nfoo3(2),nfoo3(3))
                           if (nfoo2(coord).gt.1) then
                              nfoo3(coord)=nfoo2(coord)-2
                              jm2=jpntr(nfoo3(1),nfoo3(2),nfoo3(3))
                           endif
                           sp(iprim,j,i,putmwhr(theint))=
     $                          sp(iprim,j,i,putmwhr(theint))+
     $                          twoj*twoj*scr(iprim,jp2,i0)-
     $                          twoj*float(2*nfoo2(coord)+1)*
     $                          scr(iprim,j0,i0)
                           if (nfoo2(coord).gt.1) then
                              sp(iprim,j,i,putmwhr(theint))=
     $                             sp(iprim,j,i,putmwhr(theint))+
     $                             float(nfoo2(coord)*
     $                             (nfoo2(coord)-1))*scr(iprim,jm2,i0)
                           endif
                        endif
                        theint=theint+1
 6664                continue 
c
c next x'y',x'z',y'z'
c
                     do 6665 coord1=1,2
                        do 6666 coord2=coord1+1,3
                           if (jatom .ne. katom) then
                              nfoo3(1)=nfoo2(1)
                              nfoo3(2)=nfoo2(2)
                              nfoo3(3)=nfoo2(3)
                              nfoo3(coord1)=nfoo3(coord1)+1
                              nfoo3(coord2)=nfoo3(coord2)+1
                              jpp=jpntr(nfoo3(1),nfoo3(2),nfoo3(3))
                              if (nfoo2(coord2).gt.0) then
                                 nfoo3(coord1)=nfoo2(coord1)+1
                                 nfoo3(coord2)=nfoo2(coord2)-1
                                 jpm=jpntr(nfoo3(1),nfoo3(2),nfoo3(3))
                              endif
                              if (nfoo2(coord1).gt.0) then
                                 nfoo3(coord1)=nfoo2(coord1)-1
                                 nfoo3(coord2)=nfoo2(coord2)+1
                                 jmp=jpntr(nfoo3(1),nfoo3(2),nfoo3(3))
                                 if (nfoo2(coord2).gt.0) then
                                    nfoo3(coord1)=nfoo2(coord1)-1
                                    nfoo3(coord2)=nfoo2(coord2)-1
                                    jmm=jpntr(nfoo3(1),nfoo3(2),
     $                                   nfoo3(3))
                                 endif
                              endif
                              sp(iprim,j,i,putmwhr(theint))=
     $                             sp(iprim,j,i,putmwhr(theint))+
     $                             twoj*twoj*scr(iprim,jpp,i0)
                              if (nfoo2(coord2).gt.0) then
                                 sp(iprim,j,i,putmwhr(theint))=
     $                                sp(iprim,j,i,putmwhr(theint))-
     $                                twoj*float(nfoo2(coord2))*
     $                                scr(iprim,jpm,i0)
                                 if (nfoo2(coord1).gt.0) then
                                    sp(iprim,j,i,putmwhr(theint))=
     $                                   sp(iprim,j,i,putmwhr(theint))+
     $                                   float(nfoo2(coord1)*
     $                                   nfoo2(coord2))*
     $                                   scr(iprim,jmm,i0)
                                 endif
                              endif
                              if (nfoo2(coord1).gt.0) then
                                 sp(iprim,j,i,putmwhr(theint))=
     $                                sp(iprim,j,i,putmwhr(theint))-
     $                                twoj*float(nfoo2(coord1))*
     $                                scr(iprim,jmp,i0)
                              endif
                           endif
                           theint=theint+1
 6666                   continue 
 6665                continue 
c
c now the other thing, the mixed guys.  Order done is
c  xx' xy' xz' yx' yy' yz' ...
c
                     if (iatom .ne. katom .and. jatom.ne.katom) then
                        do 6667 coord1=1,3
                           nfoo3(1)=nfoo1(1)
                           nfoo3(2)=nfoo1(2)
                           nfoo3(3)=nfoo1(3)
                           nfoo3(coord1)=nfoo1(coord1)+1
                           ip1=ipntr(nfoo3(1),nfoo3(2),nfoo3(3))
                           if (nfoo1(coord1).ne.0) then
                              nfoo3(coord1)=nfoo1(coord1)-1
                              im1=ipntr(nfoo3(1),nfoo3(2),nfoo3(3))
                           endif
                           do 6668 coord2=1,3
                              nfoo3(1)=nfoo2(1)
                              nfoo3(2)=nfoo2(2)
                              nfoo3(3)=nfoo2(3)
                              nfoo3(coord2)=nfoo2(coord2)+1
                              jp1=jpntr(nfoo3(1),nfoo3(2),nfoo3(3))
                              if (nfoo2(coord2).ne.0) then
                                 nfoo3(coord2)=nfoo2(coord2)-1
                                 jm1=jpntr(nfoo3(1),nfoo3(2),nfoo3(3))
                              endif
                              sp(iprim,j,i,putmwhr(theint))=
     $                             sp(iprim,j,i,putmwhr(theint))+
     $                             twoi*twoj*scr(iprim,jp1,ip1)
                              if (nfoo2(coord2).ne.0) then
                                 sp(iprim,j,i,putmwhr(theint))=
     $                                sp(iprim,j,i,putmwhr(theint))-
     $                                float(nfoo2(coord2))*twoi*
     $                                scr(iprim,jm1,ip1)
                                 if (nfoo1(coord1).ne.0) then
                                    sp(iprim,j,i,putmwhr(theint))=
     $                                   sp(iprim,j,i,putmwhr(theint))+
     $                                   float(nfoo1(coord1)*
     $                                   nfoo2(coord2))*
     $                                   scr(iprim,jm1,im1)
                                 endif
                              endif
                              if (nfoo1(coord1).ne.0) then
                                 sp(iprim,j,i,putmwhr(theint))=
     $                                sp(iprim,j,i,putmwhr(theint))-
     $                                twoj*float(nfoo1(coord1))*
     $                                scr(iprim,jp1,im1)
                              endif
                              theint=theint+1
 6668                      continue 
 6667                   continue 
                     endif
                  endif
   70          continue
   80       continue
   90    continue
  100 continue
c
c
      return
      end
