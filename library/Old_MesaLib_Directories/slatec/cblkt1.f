*deck cblkt1
      subroutine cblkt1 (n, an, bn, cn, m, am, bm, cm, idimy, y, b, w1,
     +   w2, w3, wd, ww, wu, prdct, cprdct)
c***begin prologue  cblkt1
c***subsidiary
c***purpose  subsidiary to cblktr
c***library   slatec
c***type      complex (blktr1-s, cblkt1-c)
c***author  (unknown)
c***description
c
c cblkt1 solves the linear system of routine cblktr.
c
c b  contains the roots of all the b polynomials.
c w1,w2,w3,wd,ww,wu  are all working arrays.
c prdct is either procp or proc depending on whether the boundary
c conditions in the m direction are periodic or not.
c cprdct is either cprocp or cproc which are called if some of the zeros
c of the b polynomials are complex.
c
c***see also  cblktr
c***routines called  inxca, inxcb, inxcc
c***common blocks    ccblk
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  cblkt1
c
      dimension       an(*)      ,bn(*)      ,cn(*)      ,am(*)      ,
     1                bm(*)      ,cm(*)      ,b(*)       ,w1(*)      ,
     2                w2(*)      ,w3(*)      ,wd(*)      ,ww(*)      ,
     3                wu(*)      ,y(idimy,*)
      common /ccblk/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
      complex         am         ,bm         ,cm         ,y          ,
     1                w1         ,w2         ,w3         ,wd         ,
     2                ww         ,wu
c***first executable statement  cblkt1
      kdo = k-1
      do 109 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i3 = i2+i1
         i4 = i2+i2
         irm1 = ir-1
         call inxcb (i2,ir,im2,nm2)
         call inxcb (i1,irm1,im3,nm3)
         call inxcb (i3,irm1,im1,nm1)
         call prdct (nm2,b(im2),nm3,b(im3),nm1,b(im1),0,dum,y(1,i2),w3,
     1               m,am,bm,cm,wd,ww,wu)
         if = 2**k
         do 108 i=i4,if,i4
            if (i-nm) 101,101,108
  101       ipi1 = i+i1
            ipi2 = i+i2
            ipi3 = i+i3
            call inxcc (i,ir,idxc,nc)
            if (i-if) 102,108,108
  102       call inxca (i,ir,idxa,na)
            call inxcb (i-i1,irm1,im1,nm1)
            call inxcb (ipi2,ir,ip2,np2)
            call inxcb (ipi1,irm1,ip1,np1)
            call inxcb (ipi3,irm1,ip3,np3)
            call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w3,w1,m,am,
     1                  bm,cm,wd,ww,wu)
            if (ipi2-nm) 105,105,103
  103       do 104 j=1,m
               w3(j) = (0.,0.)
               w2(j) = (0.,0.)
  104       continue
            go to 106
  105       call prdct (np2,b(ip2),np1,b(ip1),np3,b(ip3),0,dum,
     1                  y(1,ipi2),w3,m,am,bm,cm,wd,ww,wu)
            call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w3,w2,m,am,
     1                  bm,cm,wd,ww,wu)
  106       do 107 j=1,m
               y(j,i) = w1(j)+w2(j)+y(j,i)
  107       continue
  108    continue
  109 continue
      if (npp) 132,110,132
c
c     the periodic case is treated using the capacitance matrix method
c
  110 if = 2**k
      i = if/2
      i1 = i/2
      call inxcb (i-i1,k-2,im1,nm1)
      call inxcb (i+i1,k-2,ip1,np1)
      call inxcb (i,k-1,iz,nz)
      call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,y(1,i),w1,m,am,
     1            bm,cm,wd,ww,wu)
      izr = i
      do 111 j=1,m
         w2(j) = w1(j)
  111 continue
      do 113 ll=2,k
         l = k-ll+1
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i = i2
         call inxcc (i,ir,idxc,nc)
         call inxcb (i,ir,iz,nz)
         call inxcb (i-i1,ir-1,im1,nm1)
         call inxcb (i+i1,ir-1,ip1,np1)
         call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w1,w1,m,am,bm,
     1               cm,wd,ww,wu)
         do 112 j=1,m
            w1(j) = y(j,i)+w1(j)
  112    continue
         call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w1,w1,m,am,
     1               bm,cm,wd,ww,wu)
  113 continue
      do 118 ll=2,k
         l = k-ll+1
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i4 = i2+i2
         ifd = if-i2
         do 117 i=i2,ifd,i4
            if (i-i2-izr) 117,114,117
  114       if (i-nm) 115,115,118
  115       call inxca (i,ir,idxa,na)
            call inxcb (i,ir,iz,nz)
            call inxcb (i-i1,ir-1,im1,nm1)
            call inxcb (i+i1,ir-1,ip1,np1)
            call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w2,w2,m,am,
     1                  bm,cm,wd,ww,wu)
            do 116 j=1,m
               w2(j) = y(j,i)+w2(j)
  116       continue
            call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w2,w2,m,
     1                  am,bm,cm,wd,ww,wu)
            izr = i
            if (i-nm) 117,119,117
  117    continue
  118 continue
  119 do 120 j=1,m
         y(j,nm+1) = y(j,nm+1)-cn(nm+1)*w1(j)-an(nm+1)*w2(j)
  120 continue
      call inxcb (if/2,k-1,im1,nm1)
      call inxcb (if,k-1,ip,np)
      if (ncmplx) 121,122,121
  121 call cprdct (nm+1,b(ip),nm1,b(im1),0,dum,0,dum,y(1,nm+1),
     1             y(1,nm+1),m,am,bm,cm,w1,w3,ww)
      go to 123
  122 call prdct (nm+1,b(ip),nm1,b(im1),0,dum,0,dum,y(1,nm+1),
     1            y(1,nm+1),m,am,bm,cm,wd,ww,wu)
  123 do 124 j=1,m
         w1(j) = an(1)*y(j,nm+1)
         w2(j) = cn(nm)*y(j,nm+1)
         y(j,1) = y(j,1)-w1(j)
         y(j,nm) = y(j,nm)-w2(j)
  124 continue
      do 126 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i4 = i2+i2
         i1 = i2/2
         i = i4
         call inxca (i,ir,idxa,na)
         call inxcb (i-i2,ir,im2,nm2)
         call inxcb (i-i2-i1,ir-1,im3,nm3)
         call inxcb (i-i1,ir-1,im1,nm1)
         call prdct (nm2,b(im2),nm3,b(im3),nm1,b(im1),0,dum,w1,w1,m,am,
     1               bm,cm,wd,ww,wu)
         call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w1,w1,m,am,bm,
     1               cm,wd,ww,wu)
         do 125 j=1,m
            y(j,i) = y(j,i)-w1(j)
  125    continue
  126 continue
c
      izr = nm
      do 131 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i3 = i2+i1
         i4 = i2+i2
         irm1 = ir-1
         do 130 i=i4,if,i4
            ipi1 = i+i1
            ipi2 = i+i2
            ipi3 = i+i3
            if (ipi2-izr) 127,128,127
  127       if (i-izr) 130,131,130
  128       call inxcc (i,ir,idxc,nc)
            call inxcb (ipi2,ir,ip2,np2)
            call inxcb (ipi1,irm1,ip1,np1)
            call inxcb (ipi3,irm1,ip3,np3)
            call prdct (np2,b(ip2),np1,b(ip1),np3,b(ip3),0,dum,w2,w2,m,
     1                  am,bm,cm,wd,ww,wu)
            call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w2,w2,m,am,
     1                  bm,cm,wd,ww,wu)
            do 129 j=1,m
               y(j,i) = y(j,i)-w2(j)
  129       continue
            izr = i
            go to 131
  130    continue
  131 continue
c
c begin back substitution phase
c
  132 do 144 ll=1,k
         l = k-ll+1
         ir = l-1
         irm1 = ir-1
         i2 = 2**ir
         i1 = i2/2
         i4 = i2+i2
         ifd = if-i2
         do 143 i=i2,ifd,i4
            if (i-nm) 133,133,143
  133       imi1 = i-i1
            imi2 = i-i2
            ipi1 = i+i1
            ipi2 = i+i2
            call inxca (i,ir,idxa,na)
            call inxcc (i,ir,idxc,nc)
            call inxcb (i,ir,iz,nz)
            call inxcb (imi1,irm1,im1,nm1)
            call inxcb (ipi1,irm1,ip1,np1)
            if (i-i2) 134,134,136
  134       do 135 j=1,m
               w1(j) = (0.,0.)
  135       continue
            go to 137
  136       call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),y(1,imi2),
     1                  w1,m,am,bm,cm,wd,ww,wu)
  137       if (ipi2-nm) 140,140,138
  138       do 139 j=1,m
               w2(j) = (0.,0.)
  139       continue
            go to 141
  140       call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),y(1,ipi2),
     1                  w2,m,am,bm,cm,wd,ww,wu)
  141       do 142 j=1,m
               w1(j) = y(j,i)+w1(j)+w2(j)
  142       continue
            call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w1,y(1,i),
     1                  m,am,bm,cm,wd,ww,wu)
  143    continue
  144 continue
      return
      end
