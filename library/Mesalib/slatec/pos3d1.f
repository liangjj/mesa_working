*deck pos3d1
      subroutine pos3d1 (lp, l, mp, m, n, a, b, c, ldimf, mdimf, f, xrt,
     +   yrt, t, d, wx, wy, c1, c2, bb)
c***begin prologue  pos3d1
c***subsidiary
c***purpose  subsidiary to pois3d
c***library   slatec
c***type      single precision (pos3d1-s)
c***author  (unknown)
c***see also  pois3d
c***routines called  cosqb, cosqf, cosqi, cost, costi, pimach, rfftb,
c                    rfftf, rffti, sinqb, sinqf, sinqi, sint, sinti,
c                    tridq
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891009  removed unreferenced variable.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900308  changed call to trid to call to tridq.  (wrb)
c   900402  added type section.  (wrb)
c***end prologue  pos3d1
      dimension       a(*)       ,b(*)       ,c(*)       ,
     1                f(ldimf,mdimf,*)       ,xrt(*)     ,yrt(*)     ,
     2                t(*)       ,d(*)       ,wx(*)      ,wy(*)      ,
     3                bb(*)
c***first executable statement  pos3d1
      pi = pimach(dum)
      lr = l
      mr = m
      nr = n
c
c     generate transform roots
c
      lrdel = ((lp-1)*(lp-3)*(lp-5))/3
      scalx = lr+lrdel
      dx = pi/(2.*scalx)
      go to (108,103,101,102,101),lp
  101 di = 0.5
      scalx = 2.*scalx
      go to 104
  102 di = 1.0
      go to 104
  103 di = 0.0
  104 do 105 i=1,lr
         xrt(i) = -4.*c1*(sin((i-di)*dx))**2
  105 continue
      scalx = 2.*scalx
      go to (112,106,110,107,111),lp
  106 call sinti (lr,wx)
      go to 112
  107 call costi (lr,wx)
      go to 112
  108 xrt(1) = 0.
      xrt(lr) = -4.*c1
      do 109 i=3,lr,2
         xrt(i-1) = -4.*c1*(sin((i-1)*dx))**2
         xrt(i) = xrt(i-1)
  109 continue
      call rffti (lr,wx)
      go to 112
  110 call sinqi (lr,wx)
      go to 112
  111 call cosqi (lr,wx)
  112 continue
      mrdel = ((mp-1)*(mp-3)*(mp-5))/3
      scaly = mr+mrdel
      dy = pi/(2.*scaly)
      go to (120,115,113,114,113),mp
  113 dj = 0.5
      scaly = 2.*scaly
      go to 116
  114 dj = 1.0
      go to 116
  115 dj = 0.0
  116 do 117 j=1,mr
         yrt(j) = -4.*c2*(sin((j-dj)*dy))**2
  117 continue
      scaly = 2.*scaly
      go to (124,118,122,119,123),mp
  118 call sinti (mr,wy)
      go to 124
  119 call costi (mr,wy)
      go to 124
  120 yrt(1) = 0.
      yrt(mr) = -4.*c2
      do 121 j=3,mr,2
         yrt(j-1) = -4.*c2*(sin((j-1)*dy))**2
         yrt(j) = yrt(j-1)
  121 continue
      call rffti (mr,wy)
      go to 124
  122 call sinqi (mr,wy)
      go to 124
  123 call cosqi (mr,wy)
  124 continue
      ifwrd = 1
  125 continue
c
c     transform x
c
      do 141 j=1,mr
         do 140 k=1,nr
            do 126 i=1,lr
               t(i) = f(i,j,k)
  126       continue
            go to (127,130,131,134,135),lp
  127       go to (128,129),ifwrd
  128       call rfftf (lr,t,wx)
            go to 138
  129       call rfftb (lr,t,wx)
            go to 138
  130       call sint (lr,t,wx)
            go to 138
  131       go to (132,133),ifwrd
  132       call sinqf (lr,t,wx)
            go to 138
  133       call sinqb (lr,t,wx)
            go to 138
  134       call cost (lr,t,wx)
            go to 138
  135       go to (136,137),ifwrd
  136       call cosqf (lr,t,wx)
            go to 138
  137       call cosqb (lr,t,wx)
  138       continue
            do 139 i=1,lr
               f(i,j,k) = t(i)
  139       continue
  140    continue
  141 continue
      go to (142,164),ifwrd
c
c     transform y
c
  142 continue
      do 158 i=1,lr
         do 157 k=1,nr
            do 143 j=1,mr
               t(j) = f(i,j,k)
  143       continue
            go to (144,147,148,151,152),mp
  144       go to (145,146),ifwrd
  145       call rfftf (mr,t,wy)
            go to 155
  146       call rfftb (mr,t,wy)
            go to 155
  147       call sint (mr,t,wy)
            go to 155
  148       go to (149,150),ifwrd
  149       call sinqf (mr,t,wy)
            go to 155
  150       call sinqb (mr,t,wy)
            go to 155
  151       call cost (mr,t,wy)
            go to 155
  152       go to (153,154),ifwrd
  153       call cosqf (mr,t,wy)
            go to 155
  154       call cosqb (mr,t,wy)
  155       continue
            do 156 j=1,mr
               f(i,j,k) = t(j)
  156       continue
  157    continue
  158 continue
      go to (159,125),ifwrd
  159 continue
c
c     solve tridiagonal systems in z
c
      do 163 i=1,lr
         do 162 j=1,mr
            do 160 k=1,nr
               bb(k) = b(k)+xrt(i)+yrt(j)
               t(k) = f(i,j,k)
  160       continue
            call tridq (nr,a,bb,c,t,d)
            do 161 k=1,nr
               f(i,j,k) = t(k)
  161       continue
  162    continue
  163 continue
      ifwrd = 2
      go to 142
  164 continue
      do 167 i=1,lr
         do 166 j=1,mr
            do 165 k=1,nr
               f(i,j,k) = f(i,j,k)/(scalx*scaly)
  165       continue
  166    continue
  167 continue
      return
      end
