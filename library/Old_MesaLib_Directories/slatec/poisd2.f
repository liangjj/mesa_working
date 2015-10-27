*deck poisd2
      subroutine poisd2 (mr, nr, istag, ba, bb, bc, q, idimq, b, w, d,
     +   tcos, p)
c***begin prologue  poisd2
c***subsidiary
c***purpose  subsidiary to genbun
c***library   slatec
c***type      single precision (poisd2-s, cmposd-c)
c***author  (unknown)
c***description
c
c     subroutine to solve poisson's equation for dirichlet boundary
c     conditions.
c
c     istag = 1 if the last diagonal block is the matrix a.
c     istag = 2 if the last diagonal block is the matrix a+i.
c
c***see also  genbun
c***routines called  cosgen, s1merg, trix
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   920130  modified to use merge routine s1merg rather than deleted
c           routine merge.  (wrb)
c***end prologue  poisd2
c
      dimension       q(idimq,*) ,ba(*)      ,bb(*)      ,bc(*)      ,
     1                tcos(*)    ,b(*)       ,d(*)       ,w(*)       ,
     2                p(*)
c***first executable statement  poisd2
      m = mr
      n = nr
      jsh = 0
      fi = 1./istag
      ip = -m
      ipstor = 0
      go to (101,102),istag
  101 kr = 0
      irreg = 1
      if (n .gt. 1) go to 106
      tcos(1) = 0.
      go to 103
  102 kr = 1
      jstsav = 1
      irreg = 2
      if (n .gt. 1) go to 106
      tcos(1) = -1.
  103 do 104 i=1,m
         b(i) = q(i,1)
  104 continue
      call trix (1,0,m,ba,bb,bc,b,tcos,d,w)
      do 105 i=1,m
         q(i,1) = b(i)
  105 continue
      go to 183
  106 lr = 0
      do 107 i=1,m
         p(i) = 0.
  107 continue
      nun = n
      jst = 1
      jsp = n
c
c     irreg = 1 when no irregularities have occurred, otherwise it is 2.
c
  108 l = 2*jst
      nodd = 2-2*((nun+1)/2)+nun
c
c     nodd = 1 when nun is odd, otherwise it is 2.
c
      go to (110,109),nodd
  109 jsp = jsp-l
      go to 111
  110 jsp = jsp-jst
      if (irreg .ne. 1) jsp = jsp-l
  111 continue
c
c     regular reduction
c
      call cosgen (jst,1,0.5,0.0,tcos)
      if (l .gt. jsp) go to 118
      do 117 j=l,jsp,l
         jm1 = j-jsh
         jp1 = j+jsh
         jm2 = j-jst
         jp2 = j+jst
         jm3 = jm2-jsh
         jp3 = jp2+jsh
         if (jst .ne. 1) go to 113
         do 112 i=1,m
            b(i) = 2.*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  112    continue
         go to 115
  113    do 114 i=1,m
            t = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = t+q(i,j)-q(i,jm3)-q(i,jp3)
            q(i,j) = t
  114    continue
  115    continue
         call trix (jst,0,m,ba,bb,bc,b,tcos,d,w)
         do 116 i=1,m
            q(i,j) = q(i,j)+b(i)
  116    continue
  117 continue
c
c     reduction for last unknown
c
  118 go to (119,136),nodd
  119 go to (152,120),irreg
c
c     odd number of unknowns
c
  120 jsp = jsp+l
      j = jsp
      jm1 = j-jsh
      jp1 = j+jsh
      jm2 = j-jst
      jp2 = j+jst
      jm3 = jm2-jsh
      go to (123,121),istag
  121 continue
      if (jst .ne. 1) go to 123
      do 122 i=1,m
         b(i) = q(i,j)
         q(i,j) = 0.
  122 continue
      go to 130
  123 go to (124,126),noddpr
  124 do 125 i=1,m
         ip1 = ip+i
         b(i) = .5*(q(i,jm2)-q(i,jm1)-q(i,jm3))+p(ip1)+q(i,j)
  125 continue
      go to 128
  126 do 127 i=1,m
         b(i) = .5*(q(i,jm2)-q(i,jm1)-q(i,jm3))+q(i,jp2)-q(i,jp1)+q(i,j)
  127 continue
  128 do 129 i=1,m
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
  129 continue
  130 call trix (jst,0,m,ba,bb,bc,b,tcos,d,w)
      ip = ip+m
      ipstor = max(ipstor,ip+m)
      do 131 i=1,m
         ip1 = ip+i
         p(ip1) = q(i,j)+b(i)
         b(i) = q(i,jp2)+p(ip1)
  131 continue
      if (lr .ne. 0) go to 133
      do 132 i=1,jst
         krpi = kr+i
         tcos(krpi) = tcos(i)
  132 continue
      go to 134
  133 continue
      call cosgen (lr,jstsav,0.,fi,tcos(jst+1))
      call s1merg (tcos,0,jst,jst,lr,kr)
  134 continue
      call cosgen (kr,jstsav,0.0,fi,tcos)
      call trix (kr,kr,m,ba,bb,bc,b,tcos,d,w)
      do 135 i=1,m
         ip1 = ip+i
         q(i,j) = q(i,jm2)+b(i)+p(ip1)
  135 continue
      lr = kr
      kr = kr+l
      go to 152
c
c     even number of unknowns
c
  136 jsp = jsp+l
      j = jsp
      jm1 = j-jsh
      jp1 = j+jsh
      jm2 = j-jst
      jp2 = j+jst
      jm3 = jm2-jsh
      go to (137,138),irreg
  137 continue
      jstsav = jst
      ideg = jst
      kr = l
      go to 139
  138 call cosgen (kr,jstsav,0.0,fi,tcos)
      call cosgen (lr,jstsav,0.0,fi,tcos(kr+1))
      ideg = kr
      kr = kr+jst
  139 if (jst .ne. 1) go to 141
      irreg = 2
      do 140 i=1,m
         b(i) = q(i,j)
         q(i,j) = q(i,jm2)
  140 continue
      go to 150
  141 do 142 i=1,m
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  142 continue
      go to (143,145),irreg
  143 do 144 i=1,m
         q(i,j) = q(i,jm2)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
  144 continue
      irreg = 2
      go to 150
  145 continue
      go to (146,148),noddpr
  146 do 147 i=1,m
         ip1 = ip+i
         q(i,j) = q(i,jm2)+p(ip1)
  147 continue
      ip = ip-m
      go to 150
  148 do 149 i=1,m
         q(i,j) = q(i,jm2)+q(i,j)-q(i,jm1)
  149 continue
  150 call trix (ideg,lr,m,ba,bb,bc,b,tcos,d,w)
      do 151 i=1,m
         q(i,j) = q(i,j)+b(i)
  151 continue
  152 nun = nun/2
      noddpr = nodd
      jsh = jst
      jst = 2*jst
      if (nun .ge. 2) go to 108
c
c     start solution.
c
      j = jsp
      do 153 i=1,m
         b(i) = q(i,j)
  153 continue
      go to (154,155),irreg
  154 continue
      call cosgen (jst,1,0.5,0.0,tcos)
      ideg = jst
      go to 156
  155 kr = lr+jst
      call cosgen (kr,jstsav,0.0,fi,tcos)
      call cosgen (lr,jstsav,0.0,fi,tcos(kr+1))
      ideg = kr
  156 continue
      call trix (ideg,lr,m,ba,bb,bc,b,tcos,d,w)
      jm1 = j-jsh
      jp1 = j+jsh
      go to (157,159),irreg
  157 do 158 i=1,m
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
  158 continue
      go to 164
  159 go to (160,162),noddpr
  160 do 161 i=1,m
         ip1 = ip+i
         q(i,j) = p(ip1)+b(i)
  161 continue
      ip = ip-m
      go to 164
  162 do 163 i=1,m
         q(i,j) = q(i,j)-q(i,jm1)+b(i)
  163 continue
  164 continue
c
c     start back substitution.
c
      jst = jst/2
      jsh = jst/2
      nun = 2*nun
      if (nun .gt. n) go to 183
      do 182 j=jst,n,l
         jm1 = j-jsh
         jp1 = j+jsh
         jm2 = j-jst
         jp2 = j+jst
         if (j .gt. jst) go to 166
         do 165 i=1,m
            b(i) = q(i,j)+q(i,jp2)
  165    continue
         go to 170
  166    if (jp2 .le. n) go to 168
         do 167 i=1,m
            b(i) = q(i,j)+q(i,jm2)
  167    continue
         if (jst .lt. jstsav) irreg = 1
         go to (170,171),irreg
  168    do 169 i=1,m
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  169    continue
  170    continue
         call cosgen (jst,1,0.5,0.0,tcos)
         ideg = jst
         jdeg = 0
         go to 172
  171    if (j+l .gt. n) lr = lr-jst
         kr = jst+lr
         call cosgen (kr,jstsav,0.0,fi,tcos)
         call cosgen (lr,jstsav,0.0,fi,tcos(kr+1))
         ideg = kr
         jdeg = lr
  172    continue
         call trix (ideg,jdeg,m,ba,bb,bc,b,tcos,d,w)
         if (jst .gt. 1) go to 174
         do 173 i=1,m
            q(i,j) = b(i)
  173    continue
         go to 182
  174    if (jp2 .gt. n) go to 177
  175    do 176 i=1,m
            q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
  176    continue
         go to 182
  177    go to (175,178),irreg
  178    if (j+jsh .gt. n) go to 180
         do 179 i=1,m
            ip1 = ip+i
            q(i,j) = b(i)+p(ip1)
  179    continue
         ip = ip-m
         go to 182
  180    do 181 i=1,m
            q(i,j) = b(i)+q(i,j)-q(i,jm1)
  181    continue
  182 continue
      l = l/2
      go to 164
  183 continue
c
c     return storage requirements for p vectors.
c
      w(1) = ipstor
      return
      end
