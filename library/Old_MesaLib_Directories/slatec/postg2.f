*deck postg2
      subroutine postg2 (nperod, n, m, a, bb, c, idimq, q, b, b2, b3, w,
     +   w2, w3, d, tcos, p)
c***begin prologue  postg2
c***subsidiary
c***purpose  subsidiary to poistg
c***library   slatec
c***type      single precision (postg2-s)
c***author  (unknown)
c***description
c
c     subroutine to solve poisson's equation on a staggered grid.
c
c***see also  poistg
c***routines called  cosgen, s1merg, tri3, trix
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   920130  modified to use merge routine s1merg rather than deleted
c           routine merge.  (wrb)
c***end prologue  postg2
c
      dimension       a(*)       ,bb(*)      ,c(*)       ,q(idimq,*) ,
     1                b(*)       ,b2(*)      ,b3(*)      ,w(*)       ,
     2                w2(*)      ,w3(*)      ,d(*)       ,tcos(*)    ,
     3                k(4)       ,p(*)
      equivalence     (k(1),k1)  ,(k(2),k2)  ,(k(3),k3)  ,(k(4),k4)
c***first executable statement  postg2
      np = nperod
      fnum = 0.5*(np/3)
      fnum2 = 0.5*(np/2)
      mr = m
      ip = -mr
      ipstor = 0
      i2r = 1
      jr = 2
      nr = n
      nlast = n
      kr = 1
      lr = 0
      if (nr .le. 3) go to 142
  101 continue
      jr = 2*i2r
      nrod = 1
      if ((nr/2)*2 .eq. nr) nrod = 0
      jstart = 1
      jstop = nlast-jr
      if (nrod .eq. 0) jstop = jstop-i2r
      i2rby2 = i2r/2
      if (jstop .ge. jstart) go to 102
      j = jr
      go to 115
  102 continue
c
c     regular reduction.
c
      ijump = 1
      do 114 j=jstart,jstop,jr
         jp1 = j+i2rby2
         jp2 = j+i2r
         jp3 = jp2+i2rby2
         jm1 = j-i2rby2
         jm2 = j-i2r
         jm3 = jm2-i2rby2
         if (j .ne. 1) go to 106
         call cosgen (i2r,1,fnum,0.5,tcos)
         if (i2r .ne. 1) go to 104
         do 103 i=1,mr
            b(i) = q(i,1)
            q(i,1) = q(i,2)
  103    continue
         go to 112
  104    do 105 i=1,mr
            b(i) = q(i,1)+0.5*(q(i,jp2)-q(i,jp1)-q(i,jp3))
            q(i,1) = q(i,jp2)+q(i,1)-q(i,jp1)
  105    continue
         go to 112
  106    continue
         go to (107,108),ijump
  107    continue
         ijump = 2
         call cosgen (i2r,1,0.5,0.0,tcos)
  108    continue
         if (i2r .ne. 1) go to 110
         do 109 i=1,mr
            b(i) = 2.*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  109    continue
         go to 112
  110    do 111 i=1,mr
            fi = q(i,j)
            q(i,j) = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = fi+q(i,j)-q(i,jm3)-q(i,jp3)
  111    continue
  112    continue
         call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
         do 113 i=1,mr
            q(i,j) = q(i,j)+b(i)
  113    continue
c
c     end of reduction for regular unknowns.
c
  114 continue
c
c     begin special reduction for last unknown.
c
      j = jstop+jr
  115 nlast = j
      jm1 = j-i2rby2
      jm2 = j-i2r
      jm3 = jm2-i2rby2
      if (nrod .eq. 0) go to 125
c
c     odd number of unknowns
c
      if (i2r .ne. 1) go to 117
      do 116 i=1,mr
         b(i) = q(i,j)
         q(i,j) = q(i,jm2)
  116 continue
      go to 123
  117 do 118 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  118 continue
      if (nrodpr .ne. 0) go to 120
      do 119 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)
  119 continue
      ip = ip-mr
      go to 122
  120 continue
      do 121 i=1,mr
         q(i,j) = q(i,j)-q(i,jm1)+q(i,jm2)
  121 continue
  122 if (lr .eq. 0) go to 123
      call cosgen (lr,1,fnum2,0.5,tcos(kr+1))
  123 continue
      call cosgen (kr,1,fnum2,0.5,tcos)
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 124 i=1,mr
         q(i,j) = q(i,j)+b(i)
  124 continue
      kr = kr+i2r
      go to 141
  125 continue
c
c     even number of unknowns
c
      jp1 = j+i2rby2
      jp2 = j+i2r
      if (i2r .ne. 1) go to 129
      do 126 i=1,mr
         b(i) = q(i,j)
  126 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      ip = 0
      ipstor = mr
      do 127 i=1,mr
         p(i) = b(i)
         b(i) = b(i)+q(i,n)
  127 continue
      tcos(1) = -1.+2*(np/2)
      tcos(2) = 0.
      call trix (1,1,mr,a,bb,c,b,tcos,d,w)
      do 128 i=1,mr
         q(i,j) = q(i,jm2)+p(i)+b(i)
  128 continue
      go to 140
  129 continue
      do 130 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  130 continue
      if (nrodpr .ne. 0) go to 132
      do 131 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  131 continue
      go to 134
  132 continue
      do 133 i=1,mr
         b(i) = b(i)+q(i,jp2)-q(i,jp1)
  133 continue
  134 continue
      call cosgen (i2r,1,0.5,0.0,tcos)
      call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
      ip = ip+mr
      ipstor = max(ipstor,ip+mr)
      do 135 i=1,mr
         ii = ip+i
         p(ii) = b(i)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
         b(i) = p(ii)+q(i,jp2)
  135 continue
      if (lr .eq. 0) go to 136
      call cosgen (lr,1,fnum2,0.5,tcos(i2r+1))
      call s1merg (tcos,0,i2r,i2r,lr,kr)
      go to 138
  136 do 137 i=1,i2r
         ii = kr+i
         tcos(ii) = tcos(i)
  137 continue
  138 call cosgen (kr,1,fnum2,0.5,tcos)
      call trix (kr,kr,mr,a,bb,c,b,tcos,d,w)
      do 139 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)+b(i)
  139 continue
  140 continue
      lr = kr
      kr = kr+jr
  141 continue
      nr = (nlast-1)/jr+1
      if (nr .le. 3) go to 142
      i2r = jr
      nrodpr = nrod
      go to 101
  142 continue
c
c      begin solution
c
      j = 1+jr
      jm1 = j-i2r
      jp1 = j+i2r
      jm2 = nlast-i2r
      if (nr .eq. 2) go to 180
      if (lr .ne. 0) go to 167
      if (n .ne. 3) go to 156
c
c     case n = 3.
c
      go to (143,148,143),np
  143 do 144 i=1,mr
         b(i) = q(i,2)
         b2(i) = q(i,1)+q(i,3)
         b3(i) = 0.
  144 continue
      go to (146,146,145),np
  145 tcos(1) = -1.
      tcos(2) = 1.
      k1 = 1
      go to 147
  146 tcos(1) = -2.
      tcos(2) = 1.
      tcos(3) = -1.
      k1 = 2
  147 k2 = 1
      k3 = 0
      k4 = 0
      go to 150
  148 do 149 i=1,mr
         b(i) = q(i,2)
         b2(i) = q(i,3)
         b3(i) = q(i,1)
  149 continue
      call cosgen (3,1,0.5,0.0,tcos)
      tcos(4) = -1.
      tcos(5) = 1.
      tcos(6) = -1.
      tcos(7) = 1.
      k1 = 3
      k2 = 2
      k3 = 1
      k4 = 1
  150 call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 151 i=1,mr
         b(i) = b(i)+b2(i)+b3(i)
  151 continue
      go to (153,153,152),np
  152 tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
  153 do 154 i=1,mr
         q(i,2) = b(i)
         b(i) = q(i,1)+b(i)
  154 continue
      tcos(1) = -1.+4.*fnum
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 155 i=1,mr
         q(i,1) = b(i)
  155 continue
      jr = 1
      i2r = 0
      go to 188
c
c     case n = 2**p+1
c
  156 continue
      do 157 i=1,mr
         b(i) = q(i,j)+q(i,1)-q(i,jm1)+q(i,nlast)-q(i,jm2)
  157 continue
      go to (158,160,158),np
  158 do 159 i=1,mr
         b2(i) = q(i,1)+q(i,nlast)+q(i,j)-q(i,jm1)-q(i,jp1)
         b3(i) = 0.
  159 continue
      k1 = nlast-1
      k2 = nlast+jr-1
      call cosgen (jr-1,1,0.0,1.0,tcos(nlast))
      tcos(k2) = 2*np-4
      call cosgen (jr,1,0.5-fnum,0.5,tcos(k2+1))
      k3 = (3-np)/2
      call s1merg (tcos,k1,jr-k3,k2-k3,jr+k3,0)
      k1 = k1-1+k3
      call cosgen (jr,1,fnum,0.5,tcos(k1+1))
      k2 = jr
      k3 = 0
      k4 = 0
      go to 162
  160 do 161 i=1,mr
         fi = (q(i,j)-q(i,jm1)-q(i,jp1))/2.
         b2(i) = q(i,1)+fi
         b3(i) = q(i,nlast)+fi
  161 continue
      k1 = nlast+jr-1
      k2 = k1+jr-1
      call cosgen (jr-1,1,0.0,1.0,tcos(k1+1))
      call cosgen (nlast,1,0.5,0.0,tcos(k2+1))
      call s1merg (tcos,k1,jr-1,k2,nlast,0)
      k3 = k1+nlast-1
      k4 = k3+jr
      call cosgen (jr,1,0.5,0.5,tcos(k3+1))
      call cosgen (jr,1,0.0,0.5,tcos(k4+1))
      call s1merg (tcos,k3,jr,k4,jr,k1)
      k2 = nlast-1
      k3 = jr
      k4 = jr
  162 call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 163 i=1,mr
         b(i) = b(i)+b2(i)+b3(i)
  163 continue
      if (np .ne. 3) go to 164
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
  164 do 165 i=1,mr
         q(i,j) = b(i)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
         b(i) = q(i,j)+q(i,1)
  165 continue
      call cosgen (jr,1,fnum,0.5,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 166 i=1,mr
         q(i,1) = q(i,1)-q(i,jm1)+b(i)
  166 continue
      go to 188
c
c     case of general n with nr = 3 .
c
  167 continue
      do 168 i=1,mr
         b(i) = q(i,1)-q(i,jm1)+q(i,j)
  168 continue
      if (nrod .ne. 0) go to 170
      do 169 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  169 continue
      go to 172
  170 do 171 i=1,mr
         b(i) = b(i)+q(i,nlast)-q(i,jm2)
  171 continue
  172 continue
      do 173 i=1,mr
         t = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
         q(i,j) = t
         b2(i) = q(i,nlast)+t
         b3(i) = q(i,1)+t
  173 continue
      k1 = kr+2*jr
      call cosgen (jr-1,1,0.0,1.0,tcos(k1+1))
      k2 = k1+jr
      tcos(k2) = 2*np-4
      k4 = (np-1)*(3-np)
      k3 = k2+1-k4
      call cosgen (kr+jr+k4,1,k4/2.,1.-k4,tcos(k3))
      k4 = 1-np/3
      call s1merg (tcos,k1,jr-k4,k2-k4,kr+jr+k4,0)
      if (np .eq. 3) k1 = k1-1
      k2 = kr+jr
      k4 = k1+k2
      call cosgen (kr,1,fnum2,0.5,tcos(k4+1))
      k3 = k4+kr
      call cosgen (jr,1,fnum,0.5,tcos(k3+1))
      call s1merg (tcos,k4,kr,k3,jr,k1)
      k4 = k3+jr
      call cosgen (lr,1,fnum2,0.5,tcos(k4+1))
      call s1merg (tcos,k3,jr,k4,lr,k1+k2)
      call cosgen (kr,1,fnum2,0.5,tcos(k3+1))
      k3 = kr
      k4 = kr
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 174 i=1,mr
         b(i) = b(i)+b2(i)+b3(i)
  174 continue
      if (np .ne. 3) go to 175
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
  175 do 176 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+q(i,j)
  176 continue
      call cosgen (jr,1,fnum,0.5,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      if (jr .ne. 1) go to 178
      do 177 i=1,mr
         q(i,1) = b(i)
  177 continue
      go to 188
  178 continue
      do 179 i=1,mr
         q(i,1) = q(i,1)-q(i,jm1)+b(i)
  179 continue
      go to 188
  180 continue
c
c     case of general n and nr = 2 .
c
      do 181 i=1,mr
         ii = ip+i
         b3(i) = 0.
         b(i) = q(i,1)+p(ii)
         q(i,1) = q(i,1)-q(i,jm1)
         b2(i) = q(i,1)+q(i,nlast)
  181 continue
      k1 = kr+jr
      k2 = k1+jr
      call cosgen (jr-1,1,0.0,1.0,tcos(k1+1))
      go to (182,183,182),np
  182 tcos(k2) = 2*np-4
      call cosgen (kr,1,0.0,1.0,tcos(k2+1))
      go to 184
  183 call cosgen (kr+1,1,0.5,0.0,tcos(k2))
  184 k4 = 1-np/3
      call s1merg (tcos,k1,jr-k4,k2-k4,kr+k4,0)
      if (np .eq. 3) k1 = k1-1
      k2 = kr
      call cosgen (kr,1,fnum2,0.5,tcos(k1+1))
      k4 = k1+kr
      call cosgen (lr,1,fnum2,0.5,tcos(k4+1))
      k3 = lr
      k4 = 0
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 185 i=1,mr
         b(i) = b(i)+b2(i)
  185 continue
      if (np .ne. 3) go to 186
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
  186 do 187 i=1,mr
         q(i,1) = q(i,1)+b(i)
  187 continue
  188 continue
c
c     start back substitution.
c
      j = nlast-jr
      do 189 i=1,mr
         b(i) = q(i,nlast)+q(i,j)
  189 continue
      jm2 = nlast-i2r
      if (jr .ne. 1) go to 191
      do 190 i=1,mr
         q(i,nlast) = 0.
  190 continue
      go to 195
  191 continue
      if (nrod .ne. 0) go to 193
      do 192 i=1,mr
         ii = ip+i
         q(i,nlast) = p(ii)
  192 continue
      ip = ip-mr
      go to 195
  193 do 194 i=1,mr
         q(i,nlast) = q(i,nlast)-q(i,jm2)
  194 continue
  195 continue
      call cosgen (kr,1,fnum2,0.5,tcos)
      call cosgen (lr,1,fnum2,0.5,tcos(kr+1))
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 196 i=1,mr
         q(i,nlast) = q(i,nlast)+b(i)
  196 continue
      nlastp = nlast
  197 continue
      jstep = jr
      jr = i2r
      i2r = i2r/2
      if (jr .eq. 0) go to 210
      jstart = 1+jr
      kr = kr-jr
      if (nlast+jr .gt. n) go to 198
      kr = kr-jr
      nlast = nlast+jr
      jstop = nlast-jstep
      go to 199
  198 continue
      jstop = nlast-jr
  199 continue
      lr = kr-jr
      call cosgen (jr,1,0.5,0.0,tcos)
      do 209 j=jstart,jstop,jstep
         jm2 = j-jr
         jp2 = j+jr
         if (j .ne. jr) go to 201
         do 200 i=1,mr
            b(i) = q(i,j)+q(i,jp2)
  200    continue
         go to 203
  201    continue
         do 202 i=1,mr
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  202    continue
  203    continue
         if (jr .ne. 1) go to 205
         do 204 i=1,mr
            q(i,j) = 0.
  204    continue
         go to 207
  205    continue
         jm1 = j-i2r
         jp1 = j+i2r
         do 206 i=1,mr
            q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
  206    continue
  207    continue
         call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
         do 208 i=1,mr
            q(i,j) = q(i,j)+b(i)
  208    continue
  209 continue
      nrod = 1
      if (nlast+i2r .le. n) nrod = 0
      if (nlastp .ne. nlast) go to 188
      go to 197
  210 continue
c
c     return storage requirements for p vectors.
c
      w(1) = ipstor
      return
      end
