*deck poisn2
      subroutine poisn2 (m, n, istag, mixbnd, a, bb, c, q, idimq, b, b2,
     +   b3, w, w2, w3, d, tcos, p)
c***begin prologue  poisn2
c***subsidiary
c***purpose  subsidiary to genbun
c***library   slatec
c***type      single precision (poisn2-s, cmposn-c)
c***author  (unknown)
c***description
c
c     subroutine to solve poisson's equation with neumann boundary
c     conditions.
c
c     istag = 1 if the last diagonal block is a.
c     istag = 2 if the last diagonal block is a-i.
c     mixbnd = 1 if have neumann boundary conditions at both boundaries.
c     mixbnd = 2 if have neumann boundary conditions at bottom and
c     dirichlet condition at top.  (for this case, must have istag = 1.)
c
c***see also  genbun
c***routines called  cosgen, s1merg, tri3, trix
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   920130  modified to use merge routine s1merg rather than deleted
c           routine merge.  (wrb)
c***end prologue  poisn2
c
      dimension       a(*)       ,bb(*)      ,c(*)       ,q(idimq,*) ,
     1                b(*)       ,b2(*)      ,b3(*)      ,w(*)       ,
     2                w2(*)      ,w3(*)      ,d(*)       ,tcos(*)    ,
     3                k(4)       ,p(*)
      equivalence     (k(1),k1)  ,(k(2),k2)  ,(k(3),k3)  ,(k(4),k4)
c***first executable statement  poisn2
      fistag = 3-istag
      fnum = 1./istag
      fden = 0.5*(istag-1)
      mr = m
      ip = -mr
      ipstor = 0
      i2r = 1
      jr = 2
      nr = n
      nlast = n
      kr = 1
      lr = 0
      go to (101,103),istag
  101 continue
      do 102 i=1,mr
         q(i,n) = .5*q(i,n)
  102 continue
      go to (103,104),mixbnd
  103 if (n .le. 3) go to 155
  104 continue
      jr = 2*i2r
      nrod = 1
      if ((nr/2)*2 .eq. nr) nrod = 0
      go to (105,106),mixbnd
  105 jstart = 1
      go to 107
  106 jstart = jr
      nrod = 1-nrod
  107 continue
      jstop = nlast-jr
      if (nrod .eq. 0) jstop = jstop-i2r
      call cosgen (i2r,1,0.5,0.0,tcos)
      i2rby2 = i2r/2
      if (jstop .ge. jstart) go to 108
      j = jr
      go to 116
  108 continue
c
c     regular reduction.
c
      do 115 j=jstart,jstop,jr
         jp1 = j+i2rby2
         jp2 = j+i2r
         jp3 = jp2+i2rby2
         jm1 = j-i2rby2
         jm2 = j-i2r
         jm3 = jm2-i2rby2
         if (j .ne. 1) go to 109
         jm1 = jp1
         jm2 = jp2
         jm3 = jp3
  109    continue
         if (i2r .ne. 1) go to 111
         if (j .eq. 1) jm2 = jp2
         do 110 i=1,mr
            b(i) = 2.*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  110    continue
         go to 113
  111    continue
         do 112 i=1,mr
            fi = q(i,j)
            q(i,j) = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = fi+q(i,j)-q(i,jm3)-q(i,jp3)
  112    continue
  113    continue
         call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
         do 114 i=1,mr
            q(i,j) = q(i,j)+b(i)
  114    continue
c
c     end of reduction for regular unknowns.
c
  115 continue
c
c     begin special reduction for last unknown.
c
      j = jstop+jr
  116 nlast = j
      jm1 = j-i2rby2
      jm2 = j-i2r
      jm3 = jm2-i2rby2
      if (nrod .eq. 0) go to 128
c
c     odd number of unknowns
c
      if (i2r .ne. 1) go to 118
      do 117 i=1,mr
         b(i) = fistag*q(i,j)
         q(i,j) = q(i,jm2)
  117 continue
      go to 126
  118 do 119 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  119 continue
      if (nrodpr .ne. 0) go to 121
      do 120 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)
  120 continue
      ip = ip-mr
      go to 123
  121 continue
      do 122 i=1,mr
         q(i,j) = q(i,j)-q(i,jm1)+q(i,jm2)
  122 continue
  123 if (lr .eq. 0) go to 124
      call cosgen (lr,1,0.5,fden,tcos(kr+1))
      go to 126
  124 continue
      do 125 i=1,mr
         b(i) = fistag*b(i)
  125 continue
  126 continue
      call cosgen (kr,1,0.5,fden,tcos)
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 127 i=1,mr
         q(i,j) = q(i,j)+b(i)
  127 continue
      kr = kr+i2r
      go to 151
  128 continue
c
c     even number of unknowns
c
      jp1 = j+i2rby2
      jp2 = j+i2r
      if (i2r .ne. 1) go to 135
      do 129 i=1,mr
         b(i) = q(i,j)
  129 continue
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      ip = 0
      ipstor = mr
      go to (133,130),istag
  130 do 131 i=1,mr
         p(i) = b(i)
         b(i) = b(i)+q(i,n)
  131 continue
      tcos(1) = 1.
      tcos(2) = 0.
      call trix (1,1,mr,a,bb,c,b,tcos,d,w)
      do 132 i=1,mr
         q(i,j) = q(i,jm2)+p(i)+b(i)
  132 continue
      go to 150
  133 continue
      do 134 i=1,mr
         p(i) = b(i)
         q(i,j) = q(i,jm2)+2.*q(i,jp2)+3.*b(i)
  134 continue
      go to 150
  135 continue
      do 136 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  136 continue
      if (nrodpr .ne. 0) go to 138
      do 137 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  137 continue
      go to 140
  138 continue
      do 139 i=1,mr
         b(i) = b(i)+q(i,jp2)-q(i,jp1)
  139 continue
  140 continue
      call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
      ip = ip+mr
      ipstor = max(ipstor,ip+mr)
      do 141 i=1,mr
         ii = ip+i
         p(ii) = b(i)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
         b(i) = p(ii)+q(i,jp2)
  141 continue
      if (lr .eq. 0) go to 142
      call cosgen (lr,1,0.5,fden,tcos(i2r+1))
      call s1merg (tcos,0,i2r,i2r,lr,kr)
      go to 144
  142 do 143 i=1,i2r
         ii = kr+i
         tcos(ii) = tcos(i)
  143 continue
  144 call cosgen (kr,1,0.5,fden,tcos)
      if (lr .ne. 0) go to 145
      go to (146,145),istag
  145 continue
      call trix (kr,kr,mr,a,bb,c,b,tcos,d,w)
      go to 148
  146 continue
      do 147 i=1,mr
         b(i) = fistag*b(i)
  147 continue
  148 continue
      do 149 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)+b(i)
  149 continue
  150 continue
      lr = kr
      kr = kr+jr
  151 continue
      go to (152,153),mixbnd
  152 nr = (nlast-1)/jr+1
      if (nr .le. 3) go to 155
      go to 154
  153 nr = nlast/jr
      if (nr .le. 1) go to 192
  154 i2r = jr
      nrodpr = nrod
      go to 104
  155 continue
c
c      begin solution
c
      j = 1+jr
      jm1 = j-i2r
      jp1 = j+i2r
      jm2 = nlast-i2r
      if (nr .eq. 2) go to 184
      if (lr .ne. 0) go to 170
      if (n .ne. 3) go to 161
c
c     case n = 3.
c
      go to (156,168),istag
  156 continue
      do 157 i=1,mr
         b(i) = q(i,2)
  157 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 158 i=1,mr
         q(i,2) = b(i)
         b(i) = 4.*b(i)+q(i,1)+2.*q(i,3)
  158 continue
      tcos(1) = -2.
      tcos(2) = 2.
      i1 = 2
      i2 = 0
      call trix (i1,i2,mr,a,bb,c,b,tcos,d,w)
      do 159 i=1,mr
         q(i,2) = q(i,2)+b(i)
         b(i) = q(i,1)+2.*q(i,2)
  159 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 160 i=1,mr
         q(i,1) = b(i)
  160 continue
      jr = 1
      i2r = 0
      go to 194
c
c     case n = 2**p+1
c
  161 continue
      go to (162,170),istag
  162 continue
      do 163 i=1,mr
         b(i) = q(i,j)+.5*q(i,1)-q(i,jm1)+q(i,nlast)-q(i,jm2)
  163 continue
      call cosgen (jr,1,0.5,0.0,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 164 i=1,mr
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
         b(i) = q(i,1)+2.*q(i,nlast)+4.*q(i,j)
  164 continue
      jr2 = 2*jr
      call cosgen (jr,1,0.0,0.0,tcos)
      do 165 i=1,jr
         i1 = jr+i
         i2 = jr+1-i
         tcos(i1) = -tcos(i2)
  165 continue
      call trix (jr2,0,mr,a,bb,c,b,tcos,d,w)
      do 166 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+2.*q(i,j)
  166 continue
      call cosgen (jr,1,0.5,0.0,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 167 i=1,mr
         q(i,1) = .5*q(i,1)-q(i,jm1)+b(i)
  167 continue
      go to 194
c
c     case of general n with nr = 3 .
c
  168 do 169 i=1,mr
         b(i) = q(i,2)
         q(i,2) = 0.
         b2(i) = q(i,3)
         b3(i) = q(i,1)
  169 continue
      jr = 1
      i2r = 0
      j = 2
      go to 177
  170 continue
      do 171 i=1,mr
         b(i) = .5*q(i,1)-q(i,jm1)+q(i,j)
  171 continue
      if (nrod .ne. 0) go to 173
      do 172 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  172 continue
      go to 175
  173 do 174 i=1,mr
         b(i) = b(i)+q(i,nlast)-q(i,jm2)
  174 continue
  175 continue
      do 176 i=1,mr
         t = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
         q(i,j) = t
         b2(i) = q(i,nlast)+t
         b3(i) = q(i,1)+2.*t
  176 continue
  177 continue
      k1 = kr+2*jr-1
      k2 = kr+jr
      tcos(k1+1) = -2.
      k4 = k1+3-istag
      call cosgen (k2+istag-2,1,0.0,fnum,tcos(k4))
      k4 = k1+k2+1
      call cosgen (jr-1,1,0.0,1.0,tcos(k4))
      call s1merg (tcos,k1,k2,k1+k2,jr-1,0)
      k3 = k1+k2+lr
      call cosgen (jr,1,0.5,0.0,tcos(k3+1))
      k4 = k3+jr+1
      call cosgen (kr,1,0.5,fden,tcos(k4))
      call s1merg (tcos,k3,jr,k3+jr,kr,k1)
      if (lr .eq. 0) go to 178
      call cosgen (lr,1,0.5,fden,tcos(k4))
      call s1merg (tcos,k3,jr,k3+jr,lr,k3-lr)
      call cosgen (kr,1,0.5,fden,tcos(k4))
  178 k3 = kr
      k4 = kr
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 179 i=1,mr
         b(i) = b(i)+b2(i)+b3(i)
  179 continue
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 180 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+2.*q(i,j)
  180 continue
      call cosgen (jr,1,0.5,0.0,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      if (jr .ne. 1) go to 182
      do 181 i=1,mr
         q(i,1) = b(i)
  181 continue
      go to 194
  182 continue
      do 183 i=1,mr
         q(i,1) = .5*q(i,1)-q(i,jm1)+b(i)
  183 continue
      go to 194
  184 continue
      if (n .ne. 2) go to 188
c
c     case  n = 2
c
      do 185 i=1,mr
         b(i) = q(i,1)
  185 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 186 i=1,mr
         q(i,1) = b(i)
         b(i) = 2.*(q(i,2)+b(i))*fistag
  186 continue
      tcos(1) = -fistag
      tcos(2) = 2.
      call trix (2,0,mr,a,bb,c,b,tcos,d,w)
      do 187 i=1,mr
         q(i,1) = q(i,1)+b(i)
  187 continue
      jr = 1
      i2r = 0
      go to 194
  188 continue
c
c     case of general n and nr = 2 .
c
      do 189 i=1,mr
         ii = ip+i
         b3(i) = 0.
         b(i) = q(i,1)+2.*p(ii)
         q(i,1) = .5*q(i,1)-q(i,jm1)
         b2(i) = 2.*(q(i,1)+q(i,nlast))
  189 continue
      k1 = kr+jr-1
      tcos(k1+1) = -2.
      k4 = k1+3-istag
      call cosgen (kr+istag-2,1,0.0,fnum,tcos(k4))
      k4 = k1+kr+1
      call cosgen (jr-1,1,0.0,1.0,tcos(k4))
      call s1merg (tcos,k1,kr,k1+kr,jr-1,0)
      call cosgen (kr,1,0.5,fden,tcos(k1+1))
      k2 = kr
      k4 = k1+k2+1
      call cosgen (lr,1,0.5,fden,tcos(k4))
      k3 = lr
      k4 = 0
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 190 i=1,mr
         b(i) = b(i)+b2(i)
  190 continue
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 191 i=1,mr
         q(i,1) = q(i,1)+b(i)
  191 continue
      go to 194
  192 do 193 i=1,mr
         b(i) = q(i,nlast)
  193 continue
      go to 196
  194 continue
c
c     start back substitution.
c
      j = nlast-jr
      do 195 i=1,mr
         b(i) = q(i,nlast)+q(i,j)
  195 continue
  196 jm2 = nlast-i2r
      if (jr .ne. 1) go to 198
      do 197 i=1,mr
         q(i,nlast) = 0.
  197 continue
      go to 202
  198 continue
      if (nrod .ne. 0) go to 200
      do 199 i=1,mr
         ii = ip+i
         q(i,nlast) = p(ii)
  199 continue
      ip = ip-mr
      go to 202
  200 do 201 i=1,mr
         q(i,nlast) = q(i,nlast)-q(i,jm2)
  201 continue
  202 continue
      call cosgen (kr,1,0.5,fden,tcos)
      call cosgen (lr,1,0.5,fden,tcos(kr+1))
      if (lr .ne. 0) go to 204
      do 203 i=1,mr
         b(i) = fistag*b(i)
  203 continue
  204 continue
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 205 i=1,mr
         q(i,nlast) = q(i,nlast)+b(i)
  205 continue
      nlastp = nlast
  206 continue
      jstep = jr
      jr = i2r
      i2r = i2r/2
      if (jr .eq. 0) go to 222
      go to (207,208),mixbnd
  207 jstart = 1+jr
      go to 209
  208 jstart = jr
  209 continue
      kr = kr-jr
      if (nlast+jr .gt. n) go to 210
      kr = kr-jr
      nlast = nlast+jr
      jstop = nlast-jstep
      go to 211
  210 continue
      jstop = nlast-jr
  211 continue
      lr = kr-jr
      call cosgen (jr,1,0.5,0.0,tcos)
      do 221 j=jstart,jstop,jstep
         jm2 = j-jr
         jp2 = j+jr
         if (j .ne. jr) go to 213
         do 212 i=1,mr
            b(i) = q(i,j)+q(i,jp2)
  212    continue
         go to 215
  213    continue
         do 214 i=1,mr
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  214    continue
  215    continue
         if (jr .ne. 1) go to 217
         do 216 i=1,mr
            q(i,j) = 0.
  216    continue
         go to 219
  217    continue
         jm1 = j-i2r
         jp1 = j+i2r
         do 218 i=1,mr
            q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
  218    continue
  219    continue
         call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
         do 220 i=1,mr
            q(i,j) = q(i,j)+b(i)
  220    continue
  221 continue
      nrod = 1
      if (nlast+i2r .le. n) nrod = 0
      if (nlastp .ne. nlast) go to 194
      go to 206
  222 continue
c
c     return storage requirements for p vectors.
c
      w(1) = ipstor
      return
      end
