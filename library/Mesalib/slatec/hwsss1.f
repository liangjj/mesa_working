*deck hwsss1
      subroutine hwsss1 (ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n,
     +   nbdcnd, bdps, bdpf, elmbda, f, idimf, pertrb, am, bm, cm, sn,
     +   ss, sint, d)
c***begin prologue  hwsss1
c***subsidiary
c***purpose  subsidiary to hwsssp
c***library   slatec
c***type      single precision (hwsss1-s)
c***author  (unknown)
c***see also  hwsssp
c***routines called  genbun
c***revision history  (yymmdd)
c   801001  date written
c   891009  removed unreferenced variables.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  hwsss1
      dimension       f(idimf,*) ,bdts(*)    ,bdtf(*)    ,bdps(*)    ,
     1                bdpf(*)    ,am(*)      ,bm(*)      ,cm(*)      ,
     2                ss(*)      ,sn(*)      ,d(*)       ,sint(*)
c
c***first executable statement  hwsss1
      mp1 = m+1
      np1 = n+1
      fn = n
      fm = m
      dth = (tf-ts)/fm
      hdth = dth/2.
      tdt = dth+dth
      dphi = (pf-ps)/fn
      tdp = dphi+dphi
      dphi2 = dphi*dphi
      dth2 = dth*dth
      cp = 4./(fn*dth2)
      wp = fn*sin(hdth)/4.
      do 102 i=1,mp1
         fim1 = i-1
         theta = fim1*dth+ts
         sint(i) = sin(theta)
         if (sint(i)) 101,102,101
  101    t1 = 1./(dth2*sint(i))
         am(i) = t1*sin(theta-hdth)
         cm(i) = t1*sin(theta+hdth)
         bm(i) = -am(i)-cm(i)+elmbda
  102 continue
      inp = 0
      isp = 0
c
c boundary condition at theta=ts
c
      mbr = mbdcnd+1
      go to (103,104,104,105,105,106,106,104,105,106),mbr
  103 its = 1
      go to 107
  104 at = am(2)
      its = 2
      go to 107
  105 at = am(1)
      its = 1
      cm(1) = am(1)+cm(1)
      go to 107
  106 at = am(2)
      inp = 1
      its = 2
c
c boundary condition theta=tf
c
  107 go to (108,109,110,110,109,109,110,111,111,111),mbr
  108 itf = m
      go to 112
  109 ct = cm(m)
      itf = m
      go to 112
  110 ct = cm(m+1)
      am(m+1) = am(m+1)+cm(m+1)
      itf = m+1
      go to 112
  111 itf = m
      isp = 1
      ct = cm(m)
c
c compute homogeneous solution with solution at pole equal to one
c
  112 itsp = its+1
      itfm = itf-1
      wts = sint(its+1)*am(its+1)/cm(its)
      wtf = sint(itf-1)*cm(itf-1)/am(itf)
      munk = itf-its+1
      if (isp) 116,116,113
  113 d(its) = cm(its)/bm(its)
      do 114 i=itsp,m
         d(i) = cm(i)/(bm(i)-am(i)*d(i-1))
  114 continue
      ss(m) = -d(m)
      iid = m-its
      do 115 ii=1,iid
         i = m-ii
         ss(i) = -d(i)*ss(i+1)
  115 continue
      ss(m+1) = 1.
  116 if (inp) 120,120,117
  117 sn(1) = 1.
      d(itf) = am(itf)/bm(itf)
      iid = itf-2
      do 118 ii=1,iid
         i = itf-ii
         d(i) = am(i)/(bm(i)-cm(i)*d(i+1))
  118 continue
      sn(2) = -d(2)
      do 119 i=3,itf
         sn(i) = -d(i)*sn(i-1)
  119 continue
c
c boundary conditions at phi=ps
c
  120 nbr = nbdcnd+1
      wps = 1.
      wpf = 1.
      go to (121,122,122,123,123),nbr
  121 jps = 1
      go to 124
  122 jps = 2
      go to 124
  123 jps = 1
      wps = .5
c
c boundary condition at phi=pf
c
  124 go to (125,126,127,127,126),nbr
  125 jpf = n
      go to 128
  126 jpf = n
      go to 128
  127 wpf = .5
      jpf = n+1
  128 jpsp = jps+1
      jpfm = jpf-1
      nunk = jpf-jps+1
      fjj = jpfm-jpsp+1
c
c scale coefficients for subroutine genbun
c
      do 129 i=its,itf
         cf = dphi2*sint(i)*sint(i)
         am(i) = cf*am(i)
         bm(i) = cf*bm(i)
         cm(i) = cf*cm(i)
  129 continue
      am(its) = 0.
      cm(itf) = 0.
      ising = 0
      go to (130,138,138,130,138,138,130,138,130,130),mbr
  130 go to (131,138,138,131,138),nbr
  131 if (elmbda) 138,132,132
  132 ising = 1
      sum = wts*wps+wts*wpf+wtf*wps+wtf*wpf
      if (inp) 134,134,133
  133 sum = sum+wp
  134 if (isp) 136,136,135
  135 sum = sum+wp
  136 sum1 = 0.
      do 137 i=itsp,itfm
         sum1 = sum1+sint(i)
  137 continue
      sum = sum+fjj*(sum1+wts+wtf)
      sum = sum+(wps+wpf)*sum1
      hne = sum
  138 go to (146,142,142,144,144,139,139,142,144,139),mbr
  139 if (nbdcnd-3) 146,140,146
  140 yhld = f(1,jps)-4./(fn*dphi*dth2)*(bdpf(2)-bdps(2))
      do 141 j=1,np1
         f(1,j) = yhld
  141 continue
      go to 146
  142 do 143 j=jps,jpf
         f(2,j) = f(2,j)-at*f(1,j)
  143 continue
      go to 146
  144 do 145 j=jps,jpf
         f(1,j) = f(1,j)+tdt*bdts(j)*at
  145 continue
  146 go to (154,150,152,152,150,150,152,147,147,147),mbr
  147 if (nbdcnd-3) 154,148,154
  148 yhld = f(m+1,jps)-4./(fn*dphi*dth2)*(bdpf(m)-bdps(m))
      do 149 j=1,np1
         f(m+1,j) = yhld
  149 continue
      go to 154
  150 do 151 j=jps,jpf
         f(m,j) = f(m,j)-ct*f(m+1,j)
  151 continue
      go to 154
  152 do 153 j=jps,jpf
         f(m+1,j) = f(m+1,j)-tdt*bdtf(j)*ct
  153 continue
  154 go to (159,155,155,157,157),nbr
  155 do 156 i=its,itf
         f(i,2) = f(i,2)-f(i,1)/(dphi2*sint(i)*sint(i))
  156 continue
      go to 159
  157 do 158 i=its,itf
         f(i,1) = f(i,1)+tdp*bdps(i)/(dphi2*sint(i)*sint(i))
  158 continue
  159 go to (164,160,162,162,160),nbr
  160 do 161 i=its,itf
         f(i,n) = f(i,n)-f(i,n+1)/(dphi2*sint(i)*sint(i))
  161 continue
      go to 164
  162 do 163 i=its,itf
         f(i,n+1) = f(i,n+1)-tdp*bdpf(i)/(dphi2*sint(i)*sint(i))
  163 continue
  164 continue
      pertrb = 0.
      if (ising) 165,176,165
  165 sum = wts*wps*f(its,jps)+wts*wpf*f(its,jpf)+wtf*wps*f(itf,jps)+
     1      wtf*wpf*f(itf,jpf)
      if (inp) 167,167,166
  166 sum = sum+wp*f(1,jps)
  167 if (isp) 169,169,168
  168 sum = sum+wp*f(m+1,jps)
  169 do 171 i=itsp,itfm
         sum1 = 0.
         do 170 j=jpsp,jpfm
            sum1 = sum1+f(i,j)
  170    continue
         sum = sum+sint(i)*sum1
  171 continue
      sum1 = 0.
      sum2 = 0.
      do 172 j=jpsp,jpfm
         sum1 = sum1+f(its,j)
         sum2 = sum2+f(itf,j)
  172 continue
      sum = sum+wts*sum1+wtf*sum2
      sum1 = 0.
      sum2 = 0.
      do 173 i=itsp,itfm
         sum1 = sum1+sint(i)*f(i,jps)
         sum2 = sum2+sint(i)*f(i,jpf)
  173 continue
      sum = sum+wps*sum1+wpf*sum2
      pertrb = sum/hne
      do 175 j=1,np1
         do 174 i=1,mp1
            f(i,j) = f(i,j)-pertrb
  174    continue
  175 continue
c
c scale right side for subroutine genbun
c
  176 do 178 i=its,itf
         cf = dphi2*sint(i)*sint(i)
         do 177 j=jps,jpf
            f(i,j) = cf*f(i,j)
  177    continue
  178 continue
      call genbun (nbdcnd,nunk,1,munk,am(its),bm(its),cm(its),idimf,
     1             f(its,jps),ierror,d)
      if (ising) 186,186,179
  179 if (inp) 183,183,180
  180 if (isp) 181,181,186
  181 do 182 j=1,np1
         f(1,j) = 0.
  182 continue
      go to 209
  183 if (isp) 186,186,184
  184 do 185 j=1,np1
         f(m+1,j) = 0.
  185 continue
      go to 209
  186 if (inp) 193,193,187
  187 sum = wps*f(its,jps)+wpf*f(its,jpf)
      do 188 j=jpsp,jpfm
         sum = sum+f(its,j)
  188 continue
      dfn = cp*sum
      dnn = cp*((wps+wpf+fjj)*(sn(2)-1.))+elmbda
      dsn = cp*(wps+wpf+fjj)*sn(m)
      if (isp) 189,189,194
  189 cnp = (f(1,1)-dfn)/dnn
      do 191 i=its,itf
         hld = cnp*sn(i)
         do 190 j=jps,jpf
            f(i,j) = f(i,j)+hld
  190    continue
  191 continue
      do 192 j=1,np1
         f(1,j) = cnp
  192 continue
      go to 209
  193 if (isp) 209,209,194
  194 sum = wps*f(itf,jps)+wpf*f(itf,jpf)
      do 195 j=jpsp,jpfm
         sum = sum+f(itf,j)
  195 continue
      dfs = cp*sum
      dss = cp*((wps+wpf+fjj)*(ss(m)-1.))+elmbda
      dns = cp*(wps+wpf+fjj)*ss(2)
      if (inp) 196,196,200
  196 csp = (f(m+1,1)-dfs)/dss
      do 198 i=its,itf
         hld = csp*ss(i)
         do 197 j=jps,jpf
            f(i,j) = f(i,j)+hld
  197    continue
  198 continue
      do 199 j=1,np1
         f(m+1,j) = csp
  199 continue
      go to 209
  200 rtn = f(1,1)-dfn
      rts = f(m+1,1)-dfs
      if (ising) 202,202,201
  201 csp = 0.
      cnp = rtn/dnn
      go to 205
  202 if (abs(dnn)-abs(dsn)) 204,204,203
  203 den = dss-dns*dsn/dnn
      rts = rts-rtn*dsn/dnn
      csp = rts/den
      cnp = (rtn-csp*dns)/dnn
      go to 205
  204 den = dns-dss*dnn/dsn
      rtn = rtn-rts*dnn/dsn
      csp = rtn/den
      cnp = (rts-dss*csp)/dsn
  205 do 207 i=its,itf
         hld = cnp*sn(i)+csp*ss(i)
         do 206 j=jps,jpf
            f(i,j) = f(i,j)+hld
  206    continue
  207 continue
      do 208 j=1,np1
         f(1,j) = cnp
         f(m+1,j) = csp
  208 continue
  209 if (nbdcnd) 212,210,212
  210 do 211 i=1,mp1
         f(i,jpf+1) = f(i,jps)
  211 continue
  212 return
      end
