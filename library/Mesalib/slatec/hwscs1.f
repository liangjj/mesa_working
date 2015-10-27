*deck hwscs1
      subroutine hwscs1 (intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n,
     +   nbdcnd, bdrs, bdrf, elmbda, f, idimf, pertrb, w, s, an, bn, cn,
     +   r, am, bm, cm, sint, bmh)
c***begin prologue  hwscs1
c***subsidiary
c***purpose  subsidiary to hwscsp
c***library   slatec
c***type      single precision (hwscs1-s)
c***author  (unknown)
c***see also  hwscsp
c***routines called  blktri
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891009  removed unreferenced variables.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  hwscs1
      dimension       f(idimf,*) ,bdrs(*)    ,bdrf(*)    ,bdts(*)    ,
     1                bdtf(*)    ,am(*)      ,bm(*)      ,cm(*)      ,
     2                an(*)      ,bn(*)      ,cn(*)      ,s(*)       ,
     3                r(*)       ,sint(*)    ,bmh(*)     ,w(*)
c***first executable statement  hwscs1
      mp1 = m+1
      dth = (tf-ts)/m
      tdt = dth+dth
      hdth = dth/2.
      sdts = 1./(dth*dth)
      do 102 i=1,mp1
         theta = ts+(i-1)*dth
         sint(i) = sin(theta)
         if (sint(i)) 101,102,101
  101    t1 = sdts/sint(i)
         am(i) = t1*sin(theta-hdth)
         cm(i) = t1*sin(theta+hdth)
         bm(i) = -(am(i)+cm(i))
  102 continue
      np1 = n+1
      dr = (rf-rs)/n
      hdr = dr/2.
      tdr = dr+dr
      dr2 = dr*dr
      czr = 6.*dth/(dr2*(cos(ts)-cos(tf)))
      do 103 j=1,np1
         r(j) = rs+(j-1)*dr
         an(j) = (r(j)-hdr)**2/dr2
         cn(j) = (r(j)+hdr)**2/dr2
         bn(j) = -(an(j)+cn(j))
  103 continue
      mp = 1
      np = 1
c
c boundary condition at phi=ps
c
      go to (104,104,105,105,106,106,104,105,106),mbdcnd
  104 at = am(2)
      its = 2
      go to 107
  105 at = am(1)
      its = 1
      cm(1) = cm(1)+am(1)
      go to 107
  106 its = 1
      bm(1) = -4.*sdts
      cm(1) = -bm(1)
c
c boundary condition at phi=pf
c
  107 go to (108,109,109,108,108,109,110,110,110),mbdcnd
  108 ct = cm(m)
      itf = m
      go to 111
  109 ct = cm(m+1)
      am(m+1) = am(m+1)+cm(m+1)
      itf = m+1
      go to 111
  110 itf = m+1
      am(m+1) = 4.*sdts
      bm(m+1) = -am(m+1)
  111 wts = sint(its+1)*am(its+1)/cm(its)
      wtf = sint(itf-1)*cm(itf-1)/am(itf)
      itsp = its+1
      itfm = itf-1
c
c boundary condition at r=rs
c
      ictr = 0
      go to (112,112,113,113,114,114),nbdcnd
  112 ar = an(2)
      jrs = 2
      go to 118
  113 ar = an(1)
      jrs = 1
      cn(1) = cn(1)+an(1)
      go to 118
  114 jrs = 2
      ictr = 1
      s(n) = an(n)/bn(n)
      do 115 j=3,n
         l = n-j+2
         s(l) = an(l)/(bn(l)-cn(l)*s(l+1))
  115 continue
      s(2) = -s(2)
      do 116 j=3,n
         s(j) = -s(j)*s(j-1)
  116 continue
      wtnm = wts+wtf
      do 117 i=itsp,itfm
         wtnm = wtnm+sint(i)
  117 continue
      yps = czr*wtnm*(s(2)-1.)
c
c boundary condition at r=rf
c
  118 go to (119,120,120,119,119,120),nbdcnd
  119 cr = cn(n)
      jrf = n
      go to 121
  120 cr = cn(n+1)
      an(n+1) = an(n+1)+cn(n+1)
      jrf = n+1
  121 wrs = an(jrs+1)*r(jrs)**2/cn(jrs)
      wrf = cn(jrf-1)*r(jrf)**2/an(jrf)
      wrz = an(jrs)/czr
      jrsp = jrs+1
      jrfm = jrf-1
      munk = itf-its+1
      nunk = jrf-jrs+1
      do 122 i=its,itf
         bmh(i) = bm(i)
  122 continue
      ising = 0
      go to (132,132,123,132,132,123),nbdcnd
  123 go to (132,132,124,132,132,124,132,124,124),mbdcnd
  124 if (elmbda) 132,125,125
  125 ising = 1
      sum = wts*wrs+wts*wrf+wtf*wrs+wtf*wrf
      if (ictr) 126,127,126
  126 sum = sum+wrz
  127 do 129 j=jrsp,jrfm
         r2 = r(j)**2
         do 128 i=itsp,itfm
            sum = sum+r2*sint(i)
  128    continue
  129 continue
      do 130 j=jrsp,jrfm
         sum = sum+(wts+wtf)*r(j)**2
  130 continue
      do 131 i=itsp,itfm
         sum = sum+(wrs+wrf)*sint(i)
  131 continue
      hne = sum
  132 go to (133,133,133,133,134,134,133,133,134),mbdcnd
  133 bm(its) = bmh(its)+elmbda/sint(its)**2
  134 go to (135,135,135,135,135,135,136,136,136),mbdcnd
  135 bm(itf) = bmh(itf)+elmbda/sint(itf)**2
  136 do 137 i=itsp,itfm
         bm(i) = bmh(i)+elmbda/sint(i)**2
  137 continue
      go to (138,138,140,140,142,142,138,140,142),mbdcnd
  138 do 139 j=jrs,jrf
         f(2,j) = f(2,j)-at*f(1,j)/r(j)**2
  139 continue
      go to 142
  140 do 141 j=jrs,jrf
         f(1,j) = f(1,j)+tdt*bdts(j)*at/r(j)**2
  141 continue
  142 go to (143,145,145,143,143,145,147,147,147),mbdcnd
  143 do 144 j=jrs,jrf
         f(m,j) = f(m,j)-ct*f(m+1,j)/r(j)**2
  144 continue
      go to 147
  145 do 146 j=jrs,jrf
         f(m+1,j) = f(m+1,j)-tdt*bdtf(j)*ct/r(j)**2
  146 continue
  147 go to (151,151,153,153,148,148),nbdcnd
  148 if (mbdcnd-3) 155,149,155
  149 yhld = f(its,1)-czr/tdt*(sin(tf)*bdtf(2)-sin(ts)*bdts(2))
      do 150 i=1,mp1
         f(i,1) = yhld
  150 continue
      go to 155
  151 rs2 = (rs+dr)**2
      do 152 i=its,itf
         f(i,2) = f(i,2)-ar*f(i,1)/rs2
  152 continue
      go to 155
  153 do 154 i=its,itf
         f(i,1) = f(i,1)+tdr*bdrs(i)*ar/rs**2
  154 continue
  155 go to (156,158,158,156,156,158),nbdcnd
  156 rf2 = (rf-dr)**2
      do 157 i=its,itf
         f(i,n) = f(i,n)-cr*f(i,n+1)/rf2
  157 continue
      go to 160
  158 do 159 i=its,itf
         f(i,n+1) = f(i,n+1)-tdr*bdrf(i)*cr/rf**2
  159 continue
  160 continue
      pertrb = 0.
      if (ising) 161,170,161
  161 sum = wts*wrs*f(its,jrs)+wts*wrf*f(its,jrf)+wtf*wrs*f(itf,jrs)+
     1      wtf*wrf*f(itf,jrf)
      if (ictr) 162,163,162
  162 sum = sum+wrz*f(its,1)
  163 do 165 j=jrsp,jrfm
         r2 = r(j)**2
         do 164 i=itsp,itfm
            sum = sum+r2*sint(i)*f(i,j)
  164    continue
  165 continue
      do 166 j=jrsp,jrfm
         sum = sum+r(j)**2*(wts*f(its,j)+wtf*f(itf,j))
  166 continue
      do 167 i=itsp,itfm
         sum = sum+sint(i)*(wrs*f(i,jrs)+wrf*f(i,jrf))
  167 continue
      pertrb = sum/hne
      do 169 j=1,np1
         do 168 i=1,mp1
            f(i,j) = f(i,j)-pertrb
  168    continue
  169 continue
  170 do 172 j=jrs,jrf
         rsq = r(j)**2
         do 171 i=its,itf
            f(i,j) = rsq*f(i,j)
  171    continue
  172 continue
      iflg = intl
  173 call blktri (iflg,np,nunk,an(jrs),bn(jrs),cn(jrs),mp,munk,
     1             am(its),bm(its),cm(its),idimf,f(its,jrs),ierror,w)
      iflg = iflg+1
      if (iflg-1) 174,173,174
  174 if (nbdcnd) 177,175,177
  175 do 176 i=1,mp1
         f(i,jrf+1) = f(i,jrs)
  176 continue
  177 if (mbdcnd) 180,178,180
  178 do 179 j=1,np1
         f(itf+1,j) = f(its,j)
  179 continue
  180 xp = 0.
      if (ictr) 181,188,181
  181 if (ising) 186,182,186
  182 sum = wts*f(its,2)+wtf*f(itf,2)
      do 183 i=itsp,itfm
         sum = sum+sint(i)*f(i,2)
  183 continue
      yph = czr*sum
      xp = (f(its,1)-yph)/yps
      do 185 j=jrs,jrf
         xps = xp*s(j)
         do 184 i=its,itf
            f(i,j) = f(i,j)+xps
  184    continue
  185 continue
  186 do 187 i=1,mp1
         f(i,1) = xp
  187 continue
  188 return
      end
