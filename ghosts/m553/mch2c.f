      subroutine mch2c(grad,nbf,nob,nco,nao,nvirt,b,c,mix,icas)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c     implicit real*8(a-h,o-p,r-z),             integer*2(q)
cc
cmp   extended dummy grad,b,c,mix
cc
      dimension grad(nbf,2),b(2),c(2),mix(nob,2)
      common / number / zero,pt5,one,two,four,eight
c
      nocc  =nco+nao
      mstart=nocc+1
      mend  =nob
c
      if(nco.eq.0)go to 165
      do 160 i=1,nco
         i1=i-1
         if(i1.eq.0) go to 21
         do 20 j=1,i1
            xx=(grad(i,j)+grad(j,i))*pt5
            if(nao.eq.0)go to 35
            do 30 m=1,nao
               ii=mix(nco+m,i)
               jj=mix(nco+m,j)
               if(ii.eq.0.or.jj.eq.0)go to 30
               b(ii)=b(ii)-xx*c(jj)
               b(jj)=b(jj)-xx*c(ii)
 30         continue
 35         continue
            if(nvirt.eq.0)go to 20
            do 40 m=1,nvirt
               ii=mix(nocc+m,i)
               jj=mix(nocc+m,j)
               if(ii.eq.0.or.jj.eq.0)go to 40
               b(ii)=b(ii)-xx*c(jj)
               b(jj)=b(jj)-xx*c(ii)
 40         continue
 45         continue
 20      continue
 21      continue
c-------------------------------------------------c
c        diagonal contribution
c-------------------------------------------------c
         xx=grad(i,i)
         if(nao.eq.0)go to 535
         do 530 m=1,nao
            ii=mix(nco+m,i)
            if(ii.eq.0)go to 530
            b(ii)=b(ii)-xx*c(ii)
 530     continue
 535     continue
         if(nvirt.eq.0)go to 520
         do 540 m=1,nvirt
            ii=mix(nocc+m,i)
            if(ii.eq.0)go to 540
            b(ii)=b(ii)-xx*c(ii)
 540     continue
 545     continue
 520     continue
c
         if(nao.eq.0)go to 125
         do 120 j=1,nao
            xx=(grad(i,j+nco)+grad(j+nco,i))*pt5
            if(icas.ne.0)go to 145
            if(j.eq.1) go to 135
            j1=j-1
            ncoj=nco+j
            do 130 k=1,j1
               ii=mix(ncoj,nco+k)
               jj=mix(nco+k,i)
               if(ii.eq.0.or.jj.eq.0)go to 130
               b(ii)=b(ii)+xx*c(jj)
               b(jj)=b(jj)+xx*c(ii)
 130        continue
 135        continue
            if(j.eq.nao)go to 145
            ncoj=nco+j
            j1=j+1
            do 140 k=j1,nao
               ii=mix(nco+k,ncoj)
               jj=mix(nco+k,i)
               if(ii.eq.0.or.jj.eq.0)go to 140
               b(ii)=b(ii)-xx*c(jj)
               b(jj)=b(jj)-xx*c(ii)
 140        continue
 145        continue
            if(nvirt.eq.0)go to 155
            do 150 m=mstart,mend
               ii=mix(m,ncoj)
               jj=mix(m,i)
               if(ii.eq.0.or.jj.eq.0)go to 150
               b(ii)=b(ii)-xx*c(jj)
               b(jj)=b(jj)-xx*c(ii)
 150        continue
 155        continue
c
 120     continue
c-------------------------------------------------------c
c     core - virtual
c-------------------------------------------------------c
         if(nvirt.eq.0)go to 235
         do 220 j=1,nvirt
            xx=grad(nocc+j,i)*pt5
            do 230 m=1,nao
               ii=mix(nocc+m,i)
               jj=mix(nocc+j,nco+m)
               if(ii.eq.0.or.jj.eq.0)go to 230
               b(ii)=b(ii)+xx*c(jj)
               b(jj)=b(jj)+xx*c(ii)
 230        continue
 220     continue
 235     continue
c
 125     continue
 160  continue
 165  continue
c
c
      if(nao.eq.0)go to 415
      istart=nco+1
      iend  =nocc
      do 410 ia=istart,iend
         do 360 ja=istart,ia
            xx=(grad(ia,ja)+grad(ja,ia))*pt5
            if(ja.eq.istart)go to 315
            j1=ja-1
            do 310 m=istart,j1
               ii=mix(ia,m)
               jj=mix(ja,m)
               if(ii.eq.0.or.jj.eq.0)go to 310
               b(ii)=b(ii)-xx*c(jj)
               b(jj)=b(jj)-xx*c(ii)
 310        continue
 315        continue
            i1=ia-1
            j1=ja+1
            if(i1.lt.j1)go to 325
            do 320 m=j1,i1
               ii=mix(ia,m)
               jj=mix(m,ja)
               if(ii.eq.0.or.jj.eq.0)go to 320
               b(ii)=b(ii)-xx*c(jj)
               b(jj)=b(jj)-xx*c(ii)
 320        continue
 325        continue
            if(ia.eq.nocc)go to 335
            i1=ia+1
            do 330 m=i1,iend
               ii=mix(m,ia)
               jj=mix(m,ja)
               if(ii.eq.0.or.jj.eq.0)go to 330
               b(ii)=b(ii)-xx*c(jj)
               b(jj)=b(jj)-xx*c(ii)
 330        continue
 335        continue
c
            if (nco .eq. 0) go to 345
            do 340 m = 1, nco
               ii = mix(ia,m)
               jj = mix(ja,m)
               if (ii .eq. 0 .or. jj .eq. 0) go to 340
               b(ii) = b(ii) - xx * c(jj)
               b(jj) = b(jj) - xx * c(ii)
 340        continue
 345        continue
            if (nvirt .eq. 0) go to 355
            do 350 m = mstart, mend
               ii = mix(m,ia)
               jj = mix(m,ja)
               if (ii .eq. 0 .or. jj .eq. 0) go to 350
               b(ii) = b(ii) - xx * c(jj)
               b(jj) = b(jj) - xx * c(ii)
 350        continue
 355        continue
 360     continue
c
         if (nvirt .eq. 0) go to 410
c
         do 400 ja = mstart, mend
            xx=pt5*grad(ja,ia)
            if (nco .eq. 0) go to 375
            do 370 m = 1, nco
               ii = mix(ia,m)
               jj = mix(ja,m)
               if (ii .eq. 0 .or. jj .eq. 0) go to 370
               b(ii) = b(ii) - xx * c(jj)
               b(jj) = b(jj) - xx * c(ii)
 370        continue
 375        continue
c
            if (nao .eq. 1) go to 385
            if (ia .eq. istart) go to 385
            i1 = ia - 1
            do 380 m = istart, i1
               ii = mix(ia,m)
               jj = mix(ja,m)
               if (ii .eq. 0 .or. jj .eq. 0) go to 380
               b(ii) = b(ii) - xx * c(jj)
               b(jj) = b(jj) - xx * c(ii)
 380        continue
 385        continue
            if (ia .eq. iend) go to 395
            i1 = ia + 1
            do 390 m = i1, iend
               ii = mix(m,ia)
               jj = mix(ja,m)
               if (ii .eq. 0 .or. jj .eq. 0) go to 390
               b(ii) = b(ii) - xx * c(jj)
               b(jj) = b(jj) - xx * c(ii)
 390        continue
 395        continue
 400     continue
 410  continue
 415  continue
      return
      end
