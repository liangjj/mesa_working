*deck @(#)sch2h.f	1.1  11/30/90
      subroutine sch2h(grad,nbf,nob,nco,nao,nvirt,hmo,mix,icas,nmix,
     $     shmo,lshmo,lbhmo,asort,ncors)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)sch2h.f	1.1   11/30/90
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
      implicit real*8(a-h,o-z)
      dimension grad(nbf,2),mix(nob,2),hmo(2)
cc
      dimension shmo(*),lshmo(*),lbhmo(*),asort(*)
      common / number / zero,pt5,one,two,four,eight
c
      nocc  =nco+nao
      mstart=nocc+1
      mend  =nob
c
      lsort=2*nob*nob
      nob2=2*nob
c
      iix=0
      if(nco.eq.0)go to 165
      do 160 i=1,nco
c----------------------------------------------c
c     core - core
c----------------------------------------------c
c
         do 20 j=1,i
            xx=(grad(i,j)+grad(j,i))*pt5
            if(nao.eq.0)go to 35
c
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
c
            do 30 m=1,nao
               ii=mix(nco+m,i)
               jj=mix(nco+m,j)
               if(ii.eq.0.or.jj.eq.0)go to 30
cc
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=-xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=-xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)-xx
cc
 30         continue
 35         continue
            if(nvirt.eq.0)go to 20
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
            do 40 m=1,nvirt
               ii=mix(nocc+m,i)
               jj=mix(nocc+m,j)
               if(ii.eq.0.or.jj.eq.0)go to 40
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=-xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=-xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)-xx
cc
 40         continue
 45         continue
 20      continue
         if(nao.eq.0)go to 125
c-----------------------------------------------------c
c       core - active
c-----------------------------------------------------c
         do 120 j=1,nao
            xx=(grad(i,j+nco)+grad(j+nco,i))*pt5
            ncoj=nco+j
            if(icas.ne.0)go to 145
            if(j.eq.1) go to 135
            j1=j-1
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 130 k=1,j1
               ii=mix(ncoj,nco+k)
               jj=mix(nco+k,i)
               if(ii.eq.0.or.jj.eq.0)go to 130
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)+xx
cc
 130        continue
 135        continue
            if(j.eq.nao)go to 145
            j1=j+1
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 140 k=j1,nao
               ii=mix(nco+k,ncoj)
               jj=mix(nco+k,i)
               if(ii.eq.0.or.jj.eq.0)go to 140
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=-xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=-xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)-xx
cc
 140        continue
 145        continue
            if(nvirt.eq.0)go to 155
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 150 m=mstart,mend
               ii=mix(m,ncoj)
               jj=mix(m,i)
               if(ii.eq.0.or.jj.eq.0)go to 150
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=-xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=-xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)-xx
cc
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
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 230 m=1,nao
               ii=mix(nco+m,i)
               jj=mix(nocc+j,nco+m)
               if(ii.eq.0.or.jj.eq.0)go to 230
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)+xx
cc
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
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 310 m=istart,j1
               ii=mix(ia,m)
               jj=mix(ja,m)
               if(ii.eq.0.or.jj.eq.0)go to 310
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=-xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=-xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)-xx
cc
 310        continue
 315        continue
            i1=ia-1
            j1=ja+1
            if(i1.lt.j1)go to 325
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 320 m=j1,i1
               ii=mix(ia,m)
               jj=mix(m,ja)
               if(ii.eq.0.or.jj.eq.0)go to 320
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)+xx
cc
 320        continue
 325        continue
            if(ia.eq.nocc)go to 335
            i1=ia+1
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 330 m=i1,iend
               ii=mix(m,ia)
               jj=mix(m,ja)
               if(ii.eq.0.or.jj.eq.0)go to 330
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=-xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=-xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)-xx
cc
 330        continue
 335        continue
c
            if (nco .eq. 0) go to 345
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 340 m = 1, nco
               ii = mix(ia,m)
               jj = mix(ja,m)
               if (ii .eq. 0 .or. jj .eq. 0) go to 340
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=-xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=-xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)-xx
cc
 340        continue
 345        continue
            if (nvirt .eq. 0) go to 355
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 350 m = mstart, mend
               ii = mix(m,ia)
               jj = mix(m,ja)
               if (ii .eq. 0 .or. jj .eq. 0) go to 350
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=-xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=-xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)-xx
cc
 350        continue
 355        continue
 360     continue
c
         if (nvirt .eq. 0) go to 410
c
         do 400 ja = mstart, mend
            xx=pt5*grad(ja,ia)
            if (nco .eq. 0) go to 375
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 370 m = 1, nco
               ii = mix(ia,m)
               jj = mix(ja,m)
               if (ii .eq. 0 .or. jj .eq. 0) go to 370
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=-xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=-xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)-xx
cc
 370        continue
 375        continue
c
            if (nao .eq. 1) go to 385
            if (ia .eq. istart) go to 385
            i1 = ia - 1
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 380 m = istart, i1
               ii = mix(ia,m)
               jj = mix(ja,m)
               if (ii .eq. 0 .or. jj .eq. 0) go to 380
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=-xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=-xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)-xx
cc
 380        continue
 385        continue
            if (ia .eq. iend) go to 395
            i1 = ia + 1
cc
            if(iix+nob2.gt.lsort) then
               call sorter('with bin',asort,asort,0,iix,lshmo,
     $              lbhmo,shmo,0,0,0,.false.)
               iix=0
            endif
cc
            do 390 m = i1, iend
               ii = mix(m,ia)
               jj = mix(ja,m)
               if (ii .eq. 0 .or. jj .eq. 0) go to 390
               ipt=(ii-1)*nmix+jj
               iix=iix+1
               shmo(iix)=xx
               lshmo(iix)=ipt
               if(ii.ne.jj) then
                  jpt=(jj-1)*nmix+ii
                  iix=iix+1
                  shmo(iix)=xx
                  lshmo(iix)=jpt
               endif
cc
cc      hmo(ipt)=hmo(ipt)+xx
cc
 390        continue
 395        continue
 400     continue
 410  continue
 415  continue
cc
      if(iix.ne.0) then
         call sorter('with bin',asort,asort,0,iix,lshmo,
     $        lbhmo,shmo,0,0,0,.false.)
         iix=0
      endif
cc
      return
      end