*deck @(#)mcmix.f	1.1  11/30/90
      subroutine mcmix(mix,mixo,mixhes,mixinv,mixl,mixu,len,lok,
     $     nocc,nob,nco,nao,nvo,iobsym,iobcas,icas,nmix,nmixs,
     $     nbf,isym)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcmix.f	1.1   11/30/90
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
cc
cmp   extended dummy mix,mixo,mixhes,mixinv,mixl,mixu,len,lok
cc
      dimension iobsym(2),iobcas(2)
c
      dimension mixhes(nob,2),mix(2),mixl(2),mixu(2),len(2),lok(2)
      dimension mixo(2),mixinv(nbf,2)
      common / number / zero,pt5,one,two,four,eight
c
cc
c    len1 and loc1  are not used   yet
cc
      do 2 i=1,nocc
         do 1 j=1,nob
            mixhes(j,i)=0
            mixinv(j,i)=0
 1       continue
 2    continue
c
      nmixs=0
c
      if(nco.eq.0)go to 31
c
c--------------------------------------
c  core-active and core-virtual mixings
c--------------------------------------
c
      do 30 i=1,nco
         iobsi=iobsym(i)
         if(nao.eq.0)go to 11
c-------------
c  core-active
c-------------
         jstart=nco+1
         jend=nocc
c
         do 10 j=jstart,jend
            if (iobsym(j) .ne. iobsi) go to 10
            nmix=nmix+1
            nmixs=nmixs+1
            mixhes(j,i)=nmix
            mixhes(i,j)=nmix
            mixl(nmixs)=i
            mixu(nmixs)=j
 10      continue
 11      continue
c
         if(nvo.eq.0)go to 21
c--------------
c  core-virtual
c--------------
         jstart=nocc+1
         jend=nob
c
         do 20 j=jstart,jend
            if (iobsym(j) .ne. iobsi) go to 20
            nmix=nmix+1
            nmixs=nmixs+1
            mixhes(j,i)=nmix
            mixl(nmixs)=i
            mixu(nmixs)=j
 20      continue
 21      continue
c
 30   continue
 31   continue
c
c-------------------------------------------
c   active-active and active-virtual mixings
c-------------------------------------------
c
      if(nao.eq.0)go to 61
c
      istart=nco+1
      iend=nocc
      icao=0
      do 60 i=istart,iend
         iobsi=iobsym(i)
         icao=icao+1
         iobci=iobcas(icao)
         if(icas.ne.0)go to 41
         if(nao.lt.2)go to 41
         if(i.eq.iend)go to 41
c---------------
c  active-active
c---------------
c
         jstart=i+1
         jend=nocc
c
         if(nao.eq.1)go to 41
c
         jcao=icao
         do 40 j=jstart,jend
            jcao=jcao+1
            if (iobsym(j) .ne. iobsi) go to 40
            if (iobci .ne. 0 .and. iobci .eq. iobcas(jcao)) go to 40
            nmix=nmix+1
            nmixs=nmixs+1
            mixhes(j,i)=nmix
            mixhes(i,j)=nmix
            mixl(nmixs)=i
            mixu(nmixs)=j
 40      continue
 41      continue
         if(nvo.eq.0)go to 51
c----------------
c  active-virtual
c----------------
c
         jstart=nocc+1
         jend=nob
c
         do 50 j=jstart,jend
            if (iobsym(j) .ne. iobsi) go to 50
            nmix=nmix+1
            nmixs=nmixs+1
            mixhes(j,i)=nmix
            mixl(nmixs)=i
            mixu(nmixs)=j
 50      continue
 51      continue
c
 60   continue
 61   continue
c
      ix=0
      ic=0
      do 80 i=1,nocc
         lok(i)=ix
         leni=0
         do 70 j=1,nob
            if(mixhes(j,i).eq.0)go to 70
            leni=leni+1
            mixinv(j,i)=leni
            ic=ic+1
            mix(ic)=mixhes(j,i)
            mixo(ic)=j
 70      continue
         len(i)=leni
         loci=lok(i)
c     write(iout,9001) leni,lok(i)
c     write(iout,9002) (mix(loci+mm),mm=1,leni)
c     write(iout,9003) (mixo(loci+mm),mm=1,leni)
         ix=ix+leni
 80   continue
c
c9001 format(/, '  leni  loci ',2i6)
c9002 format('  *** mix  ',(/, 1x,15i4))
c9003 format('  *** mixo ',(/,1x,15i4))
c
c     write(iout,9004) ((mixhes(j,i),j=1,nob),i=1,nocc)
c9004 format(' mixhes ',8i5)
      return
      end
