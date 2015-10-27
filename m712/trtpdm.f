*deck @(#)trtpdm.f	5.1  11/6/94
c..bhl
c      subroutine trtpdm(sctpdm,prtpdm,iao,jao,kao,lao,ptdm,igrp,jgrp,
c     $     kgrp,lgrp,nfi,nfj,nfk,nfl,nci,ncj,nck,ncl,npi,npj,npk,npl,
c     $     ci,cj,ck,cl,lenblk,t1,t2,t3a,t3b,t4,t5,nckl,npij,npint,
c     $     aotpdm,nicmp,njcmp,nkcmp,nlcmp,flip,tdm,lentdm,dmoff,
c     $     n2pdm,printz,aointz,etot)
c..bhl
      subroutine trtpdm(sctpdm,prtpdm,iao,jao,kao,lao,ptdm,igrp,jgrp,
     $     kgrp,lgrp,nfi,nfj,nfk,nfl,nci,ncj,nck,ncl,npi,npj,npk,npl,
     $     ci,cj,ck,cl,lenblk,t1,t2,t3a,t3b,t4,t5,nckl,npij,npint,
     $     aotpdm,nicmp,njcmp,nkcmp,nlcmp,flip,tdm,lentdm,dmoff,
     $     n2pdm)
c

c
c***begin prologue     trtpdm
c***date written       871120   (yymmdd)
c***revision date      880105   (yymmdd)
c
c   05 january 1988    bhl at brl
c      fixed bugs in group-ordered ao density.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)trtpdm.f	5.1   11/6/94
c
c***purpose            to transform momentum groups of the two-particle
c    density matrix from the ao to primitive bases.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       trtpdm
c
      implicit integer (a-z)
c
      real*8 sctpdm(iao,jao,kao,lao)
      real*8 aotpdm(nicmp,njcmp,nkcmp,nlcmp)
      real*8 prtpdm(npint,lenblk)
c..bhl..test      real*8 aointz(lentdm),fac,egrp,etot,sdot
c..bhl..test      real*8 printz(iao,jao,kao,lao)
      real*8 ci(npi,nci)
      real*8 cj(npj,ncj)
      real*8 ck(npk,nck)
      real*8 cl(npl,ncl)
      real*8 t1(nci,ncj,nck,ncl)
      real*8 t2(npi,ncj,nckl)
      real*8 t3a(npi,npj,nckl)
      real*8 t3b(npij,nck,ncl)
      real*8 t4(npij,npk,ncl)
      real*8 t5(npint)
      real*8 tdm(lentdm)
c
      common /io/ inp,iout
c
c     ----- transfer the pertinent part of the density matrix,
c           reading from disk if necessary
c
      if (ptdm+iao*jao*kao*lao.gt.min(dmoff+lentdm,n2pdm)) then
         if (iao*jao*kao*lao.gt.lentdm) then
            call iosys('read real "group ordered 2pdm" from saoden',
     $           iao*jao*kao*lao,sctpdm,ptdm,' ')
c..bhl
c            call iosys('read real "group ordered 2e-ints" from saoden',
c     $           iao*jao*kao*lao,printz,ptdm,' ')
c..bhl
         else
            dmoff=ptdm
            lnread=min(lentdm,n2pdm-dmoff)
            call iosys('read real "group ordered 2pdm" from saoden',
     $           lnread,tdm,dmoff,' ')
            call vmove(sctpdm,tdm,iao*jao*kao*lao)
c..bhl
c            call iosys('read real "group ordered 2e-ints" from saoden',
c     $           lnread,aointz,dmoff,' ')
c            call vmove(printz,aointz,iao*jao*kao*lao)
c..bhl
         end if
      else
         call vmove(sctpdm,tdm(ptdm-dmoff+1),iao*jao*kao*lao)
c..bhl
c         call vmove(printz,aointz(ptdm-dmoff+1),iao*jao*kao*lao)
c..bhl
      end if
c
c     ----- for cases such as [ij;ij]
c
      if (igrp.eq.kgrp.and.jgrp.eq.lgrp) then
         do 12 i=1,iao
            do 11 j=1,jao
               do 10 k=1,i
c..bhl..jmxx
                  if(k.eq.i) then
                   jmxx=j
                  else
                   jmxx=jao
                  endif
                  do 9 l=1,jmxx
                     sctpdm(k,l,i,j)=sctpdm(i,j,k,l)
c..bhl..test                     printz(k,l,i,j)=printz(i,j,k,l)
 9                continue
 10            continue
 11         continue
 12      continue
      end if
c
c     ----- for cases such as [ii;kl], we only have the canonical
c           portion, and must fill out.
c
      if (igrp.eq.jgrp) then
         do 4 i=1,iao
            do 3 j=1,i-1
               do 2 k=1,kao
                  do 1 l=1,lao
                     sctpdm(j,i,k,l)=sctpdm(i,j,k,l)
c..bhl..test                     printz(j,i,k,l)=printz(i,j,k,l)
 1                continue
 2             continue
 3          continue
 4       continue
      end if
c
      if (kgrp.eq.lgrp) then
         do 8 k=1,kao
            do 7 l=1,k-1
               do 6 i=1,iao
                  do 5 j=1,jao
                     sctpdm(i,j,l,k)=sctpdm(i,j,k,l)
c..bhl..test                     printz(i,j,l,k)=printz(i,j,k,l)
 5                continue
 6             continue
 7          continue
 8       continue
      end if
c
c..bhl..old       if (igrp.eq.kgrp.and.jgrp.eq.lgrp) then
c
c
c$$$
c$$$     do 103 k=1,kao
c$$$        do 102 l=1,lao
c$$$           write (iout,101) k,l
c$$$101        format (/10x,'ao 2pdm for k,l=',2i4)
c$$$           call matout(sctpdm(1,1,k,l),iao,jao,iao,jao,iout)
c$$$           write (iout,111) k,l
c$$$111        format (/10x,'ao ints for k,l=',2i4)
c$$$           call matout(printz(1,1,k,l),iao,jao,iao,jao,iout)
c$$$102     continue
c$$$103  continue
c$$$
c
c ---- test code to trace group ordered density and integrals
c
c$$$  if (igrp.ne.jgrp) then
c$$$     fac=2.0d+00
c$$$  else
c$$$     fac=1.0d+00
c$$$  end if
c
c$$$  if (kgrp.ne.lgrp) then
c$$$     fac=fac*2.0d+00
c$$$  end if
c
c$$$  if (igrp.ne.kgrp.or.jgrp.ne.lgrp) then
c$$$     fac=fac*2.0d+00
c$$$  end if
c
c$$$  grptot=iao*jao*kao*lao
c$$$  egrp=sdot(grptot,sctpdm,1,printz,1)
c$$$  etot=etot+egrp*fac
c
c
c    ----- switch the density matrix around to correspond to
c              integral label order for derivatives
c
      if (flip.eq.1) then
         call vmove(aotpdm,sctpdm,iao*jao*kao*lao)
      else if (flip.eq.2) then
         do 54 l=1,lao
            do 53 k=1,kao
               do 52 j=1,jao
                  do 51 i=1,iao
                     aotpdm(j,i,k,l)=sctpdm(i,j,k,l)
 51               continue
 52            continue
 53         continue
 54      continue
      else if (flip.eq.3) then
         do 58 l=1,lao
            do 57 k=1,kao
               do 56 j=1,jao
                  do 55 i=1,iao
                     aotpdm(i,j,l,k)=sctpdm(i,j,k,l)
 55               continue
 56            continue
 57         continue
 58      continue
      else if (flip.eq.4) then
         do 64 l=1,lao
            do 63 k=1,kao
               do 62 j=1,jao
                  do 61 i=1,iao
                     aotpdm(j,i,l,k)=sctpdm(i,j,k,l)
 61               continue
 62            continue
 63         continue
 64      continue
      else if (flip.eq.5) then
         do 68 l=1,lao
            do 67 k=1,kao
               do 66 j=1,jao
                  do 65 i=1,iao
                     aotpdm(k,l,i,j)=sctpdm(i,j,k,l)
 65                  continue
 66            continue
 67         continue
 68      continue
      else if (flip.eq.6) then
         do 74 l=1,lao
            do 73 k=1,kao
               do 72 j=1,jao
                  do 71 i=1,iao
                     aotpdm(k,l,j,i)=sctpdm(i,j,k,l)
 71               continue
 72            continue
 73         continue
 74      continue
      else if (flip.eq.7) then
         do 78 l=1,lao
            do 77 k=1,kao
               do 76 j=1,jao
                  do 75 i=1,iao
                     aotpdm(l,k,i,j)=sctpdm(i,j,k,l)
 75               continue
 76            continue
 77         continue
 78      continue
      else if (flip.eq.8) then
         do 88 l=1,lao
            do 87 k=1,kao
               do 86 j=1,jao
                  do 85 i=1,iao
                     aotpdm(l,k,j,i)=sctpdm(i,j,k,l)
 85               continue
 86            continue
 87         continue
 88      continue
      end if
c
c     ----- reorder the density matrix from group order to integral
c      ao(fl,fk,fj,fi,ic,jc,kc,lc)=sc(fi,ic,fj,jc,fk,kc,fl,lc)
c
      n=0
      do 30 if=1,nfi
         do 29 jf=1,nfj
            do 28 kf=1,nfk
               do 27 lf=1,nfl
                  n=n+1
c
c                 ----- fetch the contracted block for these functions
c
                  do 16 lc=1,ncl
                     lcmp=(lc-1)*nfl+lf
                     do 15 kc=1,nck
                        kcmp=(kc-1)*nfk+kf
                        do 14 jc=1,ncj
                           jcmp=(jc-1)*nfj+jf
                           do 13 ic=1,nci
                              icmp=(ic-1)*nfi+if
                              t1(ic,jc,kc,lc)=
     $                             aotpdm(icmp,jcmp,kcmp,lcmp)
 13                        continue
 14                     continue
 15                  continue
 16               continue
c
c                 ----- transform this block -----
c
                  call ebc(t2,ci,t1,npi,nci,ncj*nck*ncl)
c
                  do 20 i=1,npi
                     do 19 j=1,npj
                        do 46 l=1,nck*ncl
                           t3a(i,j,l)=0.0d+00
 46                     continue
                        do 18 k=1,ncj
                           do 17 l=1,nck*ncl
                              t3a(i,j,l)=t3a(i,j,l)+t2(i,k,l)*cj(j,k)
 17                        continue
 18                     continue
 19                  continue
 20               continue
c
                  do 25 l=1,ncl
                     do 24 k=1,npk
                        do 21 i=1,npi*npj
                           t4(i,k,l)=0.0d+00
 21                     continue
                        do 23 j=1,nck
                           do 22 i=1,npi*npj
                              t4(i,k,l)=t4(i,k,l)+t3b(i,j,l)*ck(k,j)
 22                        continue
 23                     continue
 24                  continue
 25               continue
c
                  call ebct(t5,t4,cl,npi*npj*npk,ncl,npl)
c
c                 ----- and put them where they belong -----
c
                  do 26 i=1,npi*npj*npk*npl
                     prtpdm(i,n)=t5(i)
 26               continue
c
c
 27            continue
 28         continue
 29      continue
 30   continue
c
c
c$$$      write (iout,31)
c$$$ 31   format (10x,'primitive two-particle density matrix block')
c$$$      call matout(prtpdm,npint,lenblk,npint,lenblk,iout)
c
c
      return
      end
