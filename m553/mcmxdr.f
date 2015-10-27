*deck @(#)mcmxdr.f	5.1  11/6/94
      subroutine mcmxdr(mix,mixo,mixhes,mixinv,mixl,mixu,
     $     len,lok,locp,nocc,nob,nco,nao,nvo,iobsym,iobcas,icas,nmix,
     $     nsym, isymm,locsym, nmixs, nbf,prtflg)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcmxdr.f	5.1   11/6/94
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
cmp   extended dummy mix,mixo,mixhes,mixinv,mixl,mixu,len,lok,locp
cc
      character*8 prtflg
c
      dimension mix(2),mixo(2),mixhes(2),mixinv(2),len(2),lok(2),
     $     mixl(2),mixu(2),nocc(2),nob(2),nao(2),nvo(2),icas(2),
     $     nmixs(2),
     $     isymm(2),locsym(2),nco(2),locp(2),iobsym(2),iobcas(2),nbf(2)
c
      common /io/ inp,iout
c
c
      nmix=0
      if (prtflg.ne.'minimum') then
         write(iout,19)
 19      format('1************ orbital mixing section ************')
      end if
      iao = 1
      iobt = 1
      do 20 i=1,nsym
         if(nocc(i).eq.0)go to 20
         lsymm=isymm(i)
         llocs=locsym(i)
         call mcmix(mix(lsymm+1),mixo(lsymm+1),
     $        mixhes(lsymm+1),mixinv(lsymm+1),
     $        mixl(lsymm+1), mixu(lsymm+1), len(llocs+1),
     $        lok(llocs+1),
     $        nocc(i), nob(i), nco(i), nao(i), nvo(i),iobsym(iobt),
     $        iobcas(iao),icas(i), nmix, nmixss, nbf(i), i )
         nmixs(i)=nmixss
         iob1 = locsym(i) + 1
         iob2 = locsym(i) + nco(i) + nao(i)
         do 30 iob = iob1, iob2
            locp(iob) = lok(iob) + isymm(i)
 30      continue
         iobt = iobt + nob(i)
         iao = iao + nao(i)
 20   continue
c
c     write(iout,1022) (isymm(i),locsym(i),nmixs(i),i=1,nsym)
c1022 format(' isymm locsym nmixs ',3(3x,3i4))
c     write(iout,21) nmix
c  21 format(//,'  *** mixhes ***   nmix = ',i5,/)
      if (prtflg.ne.'minimum') then
         call mcprti(nsym,nob,nocc,isymm,mixhes)
      end if
c
      return
      end
