*deck @(#)khnvec.f	5.1  11/6/94
      subroutine khnvec(c,eigval,orbsym,temp,eigtmp,tmpsym,shlmin,
     1                  shlmax,nbf,nshell,nirrep,nsmall,noabort)
c
c***purpose: to put the scf vector on rwf in order required
c            for kohn calculation
c
c barry schneider lanl
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf), eigval(nbf), temp(nbf,nbf), eigtmp(nbf)
      real*8 p
      integer shlmin(nirrep), shlmax(nirrep), orbsym(nbf)
      integer tmpsym(nbf,2)
      character *1 itoc
      common /io/ inp, iout
c     
c
      count=0
      if (nsmall.ne.0) then
          call scopy(nbf*nsmall,c(1,1),1,temp(1,1),1)
          call scopy(nsmall,eigval,1,eigtmp,1)
          call icopy(orbsym,tmpsym,nsmall)
          count=nsmall
      endif
      do 100 irrep=1,nirrep
         icirep=count
         do 200 shell=1,nshell
            do 300 i=shlmin(shell),shlmax(shell)           
               if (i.gt.nsmall) then
                   if(orbsym(i).eq.irrep) then
                      count=count+1
                      call scopy(nbf,c(1,i),1,temp(1,count),1)
                      eigtmp(count)=eigval(i)
                      tmpsym(count,1)=orbsym(i)     
                   endif   
               endif
  300       continue
  200    continue
c----------------------------------------------------------------------c
c           within this irrep re-order by eigenvalue                   c
c----------------------------------------------------------------------c
         if (count.gt.icirep+1) then
             do 400 ii=icirep+2,count
                i=ii-1
                k=i
                p=eigtmp(i)
                do 500 j=ii,count
                   if (eigtmp(j).lt.p) then
                       k=j
                       p=eigtmp(j)
                   endif
  500           continue 
                if (k.ne.i) then
                    eigtmp(k)=eigtmp(i)
                    eigtmp(i)=p
                    call sswap(nbf,temp(1,i),1,temp(1,k),1)
                endif
  400        continue
         endif
  100 continue
      if (count.ne.nbf) then
          call lnkerr('error in orbital count in knvec')
      endif
      call scopy(count*nbf,temp,1,c,1)
      call scopy(count,eigtmp,1,eigval,1)
      call icopy(tmpsym,orbsym,count)
c
c----------------------------------------------------------------------c
c                 write out symmetry information                       c
c----------------------------------------------------------------------c
      write (iout,1)
      call iosys ('write integer "kohn symmetries" to rwf',1,
     1             nirrep,0,' ')
      call izero(tmpsym(1,2),nirrep)
      do 600 irrep=1,nirrep
         call izero(tmpsym(1,1),count)
         ii=0
         do 700 i=1,count
            if (orbsym(i).eq.irrep) then
                ii=ii+1
                tmpsym(ii,1)=i
            endif
  700    continue
         if (ii.gt.0) then
             call iosys ('write integer "symmos-'//itoc(irrep)//'" to'//
     1                    ' rwf',ii,tmpsym(1,i),0,' ')
             write (iout,2) irrep, (tmpsym(i,1),i=1,ii)
         endif
         tmpsym(irrep,2)=ii
         call iosys ('write integer "kohn mo numsym" to rwf',
     1                nirrep,tmpsym(1,2),0,' ')
  600 continue
      write (iout,3)
      return
c
c
    1 format (//,10x,'molecular orbital symmetry information')
    2 format (//,5x,'irreducible representation',1x,i1,1x,
     1         /,5x,'orbitals',( (/,5x,10(i4,1x))))
    3 format(///)
      end
