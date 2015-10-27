*deck @(#)tagcor.f	1.3 8/8/91 
      subroutine tagcor(nao,nco,nob,ta,
     1                  r,g,tfile,buf,lbufso,rabint,incor,ndf)
cc
      implicit integer(a-z)
      character*(*) tfile,file*16
cc
      real*8 buf(lbufso),rabint(*)
      real*8 ta(nob,nco)
      real*8 g(2)
      real*8 r(nob,nco)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      common /io/ inp,iout
c
c-----------------------------------------------------------------------
c
c  vectorized version by b.h.lengsfield 16may1986
c
c    g(ij)=g(ij)+rabint(m,n,ij)*ta(m,n)
c
c
c --- description     this routine makes contributions from a block
c                     of core coulomb or exchange integrals to g(ij).
c                     the integrals (kl ij) are stored in alchemy
c                     transformed integral order.
c
c --- input
c
c     nao              no of i active orbitals.
c     nco           no of k core orbitals.
c     nob             no of k orbitals, active and virtual.
c     lokk(nk)        lokk(k) starting postion of the kth vector in
c                      the arrays mix and cm.
c     lenk(nk)        length of the kth vector.
c     mix(--)         indices of vector components.
c     ta(--)          vector of der. overlap components.
c     nf35            fortran no for dataset containaong transformed
c                     integrals.
c
c --- working storage
c
c     rabint(nob,nco,nao*(nao+1))
c
c --- output
c
c     g(2)
c
c-----------------------------------------------------------------------
      if (nco .le. 0) return
c
      file=tfile
c
      nnao=nao*(nao+1)/2
      nr=nco*nob
c..bhl
c      write(iout,*)'  tagcor:  ta(nob,nco)  '
c      call matout(ta,nob,nco,nob,nco,iout)
c..bhl
      if(incor.eq.0) then
c
         lpass=lbufso/nr
         lpass=min(lpass,nnao)
         npass=(nnao-1)/lpass+1
c
         if(lpass.lt.1) then
            call lnkerr(' m1001: buffer size too small in tagcor')
          endif
c
         call iosys('rewind '//file//' on rwf',0,0,0,' ')
c
         nni=nnao
         ix=1
c
         do 10 i=1,npass
c
            nm=min(nni,lpass)
            nni=nni-nm
            lread=nm*nr
c
            call iosys('read real '//file//' from rwf '
     #                 //'without rewinding',lread,buf,0,' ')
c
c..cos      call mxmb(buf,nr,1,ta,1,nr,g(ix),1,nnao,nm,nr,ndf)
c
c           call sgmm(nm,nr,ndf,buf,nr,ta,nr,g(ix),nnao,4,2)
            call sgemm('t','n',nm,ndf,nr,one,buf,nr,ta,nr,
     $                  one,g(ix),nnao)
c
            ix=ix+nm
c
  10     continue
c
      else
c
c..bhl
c        mx=1
c        do 100 jm=1,nnao
c           write(iout,*)' r(nob,nco,jm)  jm= ',jm
c           call matout(rabint(mx),nob,nco,nob,nco,iout)
c           mx=mx+nob*nco
c 100    continue
c
c        write(iout,*)'  g  before apbtc '
c        write(iout,11)(g(im),im=1,nnao)
c..bhl
         call apbtc(g,rabint,ta,nnao,nr,ndf)
c..bhl
c         write(iout,*)'  g  after  apbtc '
c         write(iout,11)(g(im),im=1,nnao)
c..bhl
      end if
c
  11  format(5(1x,f12.8))
c
      return
      end
