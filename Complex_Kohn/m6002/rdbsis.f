*deck @(#)rdbsis.f	1.1 9/7/91
c***begin prologue     rdbsis
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rdbsis, link 6001, bigmolli
c***author             schneider, barry (lanl)
c***source             m6001
c***purpose            read in cartesian ao basis parameters
c***description        fills up arrays with numerical values of n, l, m
c***                   and center in old polyatom format
c*** 
c
c***references       
c
c***routines called
c***end prologue       rdbsis
      subroutine rdbsis(ntot,mxprcn,aosym)
      implicit integer (a-z)
      parameter (dimpr=300)
      real *8 alf, cont, eta
      character *(*) aosym
      character *3 ctype
      dimension aosym(dimpr)
      common /aosi/ npr, ncon, nxyzc(dimpr,4), nprc(dimpr), 
     1              sncon, smllst(dimpr,2)
      common /aosr/ alf(dimpr), cont(dimpr)
      common /ctype/ ctype(0:3,0:3,0:3)
      common /factr/ nfirst(dimpr), nlast(dimpr), ncntr(dimpr), 
     1               nx(dimpr), ny(dimpr), nz(dimpr), eta(dimpr,5)
      call iosys ('read integer "no. primitives" from kohndt',1,
     1             npr,0,' ')
      call iosys ('read integer "no. contracted" from kohndt',1,
     1             ncon,0,' ')
      call iosys ('read integer start from kohndt',ncon,nfirst,
     1             0,' ')
      call iosys ('read integer stop from kohndt',ncon,nlast,
     1             0,' ')
      call iosys ('read integer "x power" from kohndt',npr,nx,0,' ')
      call iosys ('read integer "y power" from kohndt',npr,ny,0,' ')
      call iosys ('read integer "z power" from kohndt',npr,nz,0,' ')
      call iosys ('read real eta from kohndt',1500,eta,0,' ')
      call iosys ('read integer atom from kohndt',npr,ncntr,0,' ')
      mxprcn=0
      ntot=0
      do 10 i=1,ncon
         nprc(i)=nlast(i)-nfirst(i)+1
         mxprcn=max(mxprcn,nprc(i))
         do 20 j=nfirst(i),nlast(i)
            ntot=ntot+1
            nxyzc(ntot,1)=nx(ntot)
            nxyzc(ntot,2)=ny(ntot)
            nxyzc(ntot,3)=nz(ntot)
            nxyzc(ntot,4)=ncntr(ntot)
            alf(ntot)=eta(ntot,4)
            cont(ntot)=eta(ntot,5)
            aosym(ntot)=ctype(nx(ntot),ny(ntot),nz(ntot))
   20    continue
   10 continue    
      return
      end
