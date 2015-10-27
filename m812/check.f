*deck @(#)check.f	5.1  11/6/94
      subroutine check(t1,t2,dm,ints,nbf,nnp,ntriang,mcscf)
c
c***begin prologue     check
c***date written       871118   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           testing
c***author             saxe, paul (lanl)
c***source             @(#)check.f	5.1   11/6/94
c
c***purpose            to test that the density matrices traced with
c    integrals give back the energy.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       check
c
      implicit integer (a-z)
c
      character*32 file
c
      logical mcscf
c
      real*8 t1(nbf,nbf)
      real*8 t2(nbf,nbf)
      real*8 dm(nnp,ntriang)
      real*8 ints(nnp,ntriang)
c
      real*8 sdot
      real*8 e1
      real*8 e2
      real*8 enuc
      real*8 e
c
      common /io/ inp,iout
c
c     ----- one-particle section -----
c
      call iosys('read real "kinetic integrals" from rwf',
     $     nnp,ints,0,' ')
      call trtosq(t1,ints,nbf,nnp)
      call iosys('read real "potential integrals" from rwf',
     $     nnp,ints,0,' ')
      call trtosq(t2,ints,nbf,nnp)
      call vadd(t1,t1,t2,nbf**2)
c
      if (mcscf) then
         call iosys('read real "mcscf active ao 1pdm" from rwf',
     $        nnp,dm,0,' ')
      else
         call iosys('read real "ci ao 1pdm" from rwf',nnp,dm,0,' ')
      end if
      call trtosq(t2,dm,nbf,nnp)
c
      e1=sdot(nbf**2,t1,1,t2,1)
c
c     ----- two-particle section -----
c
      if (mcscf) then
         file='mcscf ao 2pdm'
      else
         file='ci ao 2pdm'
      end if
c
      call iosys('rewind "'//file//'" on aoden',0,0,0,' ')
      call iosys('rewind "sorted ao integrals" on ints',0,0,0,' ')
c
      e2=0.0d+00
      maxkl=0
      kl=0
      do 4 k=1,nbf
         do 3 l=1,k
            kl=kl+1
c
c     ----- check that these triangles are in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
               call iosys('read real "sorted ao integrals" from ints '
     #                    //'without rewinding',lnread,ints,0,' ')
               call iosys('read real "'//file//'" from aoden '
     #                    //'without rewinding',lnread,dm,0,' ')
            end if
c
            call trtosq(t1,ints(1,kl-minkl+1),nbf,nnp)
            call trtosq(t2,dm(1,kl-minkl+1),nbf,nnp)
c
            if (k.eq.l) then
               e2=e2+sdot(nbf**2,t1,1,t2,1)
            else
               e2=e2+2.0d+00*sdot(nbf**2,t1,1,t2,1)
            end if
    3    continue
    4 continue
c
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $     1,enuc,0,' ')
c
      e=enuc+e1+e2
c
c
      write (iout,5) e1,e2,enuc,e
 5    format (5x,'checking the energy:',/,
     $     10x,'one-electron:           ',f14.8,/,
     $     10x,'two-electron:           ',f14.8,/,
     $     10x,'nuclear repulsion :     ',f14.8,/,
     $     10x,'total from ao density:  ',f14.8)
c
      call iosys('read real energy from rwf',1,e,0,' ')
c
      write (iout,6) e
 6    format (10x,'total from calculation: ',f14.8)
c
c
      return
      end
