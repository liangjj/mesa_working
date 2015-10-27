*deck @(#)mcgcor.f	5.1  11/6/94
      subroutine mcgcor(ni,ncork,nmk,lokk,lenk,mix,cm,
     $     nf35,r,g,tfile,buf,lbufso,rabint,incor)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcgcor.f	5.1   11/6/94
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
      character*(*) tfile,file*16
cc
      real*8 buf(lbufso),rabint(*)
      dimension lokk(2),lenk(2),mix(2),cm(nmk,ncork)
      dimension g(2)
      dimension r(nmk,ncork)
      common /pcpack/ ipkt, nhext, ipkes
c
      common /io/ inp,iout
c
c-----------------------------------------------------------------------
c
c --- description     this routine makes contributions from a block
c                     of core coulomb or exchange integrals to g(ij).
c                     the integrals (kl ij) are stored in alchemy
c                     transformed integral order.
c
c --- input
c
c     ni              no of i active orbitals.
c     ncork           no of k core orbitals.
c     nmk             no of k orbitals, active and virtual.
c     lokk(nk)        lokk(k) starting postion of the kth vector in
c                      the arrays mix and cm.
c     lenk(nk)        length of the kth vector.
c     mix(--)         indices of vector components.
c     cm(--)          vector components.
c     nf35            fortran no for dataset containing transformed
c                     integrals.
c
c --- working storage
c
c     r(nmk,ncork)
c
c --- output
c
c     g(2)
c
c-----------------------------------------------------------------------
      if (ncork .le. 0) return
c
      file=tfile
c
      ij = 0
      ix=0
      mij=0
      nni=ni*(ni+1)/2
c
      if(incor.eq.0) then
         call iosys('read real '//file//' from rwf without rewinding',
     $        lbufso,buf,0,' ')
      end if
c
      nr = ncork * nmk
c
      do 100 i = 1, ni
cc
         do 200 j = 1, i
c
            if(incor.eq.0) then
               call vmove(r,buf(ix+1),nr)
               mij=mij+1
               ix=ix+nr
               if (ix+nr.gt.lbufso) then
                  if (mij.lt.nni) then
                     call iosys('read real '//file//' from rwf '//
     $                    'without rewinding',lbufso,buf,0,' ')
                     ix=0
                  end if
               end if
            else
               call vmove(r,rabint(ix+1),nr)
               ix=ix+nr
            end if
c
c       write(iout,*)' mcgcor i j ',i,j
c       write(iout,19000) r
19000       format(5(1x,f12.8))
ccc
cc      call readpk(nf35,r,r,nr,nhext,ierr,iend)
cc      if (ierr .ne. 0 .or. iend .ne. 0) go to 6000
ccc
            t = 0.d0
            do 400 k = 1, ncork
               mk1 = lokk(k) + 1
               mk2 = mk1 + lenk(k) - 1
               if (mk2 .lt. mk1) go to 400
cccc
               do 550 m = mk1, mk2
                  t = t + r(mix(m),k) * cm(mix(m),k)
 550           continue
cccc
 400        continue
ccc
            g(ij+1) = g(ij+1) + t
            ij = ij + 1
 200     continue
cc
 100  continue
c
c     nij = ni * (ni + 1) / 2
c     write (iout,1002) (g(ij),ij=1,nij)
c1002 format(' *mcgcor g '//4(1x,f16.8))
      return
 6000 write (iout,9001)
 9001 format(//'0****** mcgcor',6x,' error reading ordered integrals')
      call lnkerr(' ')
      end
