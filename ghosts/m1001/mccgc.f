*deck %W%  %G%
      subroutine mccgc(ni,ncork,nmk,cg,
     $     nf35,r,g,tfile,bufix,lbufso,rabint,incor)
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
      implicit real*8(a-h,o-z)
cc
      character*(*) tfile,file*16
cc
      dimension g(2)
      dimension cg(nmk,ncork)
      real*8 bufix(lbufso),rabint(*)
      dimension r(nmk,ncork)
      common /pcpack/ ipkt, nhext, ipkes
c
      common /io/ inp,iout
c
c-----------------------------------------------------------------------
c
c --- description     this routine makes contributions from a block
c                     of core coulomb or exchange integrals and the
c                     corresponding density matrix to c*g.
c                     the integrals (kl ij) are stored in alchemy
c                     transformed integral order.
c
c --- input
c
c     ni              no of i active orbitals.
c     ncork           no of k core orbitals.
c     nmk             no of k orbitals, active and virtual.
c     g(--)           density matrix elements.
c     nf35            fortran no for dataset containing transformed
c                     integrals.
c
c --- working storage
c
c     r(nmk,ncork)
c
c --- output
c
c     cg(nmk,ncork)
c
c-----------------------------------------------------------------------
      if (ncork .le. 0) return
      file=tfile
c     do 901 i = 1, ncork
c 901 write (iout,9000) (cg(j,i),j=1,nmk)
c     nij = ni * (ni + 1) / 2
c     write (iout,9002) (g(ii),ii=1,nij)
c9002 format(' *mccgc g '//4(1x,f16.8))
c
      ij = 0
      nr = ncork * nmk
c
      if(incor.eq.0) then
         call iosys('read real '//file//' from rwf without rewinding',
     $        lbufso,bufix,0,' ')
      end if
c
      ix=0
      nni=ni*(ni+1)/2
c
      do 100 i = 1, ni
cc
         do 200 j = 1, i
            ij=ij+1
            if(incor.eq.0) then
               call scopy(nr,bufix(ix+1),1,r,1)
               ix=ix+nr
               if (ix+nr.gt.lbufso) then
                  if (ij.lt.nni) then
                     call iosys('read real '//file//' from rwf '//
     $                    'without rewinding',lbufso,bufix,0,' ')
                     ix=0
                  end if
               end if
            else
               call scopy(nr,rabint(ix+1),1,r,1)
               ix=ix+nr
            end if
cc
cc
            t = g(ij)
ccc
            do 400 k = 1, ncork
cccc
               do 550 m = 1, nmk
                  cg(m,k) = cg(m,k) + t * r(m,k)
 550           continue
cccc
 400        continue
ccc
 200     continue
cc
 100  continue
c
c     do 900 i = 1, ncork
c 900 write (iout,9000) (cg(j,i),j=1,nmk)
c9000 format(' *mccgc cg '//4(1x,f16.8))
c
      return
c6000 write (iout,9001)
c9001 format(//'0****** mccgc',6x,' error reading ordered integrals')
c     call lnkerr(' ')
      end
