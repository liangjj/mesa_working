*deck  %W% %G%
      subroutine form3(prmint,xyz,npint,lenblk,imax,jmax,kmax,
     $                  mini,maxi,minj,maxj,mink,maxk,nx,ny,nz)
c***begin prologue     form3
c***date written       850601  (yymmdd)
c***revision date      920417  (yymmdd)  .
c      april 17,1992   rlm at lanl
c         modified m302/fmonel to handle 3-center integrals
c***keywords           one-electron, three-center, integrals
c***author             saxe, paul and martin, richard (lanl).
c***source
c***purpose            assembles the two-dimensional 3-center integrals into 
c                      into primitive 3-center one-electron integrals.
c***description
c                      call form3(prmint,xyz,npint,lenblk,imax,jmax,kmax,
c                                  mini,maxi,minj,maxj,mink,maxk,nx,ny,nz)
c
c                        prmint   sink for primitive integrals.
c                        xyz      source of two-dimensional integrals.
c                        npint    number of primitive integrals to form.
c                                    the product of nprimi*nprimj*nprimk
c                        lenblk   length of the angular momentum shell product.
c                        imax     maximum l value of shell block i.
c                        jmax     maximum l value of shell block j.
c                        kmax     maximum l value of shell block k.
c                        mini     pointers into the polynomial powers for 
c                        maxi        shell block i.
c                        minj     pointers into the polynomial powers for 
c                        maxj        shell block j.
c                        mink     pointers into the polynomial powers for 
c                        maxk        shell block k.
c                        nx       polynomial powers
c                        ny
c                        nz
c
c***references
c***routines called    vmul(math)
c***end prologue       form3
      implicit integer (a-z)
c
c     input arrays(ummodified)
      real*8 xyz(npint,0:imax,0:jmax,0:kmax,3)
      integer nx(*),ny(*),nz(*)
c
c     output arrays
      real*8 prmint(npint,lenblk)
c
c     local variables
      logical debug
c
      parameter (debug=.false.)
c
      common/io/inp,iout
c
c
      intgrl=0
c
c     vary the "auxiliary" center most slowly.
      do 30 k=mink,maxk
         kx=nx(k)
         ky=ny(k)
         kz=nz(k)
         do 20 i=mini,maxi
            ix=nx(i)
            iy=ny(i)
            iz=nz(i)
            do 10 j=minj,maxj
               jx=nx(j)
               jy=ny(j)
               jz=nz(j)
               intgrl=intgrl+1
c
               call vmul(prmint(1,intgrl),xyz(1,ix,jx,kx,1),
     #                                    xyz(1,iy,jy,ky,2),npint)
               call vmul(prmint(1,intgrl),prmint(1,intgrl),
     #                                    xyz(1,iz,jz,kz,3),npint)
   10       continue
   20    continue
   30 continue
      if(debug) then
         do 40 i=1,intgrl
            write(iout,*) 'intgrl:',i
            write(iout,*) (prmint(j,i),j=1,npint)
   40    continue
      endif
c
c
      return
      end
