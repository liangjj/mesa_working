*deck @(#)cblock.f	2.1  10/10/91
      subroutine cblock(scfvec,cao,cso,bfsym,bfnum,orbsym,nbf,nblock,nso
     $     ,nsym,ibuff,eigval)
c
c***begin prologue     cblock
c***date written       871125   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           symmetry blocked scf vector
c***author             saxe, paul (lanl)
c***source             @(#)cblock.f	2.1   10/10/91
c
c***purpose            to symmetry block an ao scf vector where
c                      no atoms are symmetry related.
c
c***description        
c
c***references         
c
c***routines called    (none)
c
c***end prologue       cblock
c
      implicit integer (a-z)
c
      integer ibuff(*)
      real*8 cao(nbf,nbf)
      real*8 eigval(nbf)
      real*8 cso(nblock)
      character *8 scfvec
      integer bfsym(nbf)
      integer orbsym(nbf)
      integer bfnum(nbf)
      integer nso(nsym)
c
      common /io/ inp,iout
c
      call qassign(2,scfvec,ibuff,10000b)
      rewind 2
      read(2)
c
c***** assumes the nmo=nao, first loop over nmo
c
      do 130 i=1,nbf
         read(2) k,eigval(i)
         read(2) (cao(j,i),j=1,nbf)
 130  continue
c 
      write (iout,153)
 153  format (/,' the ao scf vector')
      call matout(cao,nbf,nbf,nbf,nbf,iout)
c      
      call putvec(cao,eigval,nbf)
c
c     ----- now to block -----
c
      pt=0
      do 10 sym=1,nsym
         do 9 orb=1,nbf
            if (orbsym(orb).eq.sym) then
               do 8 bf=1,nbf
                  if (bfsym(bf).eq.sym) then
                     pt=pt+1
                     cso(pt)=cao(bf,orb)
                  end if
 8             continue
            end if
 9       continue
 10   continue
c
      call iosys('write real "so scf vector" to rwf',nblock,cso,0,' ')
c
      write (iout,11)
 11   format (/,' the so scf vector')
      call matblk(cso,nso,nsym,iout)
c
c
      return
      end
