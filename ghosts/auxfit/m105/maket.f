*deck @(#)maket.f	1.1  11/20/92
      subroutine maket(t,npf,nbf,ptcont,natoms,nbtype,noprim,nocont,
     #                 cf,numcf,nocart,minmom,maxmom,start,pstart)
c
      implicit integer (a-z)
c
      real*8 t(npf,nbf),cf(numcf)
      integer ptcont(natoms,nbtype),noprim(natoms,nbtype)
      integer nocont(natoms,nbtype),nocart(nbtype),minmom(nbtype)
      integer maxmom(nbtype),start(natoms,nbtype),pstart(natoms,nbtype)
c
c     ----- zero the t matrix -----
c
      call rzero(t,npf*nbf)
c
c     ----- loop through atoms and types, filling in the transformation
c
      do 100 atom=1,natoms
         do 90 type=1,nbtype
            if (noprim(atom,type).gt.0) then
               call maket1(t,npf,nbf,noprim(atom,type),
     #                     nocont(atom,type),cf(ptcont(atom,type)),
     #                     minmom(type),maxmom(type),nocart,
     #                     start(atom,type),pstart(atom,type))
            end if
   90    continue
  100 continue
c
c
      return
      end
