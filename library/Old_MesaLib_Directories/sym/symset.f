*deck @(#)symset.f	5.1  11/6/94
      subroutine symset(scset,ns,maxmom,nsymat,ncart,lnsc,relatm,
     $     mcu,momatm,natoms,naords,maxsao)
c
c***begin prologue     scset
c***date written       870716   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           symmetry contraction sets.
c***author             saxe, paul (lanl)
c***source             @(#)symset.f	5.1   11/6/94
c
c***purpose            form pointers to symmetry contraction sets.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       symset
c
      implicit integer (a-z)
c
      integer scset(0:maxmom,ns)
      integer nsymat(ns)
      integer ncart(0:maxmom)
      integer relatm(mcu,ns)
      integer momatm(natoms)
c
c     ----- scset(angular momentum,symmetry related atom set) points
c           to the beginning of the symmetry contraction matrix
c
      naords=0
      maxsao=-1
      pt=0
      do 20 set=1,ns
         atom=relatm(1,set)
         do 10 angmom=0,momatm(atom)
            naords=naords+1
            maxsao=max(maxsao,ncart(angmom)*nsymat(set))
            scset(angmom,set)=pt
            pt=pt+(ncart(angmom)*nsymat(set))**2
 10      continue
         do 15 angmom=momatm(atom)+1,maxmom
            scset(angmom,set)=-9999999
 15      continue
 20   continue
c
      lnsc=pt
c
c
      return
      end
