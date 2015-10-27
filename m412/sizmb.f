*deck @(#)sizmb.f	5.1  11/6/94
      subroutine sizmb(natoms,nprim,ncont,nbas,atomno,atomz)
c***begin prologue     sizmb
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           minimum basis, basis
c***author             martin, richard (lanl)
c***source             @(#)sizmb.f	5.1   11/6/94
c***purpose            returns basis set dimensions appropriate for a given
c                      minimum basis.
c***description
c     call sizmb(natoms,nprim,ncont,nbas,atomno,atomz)
c       natoms  the number of atoms.
c       nprim   the number of primitives.
c       ncont   the number of contraction coefficients.
c       nbas    the number of basis functions.
c       atomno  the atomic number vector(natoms).
c       atomz   the nuclear charge vector(natoms).
c     given the atomic numbers of the centers in the molecule, this
c     routine returns the dimensions necessary to store information about a
c     minimum basis on the molecule.  note that those orbitals replaced by
c     effective core potential are omitted from the basis (see minfrz).
c***references         (none)
c***routines called    (none)
c***end prologue       sizmb
c
      implicit integer(a-z)
c
      real*8 atomz(natoms)
c
      parameter (pt1=0.1d00,lensto=3)
c
      dimension atomno(natoms)
      dimension minan(19),maxval(19),minfrz(19),nfunc(19)
      dimension type(19), ntype(4)
c
      data minan/0,2,2,10,10,18,18,20,36,36,38,54,54,56,56,86,86,88,88/
      data maxval/2,10,10,18,18,36,36,54,54,54,86,86,86,86,
     $            120,120,120,120,120/
      data minfrz/0,2,2,10,10,28,28,18,46,46,36,78,78,68,54,110,110,
     $            100,86/
      data nfunc/1,1,3,1,3,1,3,6,1,3,6,1,3,6,10,1,3,6,10/
      data type/1,1,2,1,2,1,2,3,1,2,3,1,2,3,4,1,2,3,4/
      save minan,maxval,minfrz,nfunc,type
c
c
      nshell=0
      nbas=0
      ncont=0
      do 20 i=1,natoms
         call izero(ntype,4)
         if(atomno(i).le.0) goto 20
         nfroz=int(float(atomno(i))-atomz(i)+pt1)
         do 10 j=1,19
            if(atomno(i).gt.minan(j).and.nfroz.le.minfrz(j)) then
               nshell=nshell+1
               nbas=nbas+nfunc(j)
               ntype(type(j))=ntype(type(j))+1
            endif
   10    continue
         do 15 j=1,4
            ncont=ncont+ntype(j)*ntype(j)*lensto
   15    continue
   20 continue
c
c
      nprim=nshell*lensto
c
c
      return
      end
