*deck @(#)dismat.f	5.1  11/6/94
      subroutine dismat(natoms,c,dis,conver)
c***begin prologue     dismat
c***date written       850601  yymmdd
c***revision date      870320  yymmdd
c     march 20, 1987   rlm at lanl
c     crowd testing and printing removed.
c***keywords           distance matrix, coordinates
c***author             martin, richard (lanl)
c***source             @(#)dismat.f	5.1   11/6/94
c***purpose            computes the distance matrix.
c***description
c     call dismat(natoms,c,dis)
c       natoms ... number of atoms.
c       c      ... coordinate array, (3,natoms).
c       dis    ... distance matrix (natoms,natoms).
c       conver ... conversion factor to apply to the matrix.
c
c***references
c***routines called    sdot(math),vsub(math),lnkerr(mdutil)
c***end prologue       dismat
      implicit integer(a-z)
      real*8 c(3,natoms),dis(natoms,natoms),conver
      real*8 dissq,sdot,temp(3)
c
c
      call rzero(dis,natoms*natoms)
      do 20 i=1,natoms
         do 10 j=1,i-1
            call vsub(temp,c(1,i),c(1,j),3)
            dissq=sdot(3,temp,1,temp,1)
            dis(i,j)=conver*sqrt(dissq)
            dis(j,i)=dis(i,j)
   10    continue
   20 continue
c
c
      return
      end
