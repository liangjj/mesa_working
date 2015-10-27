*deck @(#)ehuckl.f	5.1  11/6/94
      subroutine ehuckl(start,nprim,ncont,natoms,nbtype,nbasis,
     $                  ncart,minmom,maxmom,atomno,ediag,corflg)
c***begin prologue     ehuckl
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           huckel, initial guess
c***author             martin, richard (lanl)
c***source             @(#)ehuckl.f	5.1   11/6/94
c***purpose            sets up diagonal huckel energies and core flags.
c***description
c     call ehuckl(start,nprim,ncont,natoms,nbtype,nbasis,
c                 ncart,minmom,maxmom,atomno,ediag,corflg)
c       start   an array(natoms,nbtype) pointing to the index(-1) of the
c               first function of a given type on a specific atom.
c       nprim   the number of gaussian primitives.
c       ncont   the number of gaussian contractions.
c       natoms  the number of atoms.
c       nbtype  the number of basis function types.
c       nbasis  the number of basis functions.
c       ncart   a vector(nbtype) containing the number of cartesian
c               components in each basis function type.
c       minmom  a vector(nbtype) containing the minimum angular momentum
c               in each basis function type.
c       maxmom  a vector(nbtype) containing the maximum angular momentum
c               in each basis function type.
c       atomno  a vector(natoms) containing the atomic numbers.
c       ediag   a vector(nbasis) containing the huckel diagonal elements.
c       corflg  a vector(nbasis) containing the orbital type for
c               each element of the bas.  a negative orbital type denotes
c               a core function.
c***references         (none)
c***routines called    eneg(m401)
c***end prologue       ehuckl
      implicit integer(a-z)
      integer start(natoms,nbtype),nprim(natoms,nbtype)
      integer ncont(natoms,nbtype),ncart(nbtype)
      integer minmom(nbtype),maxmom(nbtype),atomno(natoms)
      integer corflg(nbasis)
      integer shlpnt(7,4),maxpqn(110)
      real*8 ediag(nbasis),eneg
c
c     shlpnt maps a principal quantum number and angular momentum
c     into an orbital type.
      data shlpnt/1, 2, 4, 6, 9,12,16,
     $            0, 3, 5, 7,10,13,17,
     $            0, 0, 8,11,14,18, 0,
     $            0, 0, 0,15,19, 0, 0/
      data maxpqn/2*1, 8*2, 8*3, 18*4, 18*5, 32*6, 24*7/
      save shlpnt,maxpqn
c
c
c     start timing.
c
c
      do 100 atom=1,natoms
         do 90 type=1,nbtype
            if(nprim(atom,type).le.0) goto 90
            nbf=start(atom,type)
            do 80 func=1,ncont(atom,type)
               do 75 angmom=minmom(type),maxmom(type)
c                 the next confusing little bit of code is necessary
c                 because the minimum basis which was set up may have had
c                 core orbitals deleted if an ecp is in use.
                  pqn=maxpqn(atomno(atom))
c                 account for the fact that d,f shells lag in
c                 principal quantum number during aufbau.
                  if(angmom.gt.1) then
                     pqn=pqn-angmom+1
                  endif
                  pqn=pqn-ncont(atom,type)+func
                  j=shlpnt(pqn,angmom+1)
                  if(func.lt.ncont(atom,type)) j=-j
                  do 70 k=1,ncart(angmom+1)
                     nbf=nbf+1
                     corflg(nbf)=j
                     ediag(nbf)=eneg(atomno(atom),j)
   70             continue
   75          continue
   80       continue
   90    continue
  100 continue
c
c     stop timing.
c
c
      return
      end
