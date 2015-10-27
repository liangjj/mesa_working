*deck @(#)dirac.f	1.1  11/30/90
      subroutine dirac(c,ex,cont,ptprim,noprim,nocont,ptcont,
     #                  nat,nprim,ntypes,nbtype,nnp,ncont,
     #                  start,nbasis,nocart,nobf,maxmom,minmom,mintyp,
     #                  nx,ny,nz,ops,gridx,gridy,gridz,ngridx,ngridy,
     #                  ngridz,phix,phiy,phiz,
     #                  inpf)
c***begin prologue     dirac
c***date written       891219   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           dirac, link 1990, amplitude, integrals
c***author             martin, richard (lanl)
c***source             @(#)dirac.f	1.1   11/30/90
c***purpose            evaluates the amplitude of the primitive cartesian
c                      gaussian functions at a specific point in space.
c***description
c       dirac evaluates three components for every primitive function.
c       these are returned in phix,phiy,phiz.  the total amplitude at some
c       point x,y,z is given by their product. these values do not include
c       normalization or the contraction coefficients.  those contributions
c       are implicitly folded into the transformation matrix constructed in m102
c       and stored in the iosys array t(prim,cont).
c
c***references
c
c***routines called
c
c***end prologue       dirac
c
      implicit integer (a-z)
c
      character*(*) ops
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 gridx(ngridx),gridy(ngridy),gridz(ngridz)
      real*8 expon
      logical logkey
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 phix(ngridx,nprim),phiy(ngridy,nprim),phiz(ngridz,nprim)
      real*8 rx,ry,rz
c
      common/io/inp,iout
c
c
 1000 format(1x,'primitive function polynomial amplitudes.')
 1001 format(5x,'x-amplitudes:')
 1002 format(5x,'y-amplitudes:')
 1003 format(5x,'z-amplitudes:')
 1010 format(/5x,'x:',f10.6)
 1011 format(/5x,'y:',f10.6)
 1012 format(/5x,'z:',f10.6)
 1020 format(5x,8e15.6)
c
c
c     evaluate the primitive amplitudes on the grid.
      call rzero(phix,ngridx*inpf)
      call rzero(phiy,ngridy*inpf)
      call rzero(phiz,ngridz*inpf)
c     loop over all atoms.
c     the loop structure is chosen to agree with the convention in m102.
      pf=0
      do 9 iatom=1,nat
c           loop over the types of shells on this atom; s,p,d,sp,etc.
c           disregard the ecp information in types 17-25.
         do 8 itype=1,nbtype
c           are there any primitive functions on this center of this type?
            if (noprim(iatom,itype).gt.0) then
c              loop over all the contractions in this shell.
c              do 7 ncontr=1,nocont(iatom,itype)
c                 loop over all the primitives in this shell.
                  do 6 prim=1,noprim(iatom,itype)
c                    loop over all the angular momentum components,
c                    e.g. s and p for an sp shell.
                     ptpf=ptprim(iatom,itype)+prim-1
                     expon=ex(ptpf)
                     do 5 angmom=minmom(itype),maxmom(itype)
c                       loop over all the individual cartesian components
c                       of this type.
                        do 4 cart=1,nocart(angmom)
                           ptpow=mintyp(itype)-1
                           pf=pf+1
                           n=nx(ptpow+cart)
                           l=ny(ptpow+cart)
                           m=nz(ptpow+cart)
c                          loop over the points on the grid.
c                          don't worry about the normalization;
c                            i.e. the factor of sqrt(3) for dx2,dy2,dz2,etc.
c                            that is handled in the contraction coefficients.
                           do 3 pt=1,ngridx
                              rx=gridx(pt)-c(1,iatom)
                              phix(pt,pf)=exp(-expon*rx*rx)
                              if(n.ne.0) then
                                 phix(pt,pf)=phix(pt,pf)*(rx**n)
                              endif
    3                      continue
                           do 2 pt=1,ngridy
                              ry=gridy(pt)-c(2,iatom)
                              phiy(pt,pf)=exp(-expon*ry*ry)
                              if(l.ne.0) then
                                 phiy(pt,pf)=phiy(pt,pf)*(ry**l)
                              endif
    2                      continue
                           do 1 pt=1,ngridz
                              rz=gridz(pt)-c(3,iatom)
                              phiz(pt,pf)=exp(-expon*rz*rz)
                              if(m.ne.0) then
                                 phiz(pt,pf)=phiz(pt,pf)*(rz**m)
                              endif
    1                      continue
    4                   continue
    5                continue
    6             continue
    7          continue
         endif
    8    continue
    9 continue
c
c     check the number of primitive functions processed in the loop.
      if(inpf.ne.pf) call lnkerr('npf wrong in dirac')
c
c     print the integrals.
      if(logkey(ops,'print=graphics=amplitudes',.false.,' ')) then
         write(iout,1000)
         write(iout,1001)
         do 30 pt=1,ngridx
            write(iout,1010) gridx(pt)
            write(iout,1020) (phix(pt,pf),pf=1,inpf)
   30    continue
         write(iout,1002)
         do 35 pt=1,ngridy
            write(iout,1011) gridy(pt)
            write(iout,1020) (phiy(pt,pf),pf=1,inpf)
   35    continue
         write(iout,1003)
         do 40 pt=1,ngridz
            write(iout,1012) gridz(pt)
            write(iout,1020) (phiz(pt,pf),pf=1,inpf)
   40    continue
      endif
c
c
      return
      end
