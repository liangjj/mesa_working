*deck @(#)fmstrt.f	5.1  11/6/94
      subroutine fmstrt(start,nprim,ncont,natoms,nbtype,nobf,nbf,ncbf,
     $                  ncart,npf,pstart)
c***begin prologue     fmstrt
c***date written       850601  (yymmdd)
c***revision date      851015  (yymmdd)  .
c***keywords           one-electron, integrals
c***author             saxe, paul (lanl).
c***source
c***purpose            forms the bookkeeping arrays used to compute integrals.
c***description
c                      call fmstrt(start,nprim,ncont,natoms,nbtype,nobf,nbf,
c                                  ncbf,ncart,npf,pstart)
c***revisions
c        851015 pws  added number of primitive functions (npf) and the
c                    starting location of a primitive shell (pstart)
c
c***references
c***routines called    (none)
c***end prologue       fmstrt
      implicit integer (a-z)
c
      integer start(natoms,nbtype),nprim(natoms,nbtype)
      integer ncont(natoms,nbtype),nobf(nbtype),ncart(nbtype)
      integer pstart(natoms,nbtype)
c
c     ----- start timing -----
c
c
      npf=0
      nbf=0
      ncbf=0
      do 2 atom=1,natoms
         do 1 type=1,nbtype
            if (nprim(atom,type).gt.0) then
               start(atom,type)=nbf
               nbf=nbf+ncont(atom,type)*nobf(type)
               ncbf=ncbf+ncont(atom,type)*ncart(type)
               pstart(atom,type)=npf
               npf=npf+nprim(atom,type)*ncart(type)
            else
               start(atom,type)=-99999
               pstart(atom,type)=-99999
            end if
    1    continue
    2 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
