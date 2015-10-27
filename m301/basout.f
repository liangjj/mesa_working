*deck @(#)basout.f	1.3  7/31/91
      subroutine basout(coords,ex,cont,ptprim,noprim,nocont,ptcont,
     #           nat,nprim,ntypes,nbtype,nnp,ncont,
     #           start,zan,nocart,nobf,maxmom,mintyp,
     #           nx,ny,nz,minmom,nstart,bstart,blen,bflabl,prnt)
c
c
      implicit integer (a-z)
c
      real*8 ex(nprim),cont(ncont)
      real*8 zan(nat)
      character*16 bflabl(*)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
c
      integer bstart(nat,ntypes),blen(nat,ntypes),nstart(nat)
      logical prnt
c
c     dimensions for the arrays being built for quad codes
c
      real*8 coords(3*nat)
c
c
      common/io/inp,iout
c
      ic=1
c
 
      do 9 iatom=1,nat
         nstart(iatom)=ic
         if (prnt) then
            write(iout,*)' '
            write(iout,*)' atom number ',iatom
            write(iout,*)' charge      ',zan(iatom)
            write(iout,*)' bflabl      ',bflabl(ic)
         endif
         do 7 itype=1,nbtype
            if (noprim(iatom,itype).gt.0) then
               mini=mintyp(itype)
               maxi=mintyp(itype)+nocart(itype)-1
               imax=maxmom(itype)
c
               if(prnt) then
                  write(iout,*)'  primitive type s=1 p=2 d=3 etc. ',
     $                            itype
               endif
c
               nprimi=noprim(iatom,itype)
               nconti=nocont(iatom,itype)
c
               bstart(iatom,itype)=ic
               blen(iatom,itype)=nconti*nocart(itype)
               ic=ic+nconti*nocart(itype)
c
               if(prnt) then
                   write(iout,*)' number of contracted functions',nconti
                   write(iout,*)' number of cartesian functions',
     $                           nocart(itype)
               endif
c
            endif
    7    continue
    9 continue
c
c
      return
      end
