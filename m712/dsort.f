*deck @(#)dsort.f	5.1  11/6/94
      subroutine dsort(d,dsortd,mgrpno,noprim,start,nocart,dpt,
     #                 nnprim,lend,nmgrp,nat,nbtype,atptd,atnod,
     #                 ndmat)
c
c***begin prologue
c***date written       yymmdd   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             saxe, paul    (lanl)
c***source             @(#)dsort.f	5.1   11/6/94
c***purpose
c***description
c***references         (none)
c
c***routines called    (none)
c
c***end prologue
c
      implicit integer (a-z)
c
      real*8 d(nnprim,ndmat),dsortd(lend)
      integer mgrpno(nat,nbtype),noprim(nat,nbtype)
      integer start(nat,nbtype),nocart(nbtype)
      integer dpt(nmgrp,nmgrp),atptd(nat,nat),atnod(nat,nat)
c
c     ----- now loop through all momentum groups, forming their density-
c           submatrices as we go
c
      pt=0
      do 400 iatom=1,nat
         do 300 jatom=1,nat
            atptd(iatom,jatom)=pt
            do 200 itype=1,nbtype
               imgrp=mgrpno(iatom,itype)
               if (imgrp.le.0) go to 200
               nprimi=noprim(iatom,itype)
               istart=start(iatom,itype)
               nfi=nocart(itype)
               do 100 jtype=1,nbtype
                  jmgrp=mgrpno(jatom,jtype)
                  if (jmgrp.le.0) go to 100
                  nprimj=noprim(jatom,jtype)
                  jstart=start(jatom,jtype)
                  nfj=nocart(jtype)
                  dpt(imgrp,jmgrp)=pt
                  do 50 dmat=1,ndmat
                     do 40 jf=1,nfj
                        jpos=jstart+jf-nfj
                        do 30 if=1,nfi
                           ipos=istart+if-nfi
                           do 20 jprim=1,nprimj
                              jj=jpos+jprim*nfj
                              do 10 iprim=1,nprimi
                                 ii=ipos+iprim*nfi
                                 if (ii.ge.jj) then
                                    ij=ii*(ii-1)/2+jj
                                 else
                                    ij=jj*(jj-1)/2+ii
                                 end if
                                 pt=pt+1
                                 dsortd(pt)=d(ij,dmat)
   10                         continue
   20                      continue
   30                   continue
   40                continue
   50             continue
  100          continue
  200       continue
            atnod(iatom,jatom)=pt-atptd(iatom,jatom)
  300    continue
  400 continue
c
c
      return
      end
