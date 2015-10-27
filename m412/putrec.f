*deck @(#)putrec.f	5.1  11/6/94
      subroutine putrec(s,conint,strti,strtj,iatom,jatom,itype,jtype,
     #                  nconti,ncontj,nbasi,nbasj,lenblk,nati,natj,
     #                  nbtype,nobf)
c***begin prologue     putrec
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           integrals, transfer
c***author             saxe, paul (lanl)
c***source             @(#)putrec.f	5.1   11/6/94
c***purpose            transfers an angular momentum block of 1-e integrals
c                      from the block array to the full matrix.
c***description
c     call putrec(s,conint,strti,strtj,iatom,jatom,itype,jtype,
c                 nconti,ncontj,nbasi,nbasj,lenblk,naati,naatj,
c                 nbtype,nobf)
c***references         (none)
c***routines called    traktm(mdutil)
c***end prologue       putrec
      implicit integer (a-z)
c
      real*8 s(nbasi,nbasj),conint(nconti,ncontj,lenblk)
      integer strti(nati,nbtype),strtj(natj,nbtype),nobf(nbtype)
c
c     ----- start timing -----
c
c
      numi=nobf(itype)
      numj=nobf(jtype)
      istart=strti(iatom,itype)
      jstart=strtj(jatom,jtype)
      intgrl=0
c
      do 4 if=1,numi
         do 3 jf=1,numj
            intgrl=intgrl+1
            do 2 jc=1,ncontj
               j=jstart+(jc-1)*numj+jf
               do 1 ic=1,nconti
                  i=istart+(ic-1)*numi+if
                  s(i,j)=conint(ic,jc,intgrl)
    1          continue
    2       continue
    3    continue
    4 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
