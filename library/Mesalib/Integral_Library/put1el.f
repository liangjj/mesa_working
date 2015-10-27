*deck @(#)put1el.f	5.1  11/6/94
      subroutine put1el(s,conint,start,iatom,jatom,itype,jtype,nconti,
     #                  ncontj,nnp,lenblk,nat,nbtype,nobf)
c
c***module to transfer an angular-momentum block of one-electron
c   integrals from conint to the array s.
c
c paul saxe                 23 july 1984                     lanl
c
      implicit integer (a-z)
c
      real*8 s(nnp),conint(nconti,ncontj,lenblk)
      integer start(nat,nbtype),nobf(nbtype)
c
c     ----- start timing -----
c
c
      numi=nobf(itype)
      numj=nobf(jtype)
      istart=start(iatom,itype)
      jstart=start(jatom,jtype)
      intgrl=0
c
      do 4 if=1,numi
         do 3 jf=1,numj
            intgrl=intgrl+1
            do 2 jc=1,ncontj
               j=jstart+(jc-1)*numj+jf
               do 1 ic=1,nconti
                  i=istart+(ic-1)*numi+if
                  if (i.ge.j) then
                     ij=i*(i-1)/2+j
                  else
                     ij=j*(j-1)/2+i
                  end if
                  s(ij)=conint(ic,jc,intgrl)
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
