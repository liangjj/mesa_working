*deck @(#)transd.f	5.1  11/6/94
      subroutine transd(prmint,conint,nprimi,nprimj,nconti,ncontj,
     #                  itype,jtype,a,b,t1,len1,lenblk,imin,imax,
     #                  jmin,jmax,nocart,iatom,jatom,d,nnp,nat,start,
     #                  nbtype,nobf,grad)
c
      implicit integer (a-z)
c
      real*8 prmint(nprimi,nprimj,lenblk,4),conint(nconti,ncontj,lenblk)
      real*8 a(nprimi,nconti,imin:imax),b(nprimj,ncontj,jmin:jmax)
      real*8 t1(len1),d(nnp),grad(3,nat)
      integer nocart(nbtype),start(nat,nbtype),nobf(nbtype)
c
c     ----- transform the primitive derivative integrals, and
c            then sum into the correct place in the final array
c
      do 100 coord=1,3
c
         call trans1(prmint(1,1,1,coord+1),conint,nprimi,nprimj,
     #               nconti,ncontj,a,b,t1,len1,lenblk,imin,imax,jmin,
     #               jmax,nocart)
c
c        ----- put this angular momentum block in the right places
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
c
                     grad(coord,iatom)=grad(coord,iatom)+
     #                      conint(ic,jc,intgrl)*d(ij)
                     grad(coord,jatom)=grad(coord,jatom)-
     #                      conint(ic,jc,intgrl)*d(ij)
    1             continue
    2          continue
    3       continue
    4    continue
c
  100 continue
c
c
      return
      end