*deck @(#)transf.f	5.1  11/6/94
      subroutine transf(prmint,conint,nprimi,nprimj,nconti,ncontj,
     #     itype,jtype,a,b,t1,len1,lenblk,imin,imax,
     #     jmin,jmax,nocart,iatom,jatom,d,nnp,nat,start,
     #     nbtype,nobf,nint,ld2e,nd2e)
c
c***begin prologue     transf
c***date written       861215  (yymmdd)
c***revision date      871113  (yymmdd)
c
c   13 november 1987   pws at lanl
c      changing to contract the contracted integrals directly with the
c      density (lagrangian) to form the second derivative contribution.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)transf.f	5.1   11/6/94
c***purpose            transformation of second-derivative one-electron
c                        integrals from primitive to contracted basis.
c***description
c
c***references
c***routines called
c***end prologue       transf
c
      implicit integer (a-z)
c
      real*8 prmint(nprimi,nprimj,lenblk,nint)
      real*8 conint(nconti,ncontj,lenblk)
      real*8 a(nprimi,nconti,imin:imax),b(nprimj,ncontj,jmin:jmax)
      real*8 t1(len1)
      real*8 d(nnp)
      real*8 ld2e(nd2e)
      real*8 t
      integer nocart(nbtype),start(nat,nbtype),nobf(nbtype)
c
c     ----- transform the primitive second-derivative integrals, and
c            then sum into the correct place in the final array
c
      ijc=0
      do 100 icoord=1,3
         icia=icoord+(iatom-1)*3
         icja=icoord+(jatom-1)*3
         do 90 jcoord=1,icoord
            jcia=jcoord+(iatom-1)*3
            jcja=jcoord+(jatom-1)*3
            iiji=icia*(icia-1)/2+jcia
            iijj=icia*(icia-1)/2+jcja
            ijjj=icja*(icja-1)/2+jcja
            jiii=jcia*(jcia-1)/2+icia
            jiij=jcia*(jcia-1)/2+icja
            jjij=jcja*(jcja-1)/2+icja
            ijc=ijc+1
c
            call trans1(prmint(1,1,1,ijc+4),conint,nprimi,nprimj,
     #           nconti,ncontj,a,b,t1,len1,lenblk,imin,imax,jmin,
     #           jmax,nocart)
c
c           ----- put this angular momentum block in the right places
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
                        t=conint(ic,jc,intgrl)*d(ij)
                        ld2e(iiji)=ld2e(iiji)+t
                        ld2e(iijj)=ld2e(iijj)-t
                        ld2e(ijjj)=ld2e(ijjj)+t
                        if (icoord.ne.jcoord) then
                          ld2e(jiij)=ld2e(jiij)-t
                        end if
    1                continue
    2             continue
    3          continue
    4       continue
c
   90    continue
  100 continue
c
c
      return
      end
