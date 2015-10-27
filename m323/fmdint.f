*deck @(#)fmdint.f	5.1  11/6/94
      subroutine fmdint(i2,lenblk,rv,angmom,imax,jmax,mmax,
     #                  lmax,nroots,px,py,pz,lenxyz,mintyp,nocart,
     #                   nbtype,nv,di,ndcen,leni,
     #                   pint,lnpint,ic,npi,nci,jc,npj,ncj,
     #                   kc,npk,nck,lc,npl,ncl,index,len,cdint,
     #                   ncint)
c
      implicit integer (a-z)
c
      real*8 i2(leni,3),di(leni,3,ndcen)
      real*8 pint(lnpint),ic(npi,nci),jc(npj,ncj),kc(npk,nck)
      real*8 lc(npl,ncl),cdint(lenblk,ncint,3,ndcen)
      integer index(len,6)
      integer angmom(4),px(lenxyz),py(lenxyz),pz(lenxyz)
      integer mintyp(nbtype),nocart(nbtype)
c
c     ----- and the gradient contributions -----
c
      do 1 dcen=1,ndcen
c
         call fmint(di(1,1,dcen),i2(1,2),i2(1,3),
     $               lenblk,rv,angmom,imax,jmax,
     #               mmax,lmax,nroots,px,py,pz,lenxyz,mintyp,nocart,
     #               nbtype,nv,
     #                pint,lnpint,ic,npi,nci,jc,npj,ncj,kc,npk,nck,
     #                lc,npl,ncl,index,len,cdint(1,1,1,dcen))
c
         call fmint(i2(1,1),di(1,2,dcen),i2(1,3),
     $               lenblk,rv,angmom,imax,jmax,
     #               mmax,lmax,nroots,px,py,pz,lenxyz,mintyp,nocart,
     #               nbtype,nv,
     #                pint,lnpint,ic,npi,nci,jc,npj,ncj,kc,npk,nck,
     #                lc,npl,ncl,index,len,cdint(1,1,2,dcen))
c
         call fmint(i2(1,1),i2(1,2),di(1,3,dcen),
     $               lenblk,rv,angmom,imax,jmax,
     #               mmax,lmax,nroots,px,py,pz,lenxyz,mintyp,nocart,
     #               nbtype,nv,
     #                pint,lnpint,ic,npi,nci,jc,npj,ncj,kc,npk,nck,
     #                lc,npl,ncl,index,len,cdint(1,1,3,dcen))
c
    1 continue
c
      return
      end
