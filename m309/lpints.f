*deck @(#)lpints.f	5.1  11/6/94
      subroutine lpints(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                   nprim,ncont,nat,nbtype,ntypes,maxp2,maxi2,
     $                   i1,j1,c,ex,noprim,ptprim,ptcont,cont,
     $                   maxmom,mini,maxi,minj,maxj,nx,ny,nz,sp,
     $                   ntpse,nlp,zlp,clp,qq)
      implicit real*8(a-h,o-z)
      integer ptprim,ptcont
      logical skip
c
c     controls computation of pseudopotential integrals.
c     details of the method are presented in
c         l.e. mcmurchie and e.r. davidson, j. comp. phys.,44,289(1981).
c         r.l. martin, to be published.
c
      common/center/xa,ya,za,xb,yb,zb,xc,yc,zc
      common/dist/cax,cay,caz,ca,ca2,cbx,cby,cbz,cb,cb2
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/const/zero,one,two,three,four,five,six,ten
      common/prims/igbegn,igend,jgbegn,jgend
      common/limit/istart,jstart,iend,jend
      common/angmax/lamax,lbmax
      dimension ex(nprim),c(3,nat),noprim(nat,ntypes)
      dimension ptprim(nat,ntypes),ptcont(nat,ntypes),cont(ncont)
      dimension maxmom(ntypes)
      dimension nx(*),ny(*),nz(*),nlp(*),zlp(*),clp(*)
      dimension ntpse(*)
      dimension qq(*),sp(*)
c
      parameter (half=5.0d-01)
      call rzero(sp,maxp2*maxi2)
      lamax=imax+1
      lbmax=jmax+1
      istart=mini
      iend=maxi
      jstart=minj
      jend=maxj
      igbegn=i1
      igend=i1+nprimi-1
      jgbegn=j1
      jgend=j1+nprimj-1
      xa=c(1,iatom)
      ya=c(2,iatom)
      za=c(3,iatom)
      xb=c(1,jatom)
      yb=c(2,jatom)
      zb=c(3,jatom)
c
      do 100 katom=1,nat
         xc=c(1,katom)
         yc=c(2,katom)
         zc=c(3,katom)
         cax=xc-xa
         cay=yc-ya
         caz=zc-za
         ca2=cax*cax+cay*cay+caz*caz
         cbx=xc-xb
         cby=yc-yb
         cbz=zc-zb
         cb2=cbx*cbx+cby*cby+cbz*cbz
         ca=sqrt(ca2)
         cb=sqrt(cb2)
c
         nt=0
         kk=0
         lmax=0
         skip=.true.
         do 10 ktype=nbtype+1,ntypes
            if(noprim(katom,ktype).ne.0) then
               skip=.false.
               lmax=max(lmax,maxmom(ktype))
               nt=nt+1
               ntpse(nt)=noprim(katom,ktype)
               do 5 jj=1,noprim(katom,ktype)
                  kk=kk+1
                  nlp(kk)=int(cont(ptcont(katom,ktype)
     $                             +noprim(katom,ktype)+jj-1)+half)
                  clp(kk)=cont(ptcont(katom,ktype)+jj-1)
                  zlp(kk)=ex(ptprim(katom,ktype)+jj-1)
    5          continue
            endif
   10    continue
         if(skip) goto 100
c
         call pseud1 (ntpse,nlp,zlp,clp,
     1                sp,ex,maxi2,maxp2,nx,ny,nz)
         call pseud2 (lmax,ntpse(2),nlp(ntpse(1)+1),zlp(ntpse(1)+1),
     1                clp(ntpse(1)+1),sp,ex,qq,maxi2,maxp2,nx,ny,nz)
c
  100 continue
      call smul(sp,sp,fpi,maxi2*maxp2)
      return
      end
