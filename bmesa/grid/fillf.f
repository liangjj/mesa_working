c $Header: fillf.f,v 1.1 92/12/31 15:23:15 bis Exp $
*deck @(#)fillf.f	1.1 9/9/91
      subroutine fillf(f,rpt,thpt,phpt,wtr,wtth,wtph,a,yukawa,eta,
     1                 nr,nth,nph,dimr,dimth,dimph,ncent,
     2                 dimcen,type,which)
      implicit real *8 (a-h,o-z)
      integer ctoi, lwr, upr
      integer dimr, dimth, dimph, dimcen
      dimension f(nr,nth,nph), rpt(dimr), thpt(dimth), phpt(dimph)
      dimension wtr(dimr), wtth(dimth), wtph(dimph)
      dimension a(3,dimcen),eta(dimcen)
      character*(*) type, which
      common/io/ inp,iout
      lwr=1
      upr=ncent
      if (which.ne.'all') then
          icen=ctoi(which)
          lwr=icen
          upr=icen
      endif
      call rzero(f,nr*nth*nph)
      if (type.eq.'shielded-coulomb') then
          do 10 k=1,nph
             snphi=sin(phpt(k))
             cnphi=cos(phpt(k))
             do 20 j=1,nth
                snth=sqrt(1.d0-thpt(j)*thpt(j))
                do 30 i=1,nr
                   xq=rpt(i)*snth*cnphi
                   yq=rpt(i)*snth*snphi
                   zq=rpt(i)*thpt(j)
                   do 40 l=lwr,upr
                      dist =sqrt((xq-a(1,l))**2 + (yq-a(2,l))**2
     1                                          + (zq-a(3,l))**2)
                      f(i,j,k) = f(i,j,k) + exp(-eta(l)*dist)/dist
   40              continue 
                   yukawa = yukawa + f(i,j,k)*wtr(i)*wtth(j)*wtph(k)
   30           continue
   20        continue
   10     continue
      elseif (type.eq.'exponential') then
          do 50 k=1,nph
             snphi=sin(phpt(k))
             cnphi=cos(phpt(k))
             do 60 j=1,nth
                snth=sqrt(1.d0-thpt(j)*thpt(j))
                do 70 i=1,nr
                   xq=rpt(i)*snth*cnphi
                   yq=rpt(i)*snth*snphi
                   zq=rpt(i)*thpt(j)
                   do 80 l=lwr,upr
                      dist =sqrt((xq-a(1,l))**2 + (yq-a(2,l))**2
     1                                          + (zq-a(3,l))**2)
                      f(i,j,k) = f(i,j,k) + exp(-eta(l)*dist)
   80              continue 
                   yukawa = yukawa + f(i,j,k)*wtr(i)*wtth(j)*wtph(k)
   70           continue
   60        continue
   50     continue
      endif
      return
      end
