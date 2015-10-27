      subroutine angles(ltop,lpmax,ptcf1,ptcf2,len1,len2,z)
      implicit integer(a-z)
c
c     ----- arguments unchanged -----
      integer ltop,lpmax
c     ----- arguments returned
      integer len1,len2
      integer ptcf1(0:ltop,0:ltop,0:ltop)
      integer ptcf2(0:ltop,0:ltop,0:ltop,0:lpmax)
c           the first region of the scratch array z returns
c           the coefficients acoef1,acoef2 of length 2*len1 and 2*len2
      real*8 z(*)
c
c     allocate core
c     determine length needed for acoef arrays
      len1=0
      len2=0
      do 34 i=0,ltop
         do 34 j=0,ltop
            do 34 k=0,ltop
               if((i+j+k).le.ltop) then
                  lbeg=mod(i+j+k,2)
                  do 31 l=lbeg,ltop,2
                     do 31 m=0,l
                        len1=len1+1
  31              continue
                  do 33 lp=0,lpmax
                     do 33 m=0,lp
                        lbeg=max(lp-i-j-k,0)
                        if(mod(lp+i+j+k,2).eq.0) then
                           if(mod(lbeg,2).ne.0) lbeg=lbeg+1
                        else
                           if(mod(lbeg,2).eq.0) lbeg=lbeg+1
                        endif
                        do 32 lambda=lbeg,lp+i+j+k,2
                           do 32 mu=0,lambda
                              len2=len2+1
   32                   continue
   33             continue 
               endif
   34 continue
      ntheta=max(ltop+1,4)
      nphi=ltop+1
      acoef1=1
c     the factor of two is to account for the complex nature of acoef1.
      wt=acoef1+2*len1
      ut=wt+ntheta
      up=ut+ntheta
      pl=up+nphi
      csn=pl+ntheta*(ltop+1)*(ltop+1)
      sn=csn+nphi*(ltop+1)
      fre=sn+nphi*(ltop+1)
      fim=fre+ntheta*nphi
      top=fim+ntheta*nphi
c
      call coef1(ltop,ntheta,nphi,z(acoef1),len1,ptcf1,z(wt),
     $           z(ut),z(pl),z(csn),z(sn),z(fre),z(fim))
c
      ntheta=max(ltop+lpmax+1,4)
      nphi=ltop+lpmax+1
      acoef2=acoef1+2*len1
      wt=acoef2+2*len2
      ut=wt+ntheta
      up=ut+ntheta
      pl=up+nphi
      csn=pl+ntheta*(lpmax+1)*(lpmax+1)*(ltop+lpmax+1)*(ltop+lpmax+1)
      sn=csn+nphi*(lpmax+ltop+lpmax+1)
      fre=sn+nphi*(lpmax+ltop+lpmax+1)
      fim=fre+ntheta*nphi
      top=fim+ntheta*nphi
      call coef2(ltop,lpmax,ntheta,nphi,z(acoef2),len2,ptcf2,z(wt),
     $           z(ut),z(up),z(pll),z(csn),z(sn),z(fre),z(fim))
c
c
      return
      end
