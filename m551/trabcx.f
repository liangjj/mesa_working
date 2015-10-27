*deck @(#)trabcx.f	1.1  11/30/90
      subroutine trabcx(c,valt,t2,nbf,nco,nao,nocc,k,l)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)trabcx.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
cbigmc
c modified for bigmc
cbigmc
      implicit real*8(a-h,o-z)
      dimension c(nbf,2),valt(2),t2(nocc,nocc)
c
      naonbf=nao*nbf
      if(k.eq.l) go to 28
      nkk=k
      nll=l
      do 27 m=1,nao
         do 26 n=1,m
            t2mn=t2(m,n)
            jk=nkk
            jl=nll
            do 25 j=1,nao
               valt(jk)=valt(jk)+t2mn*c(l,j)
               valt(jl)=valt(jl)+t2mn*c(k,j)
               jl=jl+nbf
               jk=jk+nbf
 25         continue
            nkk=nkk+naonbf
            nll=nll+naonbf
 26      continue
 27   continue
c
      go to 32
c
 28   continue
      nll=l
      do 31 m=1,nao
         do 30 n=1,m
            t2mn=t2(m,n)
            jl=nll
            do 29 j=1,nao
               valt(jl)=valt(jl)+t2mn*c(k,j)
               jl=jl+nbf
 29         continue
            nll=nll+naonbf
 30      continue
 31   continue
c
 32   continue
      return
      end
