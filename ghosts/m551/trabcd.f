*deck %W%  %G%
      subroutine trabcd(values,valt,c,nbf,nao)
C
C***Begin prologue
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue
C
      implicit real*8 (a-h,o-z)
      dimension values(2),valt(2),c(nbf,2)
cbigmc
c modified for bigmc
cbigmc
c$$$
c$$$      nll=(l-1)*nao+1
c$$$      do 31 m=1,nao
c$$$         do 30 n=1,m
c$$$            t2mn=t2(m,n)
c$$$            jl=nll
c$$$            do 29 j=1,nao
c$$$               valt(jl)=valt(jl)+t2mn*c(k,j)
c$$$               jl=jl+1
c$$$ 29         continue
c$$$            nll=nll+naonbf
c$$$ 30      continue
c$$$ 31   continue
      call lnkerr('trabcd')
      return
      end
