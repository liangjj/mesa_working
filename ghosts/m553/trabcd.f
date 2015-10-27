      subroutine trabcd(values,valt,c,nbf,nao)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
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
