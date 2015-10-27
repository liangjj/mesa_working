c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Fillfun}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle 

*deck filfun.f
c***begin prologue     filfun
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            function and derivatives on grid
c***                   
c***references         
c
c***routines called    
c***end prologue       filfun
      subroutine filfun(q,wt,p,dp,ddp,norm,qi,qj,wti,wtj,pi,pj,dpi,dpj,
     1                  ddpi,ddpj,ni,nj,rl,rc,rr,bridge,nfun,
     2                  n,start)
      implicit integer (a-z)
      real*8 q, wt, p, dp, ddp, norm, qi, qj, wti, wtj, pi, pj, dpi, dpj
      real*8 ddpi, ddpj, rl, rc, rr, fac1, fac2, fac3
      logical bridge
      character*80 title
      dimension q(n), wt(n), p(n,n), dp(n,n), ddp(n,n)
      dimension norm(n), qi(ni), qj(nj), wti(ni), wtj(nj)
      dimension pi(ni,ni), pj(nj,nj), dpi(ni,ni), dpj(nj,nj)
      dimension ddpi(ni,ni), ddpj(nj,nj) 
      common/io/inp, iout

c
c        get the normalizations and functions in
c        region i, including the bridge function if present.
c        if there is a bridge function, it contains pieces
c        from region i and region i + 1.
c

      cnti=start
      do 10 i=1,nfun-1
         q(i)=qi(cnti)
         wt(i)=wti(cnti)
         norm(i) = sqrt ( 1.d0/ wt(i) )
         cntj=start
         do 20 j=1,nfun
            p(j,i)=pi(cntj,cnti)*norm(i)
            dp(j,i)=dpi(cntj,cnti)*norm(i)
            ddp(j,i)=ddpi(cntj,cnti)*norm(i)
            cntj=cntj+1
 20      continue   
         cnti=cnti+1
 10   continue
      q(nfun) = qi(cnti)   

c
c     check if the last function in region i is a bridge function
c
      
      if(bridge) then      

c
c        yes it is. 
c

         wt(nfun) = wti(cnti) + wtj(1)
         norm(nfun) = sqrt ( 1.d0 / wt(nfun) )
         cntj=start
         do 30 j=1,nfun
            p(j,nfun)=pi(cntj,cnti)*norm(nfun)
            cntj=cntj+1
 30      continue   

c

      else

c
c        no it is not.
c

         wt(nfun)=wti(cnti) 
         norm(nfun)= sqrt ( 1.d0/ wt(nfun) )
         cntj=start
         do 40 j=1,nfun
            p(j,nfun)=pi(cntj,cnti)*norm(nfun)
            dp(j,nfun)=dpi(cntj,cnti)*norm(nfun)
            ddp(j,nfun)=ddpi(cntj,cnti)*norm(nfun)
            cntj=cntj+1
 40      continue   
      endif
      return
      end       



