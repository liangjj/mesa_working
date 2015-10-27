*deck @(#)maket1.f	1.1  11/20/92
      subroutine maket1(t,npf,nbf,nprim,ncont,coef,minmom,maxmom,
     #                  nocart,cstart,pstart)
c
      implicit integer (a-z)
c
      real*8 t(npf,nbf),coef(nprim,ncont,minmom:maxmom)
      integer nocart(0:*)
c
      csave=cstart
      do 90 cont=1,ncont
         p=pstart
         do 80 prim=1,nprim
            c=csave
            do 75 angmom=minmom,maxmom
               ncart=nocart(angmom)
               do 70 cart=1,ncart
                  t(p+cart,c+cart)=coef(prim,cont,angmom)
   70          continue
               p=p+ncart
               c=c+ncart
               if (p.gt.npf) call lnkerr('problems with p')
               if (c.gt.nbf) call lnkerr('problems with c')
   75       continue
   80    continue
         csave=c
   90 continue
c
c
      return
      end
