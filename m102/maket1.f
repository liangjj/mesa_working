*deck @(#)maket1.f	5.1  11/6/94
      subroutine maket1(t,npf,nbf,nprim,ncont,coef,minmom,maxmom,
     $                  nocart,cstart,pstart)
c***begin prologue     maket1.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe,paul (lanl)
c***source             @(#)maket1.f	5.1   11/6/94
c***purpose            generates the primitive to contraction 
c                      transformation matrix.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       maket1.f
      implicit none
c     --- input variables -----
      integer npf,nbf,nprim,ncont,minmom,maxmom
      integer cstart,pstart
c     --- input arrays (unmodified) ---
      integer nocart(0:*)
      real*8 coef(nprim,ncont,minmom:maxmom)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 t(npf,nbf)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer csave,cont,p,prim,angmom,cart,c,ncart
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
