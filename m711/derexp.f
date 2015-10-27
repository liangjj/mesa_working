*deck @(#)derexp.f	5.1  11/6/94
      subroutine derexp(i2,di,imax,jmax,kmax,lmax,mmax,nv,ndcen,d1exp)
c***begin prologue     derexp.f
c***date written       851113  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe,paul(lanl) 
c***source             @(#)derexp.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       derexp.f
      implicit none
c     --- input variables -----
      integer imax,jmax,kmax,lmax,mmax
      integer nv,ndcen
c     --- input arrays (unmodified) ---
      real*8 i2(nv,3,0:imax,0:jmax,0:mmax,0:lmax)
      real*8 d1exp(nv,3,ndcen)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 di(nv,3,0:imax,0:jmax,0:mmax,0:lmax,ndcen)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,j,k,l,coord,centre,m
c
c     --- add derivatives of exponential terms ---
      do 7 i=0,imax
         do 6 j=0,jmax
            do 5 k=0,kmax
               do 4 l=0,lmax
                  do 3 coord=1,3
                     do 2 centre=1,ndcen
                        do 1 m=1,nv
                           di(m,coord,i,j,k,l,centre)=
     $                              di(m,coord,i,j,k,l,centre)+
     $                       d1exp(m,coord,centre)*i2(m,coord,i,j,k,l)
    1                   continue
    2                continue
    3             continue
    4          continue
    5       continue
    6    continue
    7 continue
c
c
      return
      end
