*deck @(#)baspr1.f	5.1  11/6/94
      subroutine baspr1(ex,nprim,cf,ncont,nctype,minmom,maxmom,typnam,
     $                  nbtype,iout,pi32)
c***begin prologue     baspr1.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe,paul (lanl)
c***source             @(#)baspr1.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       baspr1.f
      implicit none
c     --- input variables -----
      integer nprim,ncont,minmom,maxmom
      integer nbtype,nctype
      integer iout
      real*8 pi32
c     --- input arrays (unmodified) ---
      real*8 ex(nprim),cf(nprim,ncont,minmom:maxmom)
      character*(*) typnam(nbtype)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer angmom,pr,c
c
 1000 format(' subshell:',a2)
 1010 format(13x,f12.6,5x,(10f10.6))
c
c     --- print the basis information
      do 10 angmom=minmom,maxmom
         if(minmom.lt.maxmom) write(iout,1000) typnam(angmom+1)
         do 5 pr=1,nprim
            write(iout,1010) ex(pr),(cf(pr,c,angmom)*sqrt(pi32/
     $               ((2**angmom)*(2*ex(pr))**(angmom+1.5))),c=1,ncont)
    5    continue
   10 continue
c
c
      return
      end
