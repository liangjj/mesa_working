*deck @(#)trans3.f	1.1 4/28/92
      subroutine trans3(orig,final,ia,ib,ic,fa,fb,fc,a,b,c,
     $                  t1,t2,lenblk,
     $                  imin,imax,jmin,jmax,kmin,kmax,ncart)
c***begin prologue     trans3
c***date written       850601  (yymmdd)
c***revision date      920417  (yymmdd)
c      april 17,1992   rlm at lanl
c         modified to transform 3-center primitive blocks.
c***keywords           one-electron, integrals, trnasformation
c***author             saxe, paul and martin, richard (lanl)
c***source
c***purpose            
c                      apply contraction coefficients to three-center primitive
c                      integrals
c***description
c                      call trans3(orig,final,ia,ib,ic,fa,fb,fc,
c                                  a,b,c,t1,t2,lenblk,
c                                  imin,imax,jmin,jmax,kmin,kmax,ncart)
c
c***references
c***routines called    ebtc(math),ebc(math)
c***end prologue       trans3
      implicit integer (a-z)
c
c     ----- input arrays(unmodified) -----
      real*8 orig(ia,ib,ic,lenblk)
      real*8 a(ia,fa,imin:imax),b(ib,fb,jmin:jmax),c(ic,fc,kmin:kmax)
      integer ncart(0:*)
c
c     ----- output arrays -----
      real*8 final(fa,fb,fc,lenblk)
c
c     ----- scratch arrays -----
      real*8 t1(fa*ib),t2(fa*fb,ic)
c
      common/io/inp,iout
c
c     ----- first transform the ij indices for all auxiliary functions k.
c           then transform all auxiliary functions k.
      ij=0
      do 100 kmom=kmin,kmax
         do 90 k=1,ncart(kmom)
            do 80 imom=imin,imax
               do 70 i=1,ncart(imom)
                  do 60 jmom=jmin,jmax
                     do 50 j=1,ncart(jmom)
                        ij=ij+1
                        do 40 kprim=1,ic
c
c     -----                transform the ij indices for all auxiliary 
c                          primitives k.
                           call ebtc(t1,a(1,1,imom),orig(1,1,kprim,ij),
     $                               fa,ia,ib)
                           call ebc(t2(1,kprim),t1,b(1,1,jmom),
     $                              fa,ib,fb)
   40                   continue
c
c     -----             transform all auxiliary primitives k.
                        call ebc(final(1,1,1,ij),t2(1,1),
     $                           c(1,1,kmom),fa*fb,ic,fc)
   50                continue
   60             continue
   70          continue
   80       continue
   90    continue
  100 continue
c
c
      return
      end
