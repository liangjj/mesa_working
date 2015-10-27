*deck @(#)norm1.f	5.1  11/6/94
      subroutine norm1(ex,cf,nprim,ncont,nmom,minmom,maxmom)
c***begin prologue     norm1
c***date written       850601  (yymmdd)
c***keywords           one-electron, integrals, normalize
c***author             saxe, paul (lanl).
c***source
c***purpose            normalizes the contracted basis functions.
c***description
c                      call norm1(ex,cf,nprim,ncont,nmom,minmom,maxmom)
c
c***references
c***routines called    (none)
c***end prologue       norm1
      implicit integer (a-z)
c
      real*8 ex(nprim),cf(nprim,ncont,minmom:maxmom)
      real*8 fac,pi32
c
c     ----- start timing -----
c
c
c     ----- get pi**1.5 -----
c
      call iosys('read real pi from rwf',1,pi32,0,' ')
      pi32=pi32**1.5d+00
c
      do 50 angmom=minmom,maxmom
         do 10 pr=1,nprim
            do 1 c=1,ncont
               cf(pr,c,angmom)=cf(pr,c,angmom)*sqrt((2**angmom)*
     #                    ((2*ex(pr))**(angmom+1.5))/pi32)
    1       continue
   10    continue
c
         do 20 c=1,ncont
             fac=0.0d+00
            do 12 ipr=1,nprim
               do 11 jpr=1,nprim
                  fac=fac+(cf(ipr,c,angmom)*cf(jpr,c,angmom))/
     #             ((ex(ipr)+ex(jpr))**(angmom+1.5)*2**angmom)
   11          continue
   12       continue
c
            fac=1.0d+00/sqrt(fac*pi32)
             do 13 pr=1,nprim
               cf(pr,c,angmom)=cf(pr,c,angmom)*fac
   13        continue
   20    continue
   50 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
