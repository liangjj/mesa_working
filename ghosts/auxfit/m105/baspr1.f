*deck @(#)baspr1.f	1.1  11/20/92

      subroutine baspr1(ex,nprim,cf,ncont,nctype,minmom,maxmom,typnam,
     #                  nbtype,iout,pi32)
c
      implicit integer (a-z)
c
      character*(*) typnam(nbtype)
      real*8 ex(nprim),cf(nprim,ncont,minmom:maxmom),pi32
c
 1000 format(' subshell:',a2)
 1010 format(13x,f12.6,5x,(10f10.6))
c
c     ----- start timing -----
c
c
      do 10 angmom=minmom,maxmom
         if(minmom.lt.maxmom) write(iout,1000) typnam(angmom+1)
         do 5 pr=1,nprim
            write(iout,1010) ex(pr),(cf(pr,c,angmom)*sqrt(pi32/
     $               ((2**angmom)*(2*ex(pr))**(angmom+1.5))),c=1,ncont)
    5    continue
   10 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
