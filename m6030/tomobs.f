*deck tomobs
c***begin prologue     tomobs
c***date written       890529   (yymmdd)
c***revision date      920127   (yymmdd)
c***keywords           transformation, kohn
c***author             schneider, barry (lanl)
c***source             m6006
c***purpose            bound-bound ao to mo transformation
c***
c***references         none
c
c***routines called    cebc
c***end prologue       tomobs
 
      subroutine tomobs(aomat,momat,aomatc,momatc,vec,scri,scrj,
     1                  scr,scrc,ncon,nmo,ngi,lgi,nmoi,lmoi,ngj,lgj,
     2                  nmoj,lmoj,type)
      implicit integer (a-z)
      real *8 aomat, momat, vec, scri, scrj, scr
      complex *16 aomatc, momatc, scrc
      character *8 type
      dimension aomat(ngi,ngj), momat(nmoi,nmoj), vec(ncon,nmo)
      dimension scri(ngi,nmoi), scrj(ngj,nmoj), lgi(ngi), lgj(ngj)
      dimension lmoi(nmoi), lmoj(nmoj), scr(ngi,nmoj) 
      dimension scrc(ngi,nmoj), aomatc(ngi,ngj), momatc(nmoi,nmoj)
c----------------------------------------------------------------------c
c      compact the transformation matrices for the i,j blocks          c
c      to use the standard matrix multiply routines                          c
c----------------------------------------------------------------------c
      do 10 i=1,nmoi
	 ii=lmoi(i)
	 do 20 j=1,ngi
	    jj=lgi(j)
	    scri(j,i)=vec(jj,ii)
   20    continue
   10 continue
      do 30 i=1,nmoj
	 ii=lmoj(i)
	 do 40 j=1,ngj
	    jj=lgj(j)
	    scrj(j,i)=vec(jj,ii)
   40    continue
   30 continue
c---------------------------------------------------------------------c
c                           +                                         c
c              now perform u  a  u                                    c
c                                                                     c
c---------------------------------------------------------------------c
      if (type.eq.'real') then
          call ebc(scr,aomat,scrj,ngi,ngj,nmoj)
          call ebtc(momat,scri,scr,nmoi,ngi,nmoj)
      elseif(type.eq.'complex') then
          call ecbc(scrc,aomatc,scrj,ngi,ngj,nmoj)
          call ebtcc(momat,scri,scrc,nmoi,ngi,nmoj)
      else
	  call lnkerr('wrong matrix type in tomobs')
      endif
      return
      end
 
 
 
