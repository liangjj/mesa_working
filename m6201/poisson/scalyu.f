*deck scalyu.f
c***begin prologue     scalyu
c***date written       921221   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           scalyu, link m6201, legendre integrals
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            scale an array 
c***references         none
c
c***routines called
c***end prologue       scalyu
      subroutine scalyu (fin,fout,fsc,nrtot,nthet,nphi,nang,
     1                   nonsep,noscal,prnt)
      implicit integer (a-z)
      character*80 title
      real*8 fin, fout, fsc
      logical nonsep, noscal, prnt
      dimension fin(nrtot,*), fout(nrtot,*), fsc(nrtot)
      common/io/inp,iout
      nprd=nthet*nphi
      if (nonsep) then
          nprd=nang
      endif
      if (noscal) then
          call copy(fin,fout,nrtot*nprd)
      else                    
          do 10 pt=1, nrtot
             do 20 ang=1,nprd
                fout(pt,ang) = fin(pt,ang)*fsc(pt)
   20        continue
   10     continue
      endif
      if (prnt) then
          title='scaled yukawa potential'
          call prntfm(title,fout,nrtot,nprd,nrtot,nprd,iout)
      endif           
      return
      end

