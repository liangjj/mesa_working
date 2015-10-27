*deck ylmint.f
c***begin prologue     ylmint
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           norm of two spherical harmonics
c***author             schneider, barry (nsf)
c***source             
c***references         none
c
c***routines called
c***end prologue       ylmint
      subroutine ylmint (ylm,yln,ovlp,lval,m,n,nm,nn,nang)
      implicit integer (a-z)
      real*8 ylm, yln, ovlp
      character*80 title
      character*3 itoc
      dimension ylm(nang,m:lval,nm), yln(nang,n:lval,nn)
      dimension ovlp(m:lval,n:lval)
      common /io/ inp, iout
      dimm=lval-m+1
      dimn=lval-n+1
      if (m.eq.n) then
          ii=0
          do 10 i=1,nm
             do 20 j=1,i
                ii=ii+1
                call ebtc(ovlp(m,m),ylm(1,m,i),ylm(1,m,j),dimm,
     1                    nang,dimm)
                title='normalization of spherical harmonics for '//
     1                'mi = mj = '//itoc(m)//'  matrix = '//itoc(ii)
                call prntfm(title,ovlp(m,m),dimm,dimm,dimm,dimm,iout)
   20        continue                
   10     continue
      else
          ii=0
          do 30 i=1,nm
             do 40 j=1,nn
                ii=ii+1
                call ebtc(ovlp(m,n),ylm(1,m,i),yln(1,n,j),dimm,
     1                    nang,dimn)
                title='normalization of spherical harmonics for '//
     1                'mi = '//itoc(m)//' mj = '//itoc(n)//
     2                '  matrix = '//itoc(ii)
                call prntfm(title,ovlp(m,n),dimm,dimn,dimm,dimn,iout)
   40        continue
   30     continue
      endif                                
      return
c
      end


