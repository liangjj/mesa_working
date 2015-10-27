*deck compar.f
c***begin prologue     compar
c***date written       941125   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           compare, functions
c***author             schneider, barry (nsf)
c***source             m6202
c***purpose            compare an exact and a legendre decomposed function
c***references         none
c
c***routines called
c***end prologue       compar
      subroutine compar (fnse,fse,fnsa,fsa,nr,nth,nph,nang,nonsep)
      implicit integer (a-z)
      real*8 fnse, fse, fnsa, fsa
      dimension fnse(nr,nang), fse(nr,nth,nph), fnsa(nr,nang)
      dimension fsa(nr,nth,nph)
      character*3 itoc
      character*80 title
      logical nonsep
      common /io/ inp, iout
      if (nonsep) then
          do 10 rpt=1,nr
             do 20 angpt=1,nang
                fnse(rpt,angpt)=abs(fnse(rpt,angpt)-fnsa(rpt,angpt))
   20        continue
   10     continue
          title='point by point absolute error'
          call prntfm(title,fse,nr,nang,nr,nang,iout)                    
      else
          do 30 rpt=1,nr
             do 40 thept=1,nth
                do 50 phipt=1,nph
                   fse(rpt,thept,phipt)=abs(fse(rpt,thept,phipt)
     1                                  - fsa(rpt,thept,phipt))
   50           continue
   40        continue                   
   30     continue
          do 60 phipt=1,nph
             title='point by point absolute error for phi point = '
     1              //itoc(phipt)
             call prntfm(title,fse(1,1,phipt),nr,nth,nr,nth,iout) 
   60     continue               
      endif
      return
      end

