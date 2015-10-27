*deck @(#)frelm.f
c***begin prologue     frelm
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           frelm, link 7000
c***author             schneider, barry (nsf)
c***source             m7000
c***purpose            fill free angular momentum arrays
c***references       
c
c***routines called
***end prologue       frelm
      subroutine frelm(card,nsym,n,l,m,types,ltyp,mtyp,dimsym)
      implicit integer (a-z)
      character *1600 card
      character *8 types
      logical logkey, tstsym
      dimension l(n), m(n), types(dimsym), ltyp(dimsym)
      dimension mtyp(dimsym), nsym(dimsym)
      icnt=0
      do 10 i=1,dimsym
         tstsym=logkey(card,types(i),.false.,' ')
         if (tstsym) then
             icnt=icnt+1
             l(icnt)=ltyp(i)
             m(icnt)=mtyp(i)
             nsym(i)=nsym(i)+1
         endif
   10 continue     
      return
      end























