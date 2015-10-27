*deck @(#)rdmom.f	1.1 9/7/91
c***begin prologue     rdmom
c***date written       901008   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rdmom, link 6002
c***author             schneider, barry (lanl)
c***source             m6002
c***purpose            input long range moments
c***references         none
c
c***routines called
c***end prologue       rdmom
      subroutine rdmom(mom,nomom,index,nsts,ntri,mxmom)
      implicit integer (a-z)
      real *8 mom
      character *3 itoc
      character *2 ic, jc
      character *5 code
      character*15 cpass
      character *1600 card
      logical test, logkey
      dimension mom(mxmom,ntri), index(3,mxmom,ntri), nomom(ntri)
      common /io/ inp, iout          
      call rzero(mom,mxmom*ntri)
      ij=0
      do 10 is=1,nsts
         ic=itoc(i)
         do 20 js=1,is
            jc=itoc(j)
            code=ic//','//jc
            ij=ij+1
            call posinp('$vlong('//code//')',cpass)
            call cardin(card)
            test=logkey(card,'present',.false.,' ')
            if (test) then
                write (iout,1) is, js
                nomom(ij)=intkey(card,'no-moments',0,' ')
                if (nomom(ij).ne.0) then
                    call fparr(card,'vlr',mom(1,ij),nomom(ij),' ')
                    do 30 imom=1,nomom(ij)
                       call intarr(card,'l-m-n',index(1,imom,ij),3,' ')
   30               continue
                    write (iout,2) (index(1,imom,ij), index(2,imom,ij),
     1                              index(3,imom,ij), mom(imom,ij),
     2                              imom=1,nomom(ij))
                endif
            endif
   20    continue
   10 continue
    1 format (/,1x,'long range moments state',1x,i2,1x,'state',1x,i2)
    2 format (/, 3(10x,'l',1x,i3,1x,'m',1x,i3,1x,'power',1x,i2,1x,
     1                                           'moment',1x,e15.8))
      return
      end
