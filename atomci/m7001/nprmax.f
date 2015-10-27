*deck @(#)nprmax.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001, spline
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            determine size of largest block needed 
c***                   in two electron integral evaluation of f(k)
c***
c***                 
c***                 
c***references       
c
c***routines called    none
c***end prologue       m7001
      subroutine nprmax(nsyma,nsymb,typea,typeb,dimsym,old)
      implicit integer (a-z)
      character *(*) typea, typeb
      character *8 ylm
      dimension nsyma(dimsym), nsymb(dimsym)
      common /io/ inp, iout
      lena=length(typea)      
      call iosys ('read integer "'//typea(1:lena)//' orthogonal '//
     1            'symmetry list" from atomci',dimsym,nsyma,0,' ')
      if (typea.eq.typeb) then
          maxpr=0   
          do 10 i=1,dimsym
             if (nsyma(i).ne.0) then
                 call lmtyp(i,l1,m1,ylm) 
                 do 20 j=1,i
                    if (nsyma(j).ne.0) then
                        call lmtyp(j,l2,m2,ylm)
                        npr=nsyma(i)*nsyma(j)
                        if (i.eq.j) then
                            npr=nsyma(i)*(nsyma(i)+1)/2
                        endif
                        lsum=l1+l2
                        ldiff=abs(l1-l2)
                        deltal=lsum-ldiff+1
                        tpair=npr+deltal*npr
                        maxpr=max(maxpr,tpair)
                    endif
   20            continue
             endif
   10     continue
      else
          lenb=length(typeb)
          call iosys ('read integer "'//typeb(1:lenb)//' orthogonal '//
     1                'symmetry list" from atomci',dimsym,nsymb,0,' ')
          maxpr=0   
          do 30 i=1,dimsym
             if (nsyma(i).ne.0) then
                 call lmtyp(i,l1,m1,ylm) 
                 do 40 j=1,dimsym
                    if (nsymb(j).ne.0) then
                        call lmtyp(j,l2,m2,ylm)
                        npr=nsyma(i)*nsymb(j)
                        lsum=l1+l2
                        ldiff=abs(l1-l2)
                        deltal=lsum-ldiff+1
                        tpair=npr+deltal*npr
                        maxpr=max(maxpr,tpair)
                    endif
   40            continue
             endif
   30     continue
      endif
      maxpr=max(maxpr,old)
      old=maxpr
      call iosys ('write integer "largest '//typea(1:lena)//'-'//
     1            typeb(1:lenb)//' pair block" to atomci',1,maxpr,0,' ')
      return
      end







