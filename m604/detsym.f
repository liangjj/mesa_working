*deck @(#)detsym.f	5.1  11/6/94
      subroutine detsym(salc,nsym,aosym,mosym,aolst,num)
c***begin prologue     detsym
c***date written       910718   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m604, link 604
c***author             schneider barry (lanl)
c***source             m604
c***purpose            to determine and write out the 
c***                   atomic orbitals for each irrep.
c***
c***routines called
c***end prologue       m604
      implicit integer(a-z)
      real*8 salc
      character *1 itoc
      dimension salc(num,num), mosym(nsym), aosym(nsym), aolst(num)
      common /io/ inp, iout
c
c
      write(iout,3) 
      call iosys ('write integer "kohn symmetries" to rwf',1,
     1             nsym,0,' ')
      call iosys ('write integer "kohn symmetries" to chk',1,
     1             nsym,0,' ')
      call izero(aosym,nsym)
      count=0
      do 10 i=1,nsym
         call izero(aolst,num)
         do 20 j=1,mosym(i)
            count=count+1
            do 30 k=1,num
               if (abs(salc(k,count)).gt.1.d-15) then
                   aolst(k)=k
               endif
   30       continue
   20    continue
         jcntr=0
         do 40 j=1,num
            if (aolst(j).ne.0) then
                aosym(i)=aosym(i)+1
                jcntr=jcntr+1
                aolst(jcntr)=aolst(j)   
            endif
   40    continue
         if (aosym(i).ne.0) then
             call iosys ('write integer "symaos-'//itoc(i)//'" to rwf',
     1                    aosym(i),aolst,0,' ')
             call iosys ('write integer "symaos-'//itoc(i)//'" to chk',
     1                    aosym(i),aolst,0,' ')
             write (iout,*) ' '
             write (iout,*) 'symmetry', i
             write(iout,1) aosym(i)
             write (iout,2) (aolst(j),j=1,aosym(i))
         endif
   10 continue
      call iosys ('write integer "kohn ao numsym" to rwf',nsym,
     1             aosym,0,' ')
      return
    1 format(/,10x,'number of symmetry orbitals in this irrep',1x,i4)
    2 format(/,5x,'list of aos',(/,5x,10(i4,1x)))
    3 format(/,5x,'ao orbital symmetry information')
      end
