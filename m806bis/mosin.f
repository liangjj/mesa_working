*deck @(#)mosin.f
      subroutine mosin(bfsym,mosym,molst,nsym,nbf)
      implicit integer(a-z)
      character*1 itoc
      dimension bfsym(nbf), mosym(nsym), molst(nbf,nsym)
      common /io/ inp, iout
      call iosys ('read integer "scattering mo numsym" from rwf',
     #                nsym,mosym,0,' ')
      write(iout,1) 
      do 10 i=1,nsym
         if(mosym(i).gt.0) then
            write (iout,2) i, mosym(i)
            call iosys ('read integer "symmos-'//itoc(i)//'" from rwf',
     1                   mosym(i),molst(1,i),0,' ')
            write(iout,3) (molst(j,i), j=1,mosym(i))
            do 20 j=1,mosym(i)
              bfsym(molst(j,i))=i-1
 20        continue   
         endif
 10   continue
 1    format(/,1x,'mo orbital symmetry information read in from rwf')
 2    format(/,1x,'symmetry number             = ',i2,/,1x,
     1            'number of symmetry orbitals = ',i4)
 3    format(/,5x,'list of mos',(/,5x,10(i4,1x)))
      return
      end
