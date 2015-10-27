*deck vcpl.f
c***begin prologue     vcpl
c***date written       930627   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6236, link 6236
c***author             schneider, b. i.(nsf)
c***source             m6236
c***purpose            read in coupling constants for potential matrix
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       vcpl
      subroutine vcpl(card,cpass,vcij,nc,prnt)
      implicit integer (a-z)
      real*8 vcij, fpkey
      character*3 itoc, chr1, chr2
      character*(*) card, cpass
      character*80 title
      logical prnt, posinp
      dimension vcij(nc*(nc+1)/2)
      common/io/ inp, iout
      if (posinp('$channels',cpass) ) then
          call cardin(card)
      endif
      ij=0
      do 10 i=1,nc
         call pakstr(itoc(i),ii)
         do 20 j=1,i
            ij=ij+1
            call pakstr(itoc(j),jj)
            chr1=itoc(i)
            chr2=itoc(j)
            vcij(ij)=fpkey(card,'v'//chr1(1:ii)//chr2(1:jj),-1.d0,' ')
 20      continue
 10   continue
      if (prnt) then
          title='lower triangle of coupling constants'
          write(iout,1) title
          ij=0
          do 30 i=1,nc
             j=ij+1
             k=ij+i
             write(iout,2) i
             write(iout,3) ( vcij(jj),jj=j,k)  
             ij=k
   30     continue
      endif                    
      return
 1    format (//,a80)
 2    format (/,1x,'row = ',i4)
 3    format( (5x,5e15.8) )      
      end
