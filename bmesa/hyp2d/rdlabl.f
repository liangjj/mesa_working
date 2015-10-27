*deck rdlabl.f
c***begin prologue     rdlabl
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            read and relabel input so that the order
c***                   is x,y,z,hyperangle,hyperradius,time
c***references         
c
c***routines called    
c***end prologue       rdlabl
      subroutine rdlabl(ops,coord,dim)
      implicit integer (a-z)
      character*(*) ops, coord
      character*128 dummy, chrkey
      character*16 keywrd
      dimension coord(dim), keywrd(10)
      common/io/inp, iout
      data keywrd / 'x','y','z','r1','r2','r3','r','hyperangle',
     1              'hyperradius','time' /
      dummy=chrkey(ops,'coordinates','x,y,z,hyperangle,hyperradius,'//
     1                  'time',' ')
c
c     search for possible keys and re-order
c
      count=0
      do 10 i=1,10
         n=keyloc(dummy,keywrd(i),lockey)
         if(n.ne.0) then
            count=count+1
            coord(count)=keywrd(i)
         endif
 10   continue   
      do 20 i=1,count
         write(iout,1) i, coord(i)
 20   continue    
      return
 1    format(/,1x,'coordinate position = ',i1,1x,
     1            'coordinate type = ',a16)
      end



