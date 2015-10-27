*deck @(#)c2v.f	1.1 9/7/91
c***begin prologue     c2v      
c***date written       880814   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m2703, link 2703, symmetry
c***author             schneider, barry (lanl)
c***source             m2703
c***purpose            product symmetries in c2v
c***description        compute c2v symmetry of orbital products
c***                   
c***references         none
c
c***routines called
c***end prologue       c2v
      subroutine c2v(g1,g2,g3)
      implicit integer (a-z)
      character *3 group, g1, g2, g3
      dimension group(4), number(2), table(4,4)
      data group /'a1','a2','b1','b2'/ 
      data table/1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1/
      do 10 i=1,4
         if (g1.eq.group(i)) number(1)=i
         if (g2.eq.group(i)) number(2)=i
   10 continue
      i1=number(1)
      i2=number(2)
      i3=table(i1,i2)
      g3=group(i3)
      return
      end
