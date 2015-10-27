*deck @(#)d2h.f	1.1 9/7/91
c***begin prologue     d2h
c***date written       880814   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m2703, link 2703, symmetry
c***author             schneider, barry (lanl)
c***source             m2703
c***purpose            product symmetries in d2h
c***description        compute d2h symmetry of orbital products
c***                   
c***references         none
c
c***routines called
c***end prologue       d2h
      subroutine d2h(g1,g2,g3)
      implicit integer (a-z)
      character *3 group, g1, g2, g3
      dimension group(8), number(2), table(8,8)
      data group /'a1g','b1g','b2g','b3g','a1u','b1u','b2u','b3u'/
      data table/1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,
     1           4,3,2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,
     2           7,8,5,6,3,4,1,2,8,7,6,5,4,3,2,1/
      do 10 i=1,8
         if (g1.eq.group(i)) number(1)=i
         if (g2.eq.group(i)) number(2)=i
   10 continue
      i1=number(1)
      i2=number(2)
      i3=table(i1,i2)
      g3=group(i3)
      return
      end 
