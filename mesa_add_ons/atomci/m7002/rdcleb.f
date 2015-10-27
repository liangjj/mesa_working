*deck rdcleb
      subroutine rdcleb (cleb,acleb,icleb)
      implicit integer (a-h,o-z)
      real *8 cleb
      dimension cleb(*), icleb(*), acleb(*)
      data msk1, msk6 /1,63/
      common /io/ inp, iout
      write(iout,100)
c      
      call iosys('read integer "cg buffer size" from atomci',1,
     1            bufsiz,0,' ')
      call iosys('read integer "no. cg writes" from atomci',1,nowrts,
     1            0,' ')
      call iosys('read integer "size of last cg write" from atomci',1,
     1            last,0,' ')      
      do 10 i=1,nowrts
         number=bufsiz
         if (i.eq.nowrts) then
             number=last
         endif
         call iosys('read integer "cg coef" from atomci without '//
     1              'rewinding',number,icleb,0,' ')
         call iosys('read integer "cg coef" from atomci without '//
     1              'rewinding',2*number,acleb,0,' ')
         do 20 j=1,number
c**********************************************************************c
c                   signs of m1 and m2                                 c
            signm2=and(icleb(j),msk1)
            snm2=1
            if (signm2.eq.1) then
                snm2=-1
            endif
            signm1=and(shiftr(icleb(j),1),msk1)
            snm1=1
            if (signm1.eq.1) then
                snm1=-1
            endif
c                   l1,l2,abs(m1),abs(m2) and l
c
            l=and(shiftr(icleb(j),2),msk6)
            m2=snm2*and(shiftr(icleb(j),8),msk6)
            m1=snm1*and(shiftr(icleb(j),14),msk6)
            l2=and(shiftr(icleb(j),20),msk6)
            l1=and(shiftr(icleb(j),26),msk6)
            m=m2-m1
            write(iout,200) l1,l2,l,m1,m2,m,cleb(j)
   20    continue      
   10 continue      
c**********************************************************************c
      return
  100 format(/,7x,' l1 ',1x,' l2 ',1x,' l ',1x,' m1 ',1x,' m2 ',
     1         2x,' m ',15x,' cg coef ')
  200 format(5x,i4,1x,i4,1x,i4,1x,i4,1x,i4,1x,i4,10x,e15.8) 
      end

