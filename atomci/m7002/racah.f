*deck racah
c***begin prologue     m7002
c***date written       920623   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7002,clebsch-gordan coefficients
c***author             schneider, barry (nsf)
c***source             m7002
c***purpose            to calculate clebsch-gordan coefficients.
c***description        non zero cg coefficients calculated and stored
c***                   with indices in packed form.
c
c***references         none
c
c***routines called clebgd
c***end prologue       racah
      subroutine racah (cleb,acleb,icleb,lmax,bufsiz,nowrts,nolst)
      implicit integer (a-z)
      real *8 pre, prefac, rac1, rac2, cleb
      dimension cleb(*), icleb(*), acleb(*)
      common /io/ inp, iout
c**********************************************************************c
c                    acleb and cleb are equivalenced in                c 
c                    the calling routine to facilitate                 c
c                    iosys writes of mixed integer and real            c
c                    variable on a single file.                        c
c**********************************************************************c
      nowrts=0
      count=0
      words=0
      do 10 l1=0,lmax
         l1fac=2*l1+1
         tl1=l1+l1 
         do 20 l2=0,l1
            tl2=l2+l2
            l2fac=2*l2+1
            ltmax=l1+l2
            ltmin=abs(l1-l2)
            do 30 l=ltmin,ltmax
               do 40 m1=-l1,l1 
                  tm1=m1+m1
                  m1fac=1
                  m1abs=abs(m1)
                  mtst=m1abs-2*(m1abs/2)
                  if (mtst.eq.1) then
                      m1fac=-1
                  endif
                  setsm1=0
                  if (m1.lt.0) then
                      setsm1=1
                  endif
                  do 50 m2=-l2,l2
                     m2abs=abs(m2)
                     setsm2=0
                     if (m2.lt.0) then
                         setsm2=1
                     endif
                     tm2=m2+m2
                     mtot=m2-m1
                     tm=mtot+mtot
                     mtabs=abs(mtot)
                     if (count.eq.bufsiz) then
                         nowrts=nowrts+1
                         call iosys('write integer "cg coef" to '//
     1                              'atomci without rewinding',
     2                               bufsiz,icleb,0,' ')                     
                         call iosys('write integer "cg coef" to '//
     1                              'atomci without rewinding',
     2                               2*bufsiz,acleb,0,' ')                     
                         count=0
                     endif
                     if(l.ge.mtabs) then
                        tl=l+l
                        pre=l1fac*l2fac
                        prefac=m1fac*sqrt(pre)/(2*l+1)
                        call clebgd(tl1,tl2,tl,0,0,0,1,rac1)
                        if (rac1.ne.0.d0) then
                            call clebgd(tl1,tl2,tl,-tm1,tm2,tm,1,rac2)
                            if (rac2.ne.0.d0) then
                                words=words+1
                                count=count+1
                                icleb(count)=shiftl(l1,26)
                                icleb(count)=or(shiftl(l2,20),
     1                                          icleb(count))
                                icleb(count)=or(shiftl(m1abs,14),
     1                                          icleb(count))
                                icleb(count)=or(shiftl(m2abs,8),
     1                                          icleb(count))
                                icleb(count)=or(shiftl(l,2),
     1                                          icleb(count))
                                icleb(count)=or(shiftl(setsm1,1),
     1                                          icleb(count))
                                icleb(count)=or(setsm2,icleb(count))
                                cleb(count)=prefac*rac1*rac2
                            endif
                        endif
                     endif
   50             continue
   40          continue
   30       continue      
   20    continue         
   10 continue     
      if (count.gt.0) then
          nowrts=nowrts+1
          nolst=count
          call iosys('write integer "cg coef" to '//
     1               'atomci without rewinding',count,icleb,0,' ')
          call iosys('write integer "cg coef" to '//
     1               'atomci without rewinding',2*count,acleb,0,' ')
      endif
      write(iout,100) nowrts, words
  100 format(//,5x,'number writes and number of words to cg file',1x,i5,
     1     2x,i6)         
      return
      end




