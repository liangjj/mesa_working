*deck @(#)dmatin.f	1.1 9/7/91
c***begin prologue     dmatin
c***date written       880814   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           dmatin, link 2703
c***author             schneider, barry (lanl)
c***source             m2703
c***purpose            input density matrices and transform
c***                   to ao basis
c***references         none
c
c***routines called
c***end prologue       dmatin
      subroutine dmatin (rhocn,rhomo,trans,temp,temp1,temp2,roexst,
     1                   smllst,nel,ncon,nmo,nsts,ndim,nmotri,dimpr,
     2                   ncntri,file,prnt)
      implicit integer (a-z)
      character *800 card
      character *4 itoc
      character *20 cpass
      character *(*) file, roexst
      logical prnt, test, logkey
      real *8 rhocn, rhomo, trans, temp, temp1, temp2, sumro
      dimension trans (ncon,nmo), rhomo(nmotri), rhocn(ncntri,ndim)
      dimension temp(nmo,nmo), temp1(ncon,nmo), temp2(ncon,ncon)
      dimension roexst(ndim), smllst(dimpr,2)
      common /io/ inp, iout
      call rzero(rhomo,nmotri)
      call rzero(rhocn,ncntri*ndim)
      if (file.eq.'kohndt') then
c----------------------------------------------------------------------c
c              temp and temp2 can be equivalenced in call              c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c           get density matrix from disk in ao representation          c
c----------------------------------------------------------------------c
          ij=0
          do 10 is=1,nsts
             do 20 js=1,is
                ij=ij+1
                call iosys ('read real "ao t1pdm:'//itoc(is)//itoc(js)
     1                      //'" from kohndt',ncon*ncon,temp2,0,' ')
c----------------------------------------------------------------------c
c        in general the density matrix is not symmetric when           c
c        the two states are different. however from the point of       c
c        view of the final integral, where we sum over the orbital     c
c        indices, it is sufficient to take the average of the          c
c        elements and then treat the density matrix as upper           c
c        triangular.                                                   c
c----------------------------------------------------------------------c
                count=0
                do 21 i=1,ncon
                   do 22 j=1,i
                      count=count+1
                      rhocn(count,ij)=temp2(i,j) + temp2(j,i)
   22              continue
                   rhocn(count,ij)=rhocn(count,ij)*.5d+00
   21           continue
                sumro=0.0d+00
                do 45 itest=1,ncntri
                   sumro=sumro+rhocn(itest,ij)
   45           continue
                if (sumro.eq.0.0d+00) then
                    roexst(ij)='no'
                else
                    roexst(ij)='yes'
                endif                    
   20        continue
   10     continue
c----------------------------------------------------------------------c
c          get density matrix from input in mo representation          c
c----------------------------------------------------------------------c
      else
          ij=0
          do 30 is=1,nsts
             do 40 js=1,is
                ij=ij+1
                call posinp('$denmat-'//itoc(is)//itoc(js),cpass)
                call cardin(card)
                test=logkey(card,'not-present',.false.,' ')
                roexst(ij)='yes'
                if (test) then
                    roexst(ij)='no'
                endif
                if (roexst(ij).eq.'yes') then
                    sumro=0.d+00
                    do 100 i=1,nmo
                       call fparr(card,'density-matrix-row-'//itoc(i),
     1                            rhomo,nmo,' ')
                       do 101 j=1,nmo
                          temp(i,j)=rhomo(j)
  101                  continue
                       sumro=sumro+temp(i,i)
  100               continue
c----------------------------------------------------------------------c
c                    normalize density matrix                          c
c----------------------------------------------------------------------c
                    sumro=nel/sumro
                    do 200 i=1,nmo
                       do 210 j=1,nmo
                          temp(i,j)=temp(i,j)*sumro
  210                  continue
  200               continue
c----------------------------------------------------------------------c
c               transform to ao basis                                  c
c----------------------------------------------------------------------c
                    call ebc(temp1,trans,temp,ncon,nmo,nmo)
                    call ebct(temp2,temp1,trans,ncon,nmo,ncon)
                    count=0
                    do 201 i=1,ncon
                       do 202 j=1,i
                          count=count+1
                          rhocn(count,ij)=temp2(i,j) + temp2(j,i)
  202                  continue
                       rhocn(count,ij)=rhocn(count,ij)*.5d+00
  201               continue
                endif
   40        continue  
   30     continue
      endif 
      if (prnt) then
          ii=0 
          do 50 is=1,nsts
             do 60 js=1,is
                ii=ii+1
                write (iout,70) is, js
                call print(rhocn(1,ii),ncntri,ncon,iout)
   60        continue
   50     continue
      endif
      call izero(smllst,2*dimpr)
      do 110 i=1,ncon
         sumro=0.d+00
         do 220 j=1,ncon
            itri=i*(i-1)/2+j
            if ( i.lt.j) then
                 itri=j*(j-1)/2+i
            endif
            ij=0  
            do 310 is=1,nsts
               do 410 js =1,is
                  ij=ij+1
                  sumro=sumro+abs(rhocn(itri,ij))
  410          continue
  310       continue
  220    continue
         if (sumro.gt.1.d-08) then
             smllst(i,1)=i
         endif
  110 continue
      itri=0
      do 120 i=1,ncon
         if (smllst(i,1).eq.0) then
             do 230 j=1,i
                itri=itri+1
                ij=0
                do 330 is=1,nsts
                   do 430 js=1,is
                      ij=ij+1
                      rhocn(itri,ij)=0.d+00
  430              continue
  330           continue    
  230        continue
         else
             itri=itri+i
         endif
  120 continue
      nsmall=0
      do 510 i=1,ncon
         if(smllst(i,1).ne.0) then
            nsmall=nsmall+1
            smllst(nsmall,2)=smllst(i,1)
         endif
  510 continue
   70 format(/,5x,'density matrix i='1x,i3,2x,'j=',1x,i3)     
      return
      end 
