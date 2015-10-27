*deck @(#)out64.f	5.1  11/6/94
      subroutine out64(conint,nconti,ncontj,ncontk,ncontl,lenblk,
     #                  symcen,angmom,start,nat,nbtype,nnp,ints,labels,
     #                  lenbuf,numbuf,ncint,flbl,numint,
     #                  plbls,num,nobf,nactul,drop,pkindx)
c
c***begin prologue     out64
c***date written       840807   (yymmdd)
c***revision date      910702   (yymmdd)
c
c   2 july   1991      rlm at lanl
c      adding the option to drop functions ala bhl.
c***keywords           
c***author             
c***source             @(#)out64.f	5.1   11/6/94
c***purpose            vectorized module to attach labels to integrals 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       out64
c
      implicit integer (a-z)
c
      real*8 conint(numint),ints(lenbuf)
      real*8 cutoff
      integer symcen(4),angmom(4),start(nat,nbtype),labels(lenbuf)
      integer flbl(lenblk,4),plbls(numint)
      integer nobf(nbtype)
      integer pkindx(num)
      logical drop
c
      parameter (mask10=1023)
c
      common /toler/  cutoff
c
c     ----- and the number of functions on each shell -----
c
      numi=nobf(angmom(1))
      numj=nobf(angmom(2))
      numk=nobf(angmom(3))
      numl=nobf(angmom(4))
c
c     ----- get starting function for each angular-momentum function -----
c
      istart=start(symcen(1),angmom(1))-numi
      jstart=start(symcen(2),angmom(2))-numj
      kstart=start(symcen(3),angmom(3))-numk
      lstart=start(symcen(4),angmom(4))-numl
c
c     ----- form the labels for the functions -----
c
      do 24 if=1,numi
         i=if+istart
         p=numj*numk*numl*(if-1)
         do 21 q=1,numl*numk*numj
            p=p+1
            flbl(p,1)=i
   21    continue
         do 23 jf=1,numj
            j=jf+jstart
            p=numl*numk*((jf-1)+numj*(if-1))
            do 22 q=1,numl*numk
               p=p+1
               flbl(p,2)=j
   22       continue
   23    continue
   24 continue
      do 28 lf=1,numl
         l=lf+lstart
         p=lf
         do 25 q=1,numk*numj*numi
            flbl(p,4)=l
            p=p+numl
   25    continue
         do 27 kf=1,numk
            k=kf+kstart
            p=lf+(kf-1)*numl
            do 26 q=1,numj*numi
               flbl(p,3)=k
               p=p+numl*numk
   26       continue
   27    continue
   28 continue
c
c     ----- form total labels -----
c
      if (lenblk.gt.ncontj*ncontk*ncontl) then
         n=0
         do 55 p=1,ncontj*ncontk*ncontl
            do 54 ic=1,nconti
               i=ic*numi
               do 53 q=1,lenblk
                  n=n+1
                  plbls(n)=i+flbl(q,1)
   53          continue
   54       continue
   55    continue
      else
         n=0
         do 58 q=1,lenblk
            ii=flbl(q,1)
            do 57 ic=1,nconti
               i=ii+ic*numi
               n=q+lenblk*(ic-1)
               do 56 p=1,ncontj*ncontk*ncontl
                  plbls(n)=i
                  n=n+lenblk*nconti
   56          continue
   57       continue
   58    continue
      end if
      if (lenblk.gt.ncontk*ncontl) then
         n=0
         do 48 r=1,ncontk*ncontl
            do 47 jc=1,ncontj
               j=jc*numj
               do 46 ic=1,nconti
                  do 45 q=1,lenblk
                     n=n+1
                     plbls(n)=or(shiftl(plbls(n),10),j+flbl(q,2))
   45             continue
   46          continue
   47       continue
   48    continue
      else
         n=0
         do 52 q=1,lenblk
            jj=flbl(q,2)
            do 51 jc=1,ncontj
               j=jc*numj+jj
               do 50 ic=1,nconti
                  n=q+lenblk*((ic-1)+nconti*(jc-1))
                  do 49 p=1,ncontk*ncontl
                     plbls(n)=or(shiftl(plbls(n),10),j)
                     n=n+lenblk*nconti*ncontj
   49             continue
   50          continue
   51       continue
   52    continue
      end if
      if (lenblk.gt.ncontj*nconti) then
         n=0
         do 40 lc=1,ncontl
            do 39 kc=1,ncontk
               k=kc*numk
               do 38 p=1,nconti*ncontj
                  do 37 q=1,lenblk
                     n=n+1
                     plbls(n)=or(shiftl(plbls(n),10),k+flbl(q,3))
   37             continue
   38          continue
   39       continue
   40    continue
      else
         n=0
         do 44 lc=1,ncontl
            do 43 kc=1,ncontk
               kk=kc*numk
               do 42 q=1,lenblk
                  n=lenblk*nconti*ncontj*((kc-1)+ncontk*(lc-1))+q
                  k=kk+flbl(q,3)
                  do 41 p=1,nconti*ncontj
                     plbls(n)=or(shiftl(plbls(n),10),k)
                     n=n+lenblk
   41             continue
   42          continue
   43       continue
   44    continue
      end if
      if (lenblk.gt.ncontk*ncontj*nconti) then
         n=0
         do 33 lc=1,ncontl
            l=lc*numl
            do 32 p=1,ncontk*ncontj*nconti
               do 31 q=1,lenblk
                  n=n+1
                  plbls(n)=or(shiftl(plbls(n),10),l+flbl(q,4))
   31          continue
   32       continue
   33    continue
      else
         do 36 lc=1,ncontl
            ll=lc*numl
            do 35 q=1,lenblk
               l=ll+flbl(q,4)
               n=(lc-1)*lenblk*nconti*ncontj*ncontk+q
               do 34 p=1,ncontk*ncontj*nconti
                  plbls(n)=or(shiftl(plbls(n),10),l)
                  n=n+lenblk
   34          continue
   35       continue
   36    continue
      end if
c
c     ----- sift out noncanonical integrals and transfer to buffers -----
c
      if (drop) then
c
c     drop integrals from final list according to pkindx.
         do 410 i=1,numint
         if (abs(conint(i)).lt.cutoff) go to 410
            io=shiftr(plbls(i),30)
            jo=and(shiftr(plbls(i),20),mask10)
            ko=and(shiftr(plbls(i),10),mask10)
            lo=and(plbls(i),mask10)
c
            idrop=pkindx(io)*pkindx(jo)*pkindx(ko)*pkindx(lo)
         if(idrop.eq.0) goto 410
            numbuf=numbuf+1
            if (numbuf.gt.lenbuf) then
               ints(1)=lenbuf
               call iosys('write real "unsorted ao integrals" '//
     $           'on rints without rewinding',lenbuf,labels,0,' ')
               call iosys('write real "unsorted ao integrals" '//
     $           'on rints without rewinding',lenbuf,ints,0,' ')
               nactul=nactul+lenbuf-1
               numbuf=2
            end if
            ints(numbuf)=conint(i)
            labels(numbuf)=plbls(i)
  410    continue
      else
c
c        normal processing
         do 420 i=1,numint
         if (abs(conint(i)).lt.cutoff) go to 420
c        if (lbls(i,1).lt.lbls(i,2).or.lbls(i,1).lt.lbls(i,3).or.
c    #       (lbls(i,1).eq.lbls(i,3).and.lbls(i,2).lt.lbls(i,4)).or.
c    #       lbls(i,3).lt.lbls(i,4)) go to 420
            numbuf=numbuf+1
            if (numbuf.gt.lenbuf) then
               ints(1)=lenbuf
               call iosys('write real "unsorted ao integrals" '//
     $           'on rints without rewinding',lenbuf,labels,0,' ')
               call iosys('write real "unsorted ao integrals" '//
     $           'on rints without rewinding',lenbuf,ints,0,' ')
               nactul=nactul+lenbuf-1
               numbuf=2
            end if
            ints(numbuf)=conint(i)
            labels(numbuf)=plbls(i)
  420    continue
      endif
c
c
      return
      end
