*deck @(#)dout32.f	5.1  11/6/94
      subroutine dout32(conint,nconti,ncontj,ncontk,ncontl,lenblk,
     #            symcen,angmom,start,nat,nbtype,nnp,ints,labels,
     #            lenbuf,numbuf,ncint,flbl,numint,plbls,num,
     #            nobf,nactul,cdint,ndcen,npass,last,nder,next)
c
c***vectorised module to attach labels to derivative integrals
c
c  5 july   1991    rlm at lanl
c     32 bit version made from dout64.f
c paul saxe                    7 august 1984                     lanl
c  modified 14 may 1986 by pws at lanl for derivative integrals
c
      implicit integer (a-z)
c
      real*8 conint(numint),ints(lenbuf,nder),cdint(numint,3,ndcen)
      real*8 cutoff
      integer symcen(4),angmom(4),start(nat,nbtype)
      integer labels(2,lenbuf,nder)
      integer flbl(lenblk,4),plbls(2,numint)
      integer nobf(nbtype)
      integer last(nder),numbuf(nder)
c
      common /toler/  cutoff
      common/io/inp,iout
c
c     ----- get the number of functions on each shell -----
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
                  plbls(1,n)=i+flbl(q,1)
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
                  plbls(1,n)=i
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
                     plbls(1,n)=or(shiftl(plbls(1,n),10),j+flbl(q,2))
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
                     plbls(1,n)=or(shiftl(plbls(1,n),10),j)
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
                     plbls(2,n)=k+flbl(q,3)
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
                     plbls(2,n)=k
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
                  plbls(2,n)=or(shiftl(plbls(2,n),10),l+flbl(q,4))
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
                  plbls(2,n)=or(shiftl(plbls(2,n),10),l)
                  n=n+lenblk
   34          continue
   35       continue
   36    continue
      end if
c
c     ----- now for the derivative integrals -----
c
      do 500 dcen=1,ndcen
         if (dcen.eq.1) then
            cen=symcen(1)
         else if (dcen.eq.2.and.(npass.eq.1.or.npass.eq.2)) then
            cen=symcen(2)
         else if (dcen.eq.2.and.npass.eq.3) then
            cen=symcen(3)
         else if (dcen.eq.3) then
            cen=symcen(3)
         else
            call lnkerr('error with centres in dout32')
         end if
c
         do 490 coord=1,3
            der=3*(cen-1)+coord
            do 480 i=1,numint
               if (abs(cdint(i,coord,dcen)).lt.cutoff) go to 480
               numbuf(der)=numbuf(der)+1
               if (numbuf(der).gt.lenbuf) then
                  ints(1,der)=lenbuf
                  labels(1,1,der)=last(der)
                  labels(2,1,der)=0
                 call iosys('write real "unsorted derivative integrals"'
     $                       //' on rdints  without rewinding',
     $                       lenbuf,labels(1,1,der),0,' ')
                 call iosys('write real "unsorted derivative integrals"'
     $                       //' on rdints without rewinding',
     $                       lenbuf,ints(1,der),0,' ')
                  last(der)=next
                  next=next+2*lenbuf
                  nactul=nactul+lenbuf-1
                  numbuf(der)=2
               end if
               ints(numbuf(der),der)=cdint(i,coord,dcen)
               labels(1,numbuf(der),der)=plbls(1,i)
               labels(2,numbuf(der),der)=plbls(2,i)
  480       continue
  490    continue
  500 continue
c
c        ----- and last derivatives by invariance -----
c
         cen=symcen(4)
         do 510 dcen=2,ndcen
            do 509 coord=1,3
               do 508 i=1,numint
                  cdint(i,coord,1)=cdint(i,coord,1)+cdint(i,coord,dcen)
  508          continue
  509       continue
  510    continue
         do 590 coord=1,3
            der=3*(cen-1)+coord
            do 580 i=1,numint
               if (abs(cdint(i,coord,1)).lt.cutoff) go to 580
               numbuf(der)=numbuf(der)+1
               if (numbuf(der).gt.lenbuf) then
                  ints(1,der)=lenbuf
                  labels(1,1,der)=last(der)
                  labels(2,1,der)=0
                 call iosys('write real "unsorted derivative integrals"'
     $                       //' on rdints without rewinding',
     $                       lenbuf,labels(1,1,der),0,' ')
                 call iosys('write real "unsorted derivative integrals"'
     $                      //' on rdints without rewinding',
     $                      lenbuf,ints(1,der),0,' ')
                  last(der)=next
                  next=next+2*lenbuf
                  nactul=nactul+lenbuf-1
                  numbuf(der)=2
               end if
               ints(numbuf(der),der)=-cdint(i,coord,1)
               labels(1,numbuf(der),der)=plbls(1,i)
               labels(2,numbuf(der),der)=plbls(2,i)
  580       continue
  590    continue
c
c
      return
      end
