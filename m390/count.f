*deck @(#)count.f	1.1  8/1/91
      subroutine count(label,int,lenbuf,values,ptr,lenbin,sort,
     #                  asort,lnsort,n2int,ijklpt,nnqsym,nso,nsym,
     #                  itape,bins,intfil)
c
      implicit integer (a-z)
c
      real*8 int(*),values(lenbin),sort(*)
      integer label(lenbuf),ptr(lenbin),ijklpt(nnqsym),nso(nsym)
      integer asort(lnsort),bins(lenbin)
      character*(*) intfil
c
      common /io/ inp,iout
c
      ioff(i,j)=i*(i-1)/2+j
      ioff2(i,j,k,l)=ioff(ioff(i,j),ioff(k,l))
c
      n=0
      nints=0
      do 800 ism=1,nsym
         do 700 jsm=1,ism
            do 600 ksm=1,ism
               do 500 lsm=1,ksm
c
c                 ----- check for totally symmetric integrals -----
c
                  if (symprd(symprd(ism,jsm),symprd(ksm,lsm)).ne.1)
     #                go to 500
c
c           ----- check for non-canonical symmetry labels -----
c
            if (jsm.gt.ism.or.lsm.gt.ksm.or.(ism.eq.ksm.and.lsm.gt.jsm))
     #                                                              then
               go to 500
            end if
      do 400 ior=1,nso(ism)
         do 300 jor=1,nso(jsm)
            do 200 kor=1,nso(ksm)
               do 100 lor=1,nso(lsm)
c
c           ----- and non-canonical orbital labels -----
c
            if ((ism.eq.jsm.and.jor.gt.ior).or.
     #          (ksm.eq.lsm.and.lor.gt.kor).or.
     #          (ism.eq.ksm.and.kor.gt.ior).or.
     #          (ism.eq.ksm.and.jsm.eq.lsm.and.ior.eq.kor
     #           .and.lor.gt.jor)) then
               go to 100
            end if
c
c           ----- some integrals must be placed two places -----
c
            if (ism.eq.ksm.and.jsm.eq.lsm) then
c
c              ----- ij;kl --> ij;kl & kl;ij
c
               if (ior.eq.kor.and.jor.eq.lor) then
                  n=n+1
               else
                  n=n+2
               end if
            else
c
c              ----- ij;kl --> ij;kl
c
               n=n+1
            end if
   45    continue
  100             continue
  200          continue
  300       continue
  400    continue
  500 continue
  600 continue
  700 continue
  800 continue
c
      write (iout,60) n
   60 format (/,' number of integrals possible: ',i10)
c
c
      return
      end
