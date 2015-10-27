*deck @(#)prefac.f	5.1  11/6/94
      subroutine prefac(ar,a,nij,br,b,nkl,expon,index,ijindx,klindx,
     #                  kl,len,lenv,thresh,pi252,t1)
c
      implicit integer (a-z)
c
      real*8 ar(nij),a(nij),br(nkl),b(nkl),expon(len),t1(nij)
      real*8 thresh,pi252
      real*8 brkl,bkl
      integer index(len,6),ijindx(nij,2),klindx(nkl,2)
c
      lenv=0
c
    1 continue
c
         kl=kl+1
         brkl=br(kl)
         bkl=b(kl)
         do 2 ij=1,nij
            t1(ij)=exp(-ar(ij)-brkl)*pi252/(a(ij)*bkl*sqrt(a(ij)+bkl))
    2    continue
c
         lenvsv=lenv
         do 3 ij=1,nij
            if (abs(t1(ij)).gt.thresh) then
               lenv=lenv+1
            expon(lenv)=t1(ij)
               index(lenv,5)=ij
            end if
    3    continue
c
         k=klindx(kl,1)
         l=klindx(kl,2)
         do 4 ij=lenvsv+1,lenv
            index(ij,3)=k
            index(ij,4)=l
            index(ij,6)=kl
    4    continue
c
      if (len-lenv.ge.nij.and.kl.lt.nkl) go to 1
c
c     ----- fill in the i and j indeices in index -----
c
      do 5 ij=1,lenv
         index(ij,1)=ijindx(index(ij,5),1)
         index(ij,2)=ijindx(index(ij,5),2)
    5 continue
c
c
      return
      end
