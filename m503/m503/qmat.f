*deck @(#)qmat.f	5.1  11/6/94
      subroutine qmat(jmat,kmat,h,nnp,ncoul,nexch,q,pulay,calc,shlmin,
     $     shlmax,nshell,nbf,energy,extrap,diagnl,
     $     level,level1,gamma,c,fcoef,alpha,beta,damp)
c
c***begin prologue     qmat
c***date written       850601   (yymmdd)
c***revision date      870513   (yymmdd)
c
c     13 may 1987       pws at lanl
c     adding tcscf capabilities.
c
c     20 april 1987     pws at lanl
c     adding level-shifts in the form of level and level1. the first
c     form sets an increment for the hessian divisor; the second, an
c     absolute minumum for it.
c
c***keywords            pseudocanonical scf, newton-raphson
c***author              saxe, paul (lanl)
c***source              @(#)qmat.f	5.1   11/6/94
c
c***purpose             to form the 'q', eqs. 66, 84, or 111 and
c     capital 'q', eqs 3,4 and 24 of

c     m. page and j. w. mciver, jr., jcp 79, 4985 (1983).
c
c***description
c
c***references          m. page and j. w. mciver, jr., jcp 79, 4985 (1983).
c
c***routines called     (none)
c
c***end prologue        qmat
c
c
      implicit integer (a-z)
c
      real*8 jmat(nnp,ncoul),kmat(nnp,nexch),h(nnp)
      real*8 q(nnp),energy,elast,lambda,llast,dedl
      real*8 denom,small,big,level,level1,gamma,c(*)
      real*8 fcoef(nshell),alpha(nshell,nshell),beta(nshell,nshell)
      real*8 two,half,zero,fact,damp(*)
      real*8 f,fc,fo,fpo,ko,fa,fb,ka,kb,f2a,f2b,f2c,aa,bb
      integer shlmin(nshell),shlmax(nshell)
      character*(*) calc
      logical pulay
c
      parameter (two=2.0d+00, half=0.5d+00, zero=0.0d+00)
c
      common /io/     inp,iout
c
      data call /0/
      data lambda /1.0d+00/
      save call,lambda
c
      save elast,dedl,llast
c
c     ----- statement functions -----
c
c
c     ----- closed shell fock matrices -----
c
      f(ij)=h(ij)+two*jmat(ij,1)-kmat(ij,1)
c
c     ----- mike page's open-shell cases -----
c
      fc(ij)=f(ij)+jmat(ij,2)-half*kmat(ij,2)
      fo(ij)=half*(f(ij)+jmat(ij,2)-kmat(ij,2))
c
c     ----- mike page's tcscf -----
c
      f2a(ij)=c(1)**2*(f(ij)+two*jmat(ij,2)-kmat(ij,2))+
     $     c(1)*c(2)*kmat(ij,3)
      f2b(ij)=c(2)**2*(f(ij)+two*jmat(ij,3)-kmat(ij,3))+
     $     c(1)*c(2)*kmat(ij,2)
      f2c(ij)=f2a(ij)+f2b(ij)-c(1)*c(2)*(kmat(ij,2)+kmat(ij,3))
c
c     ----- pulay's open-shell -----
c
      fpo(ij)=f(ij)+half*(two*jmat(ij,2)-kmat(ij,2))
      ko(ij)=kmat(ij,2)
c
c     ----- pulay's tcscf -----
c
      fa(ij)=c(1)**2*(f(ij)+two*jmat(ij,2)-kmat(ij,2))
      fb(ij)=c(2)**2*(f(ij)+two*jmat(ij,3)-kmat(ij,3))
      ka(ij)=c(1)*c(2)*kmat(ij,2)
      kb(ij)=c(1)*c(2)*kmat(ij,3)
c
      big=0.0d+00
      small=1.0d+30
      call=call+1
      call rzero(q,nnp)
c
      if ((calc.eq.'closed').and.(.not.pulay)) then
c
c        ----- closed-shell scf -----
c
         if (extrap.eq.1) then
c
c        ----- calculate lambda -----
c
            if (call.gt.1) then
               lambda=0.5d+00*llast**2*dedl/(elast-energy+llast*dedl)
               if (lambda.lt.0.0d+00) lambda=1.0d+00
            end if
            dedl=0.0d+00
            do 2 j=shlmin(2),shlmax(2)
               ja=j*(j-1)/2
               jj=ja+j
               do 1 i=shlmin(1),shlmax(1)
                  ij=ja+i
                  ii=(i+1)*i/2
                  denom=f(jj)-f(ii)
                  if (denom.lt.0.0d+00) then
                     write (iout,45) j,i,denom
                  end if
                  if (level1.ne.0.0d+00) denom=max(denom,level1)
c
                  dedl=dedl-4*f(ij)**2/(denom+level)
    1          continue
    2       continue
            llast=lambda
            elast=energy
         end if
c
         do 4 i=shlmin(1),shlmax(1)
            ii=(i+1)*i/2
            do 3 j=shlmin(2),shlmax(2)
               ij=j*(j-1)/2+i
               jj=(j+1)*j/2
               if (diagnl.eq.0) then
                  denom=f(jj)-f(ii)
                  if (denom.lt.0.0d+00) then
                     write (iout,45) j,i,denom
                  end if
                  if (level1.ne.0.0d+00) denom=min(denom,level1)
c
                  q(ij)=f(ij)*(j-i)/(denom+level)
               else
                  q(ij)=f(ij)
               end if
    3       continue
    4    continue
      else if ((calc.eq.'open').and.(.not.pulay)) then
c
c     ----- high-spin open-shell scf -----
c
         if (extrap.eq.1.and.call.gt.1) then
            lambda=0.5d+00*llast**2*dedl/(elast-energy+llast*dedl)
            if (lambda.lt.0.0d+00) lambda=1.0d+00
         end if
         llast=lambda
         elast=energy
c
c     ----- calculate samll-q matrix and derivative of energy wrt lambda -----
c
         dedl=0.0d+00
         do 6 j=shlmin(2),shlmax(2)
            jj=(j+1)*j/2
            ja=j*(j-1)/2
            do 5 i=shlmin(1),shlmax(1)
               ii=(i+1)*i/2
               ij=ja+i
               denom=fc(jj)-fo(jj)-fc(ii)+fo(ii)
               big=max(big,abs(denom))
               small=min(small,abs(denom))
               if (denom.lt.0.0d+00) then
                  write (iout,45) j,i,denom
 45               format ('     !!!! negative hessian element:',2i4,
     $                 f15.9)
               end if
cps   if (level1.ne.0.0d+00) denom=max(denom,level1)
c
               q(ij)=(j-i)*(fc(ij)-fo(ij))/(denom)
               dedl=dedl+(fc(ij)-fo(ij))**2/(denom)
    5       continue
    6    continue
         do 9 j=shlmin(3),shlmax(3)
            jj=(j+1)*j/2
            ja=j*(j-1)/2
            do 7 i=shlmin(1),shlmax(1)
               ii=(i+1)*i/2
               ij=ja+i
               denom=fc(jj)-fc(ii)
               big=max(big,abs(denom))
               small=min(small,abs(denom))
               if (denom.lt.0.0d+00) then
                  write (iout,45) j,i,denom
               end if
               if (level1.ne.0.0d+00) denom=max(denom,level1)
               q(ij)=(j-i)*fc(ij)/(denom+level)
               dedl=dedl+fc(ij)**2/(denom+level)
    7       continue
            do 8 i=shlmin(2),shlmax(2)
               ii=(i+1)*i/2
               ij=ja+i
               denom=fo(jj)-fo(ii)
               big=max(big,abs(denom))
               small=min(small,abs(denom))
               if (denom.lt.0.0d+00) then
                  write (iout,45) j,i,denom
               end if
               if (level1.ne.0.0d+00) denom=max(denom,level1)
               q(ij)=(j-i)*fo(ij)/(denom+level)
               dedl=dedl+fo(ij)**2/(denom+level)
    8       continue
    9    continue
c
         dedl=-4*dedl
c
      else if ((calc.eq.'gvb').and.(.not.pulay)) then
c
c     ----- tcscf -----
c
         call rzero(q,nnp)
c
c     ----- doubly occupied-first of pair ----
c
         i=shlmin(2)
         ia=i*(i-1)/2
         ii=ia+i
         do 610 j=shlmin(1),shlmax(1)
            ij=ia+j
            jj=j*(j+1)/2
            denom=f2c(ii)-f2a(ii)-f2c(jj)+f2a(jj)+level
            if (denom.lt.zero) write (iout,45) i,j,denom
            if (level1.ne.zero) denom=max(denom,level1)
            q(ij)=(i-j)*(f2c(ij)-f2a(ij))/denom
 610     continue
c
c     ----- doubly occupied-second of pair -----
c
         i=shlmin(3)
         ia=i*(i-1)/2
         ii=ia+i
         do 620 j=shlmin(1),shlmax(1)
            ij=ia+j
            jj=j*(j+1)/2
            denom=f2c(ii)-f2b(ii)-f2c(jj)+f2b(jj)+level
            if (denom.lt.zero) write (iout,45) i,j,denom
            if (level1.ne.zero) denom=max(denom,level1)
            q(ij)=(i-j)*(f2c(ij)-f2b(ij))/denom
 620     continue
c
c     ----- first of pair -- second of pair -----
c
         i=shlmin(3)
         ia=i*(i-1)/2
         ii=ia+i
         j=shlmin(2)
         ij=ia+j
         jj=j*(j+1)/2
         denom=f2a(ii)-f2b(ii)-f2a(jj)+f2b(jj)+gamma+level
         if (denom.lt.zero) write (iout,45) i,j,denom
         if (level1.ne.zero) denom=max(denom,level1)
         q(ij)=(i-j)*(f2a(ij)-f2b(ij))/denom
c
c     ----- doubly occupied -- virtuals -----
c
         do 640 i=shlmin(4),shlmax(4)
            ia=i*(i-1)/2
            ii=ia+i
            do 630 j=shlmin(1),shlmax(1)
               ij=ia+j
               jj=j*(j+1)/2
               denom=f2c(ii)-f2c(jj)+level
               if (denom.lt.zero) write (iout,45) i,j,denom
               if (level1.ne.zero) denom=max(denom,level1)
               q(ij)=(i-j)*f2c(ij)/denom
 630        continue
 640     continue
c
c     ----- first of pair -- virtuals -----
c
         j=shlmin(2)
         jj=j*(j+1)/2
         do 650 i=shlmin(4),shlmax(4)
            ia=i*(i-1)/2
            ii=ia+i
            ij=ia+j
            denom=f2a(ii)-f2a(jj)+level
            if (denom.lt.zero) write (iout,45) i,j,denom
            if (level1.ne.zero) denom=max(denom,level1)
            q(ij)=(i-j)*f2a(ij)/denom
 650     continue
c
c     ----- second of pair -- virtuals -----
c
         j=shlmin(3)
         jj=j*(j+1)/2
         do 660 i=shlmin(4),shlmax(4)
            ia=i*(i-1)/2
            ii=ia+i
            ij=ia+j
            denom=f2b(ii)-f2b(jj)+level
            if (denom.lt.zero) write (iout,45) i,j,denom
            if (level1.ne.zero) denom=max(denom,level1)
            q(ij)=(i-j)*f2b(ij)/denom
 660     continue
c
      else if ((calc.eq.'gvb').and.pulay) then
c
c     ----- pulay's tcscf fock matrix -----
c
         call rzero(q,nnp)
         do 130 ishell=1,nshell
            do 120 i=shlmin(ishell),shlmax(ishell)
               ia=i*(i-1)/2
               do 110 j=shlmin(ishell),i
                  ij=ia+j
                  q(ij)=fa(ij)+fb(ij)
 110           continue
 120        continue
 130     continue
c
c     ----- core-active1 -----
c
         do 150 i=shlmin(2),shlmax(2)
            ia=i*(i-1)/2
            do 140 j=shlmin(1),shlmax(1)
               ij=ia+j
               q(ij)=fb(ij)-kb(ij)
 140        continue
 150     continue
c
c     ----- core-active2 -----
c
         do 170 i=shlmin(3),shlmax(3)
            ia=i*(i-1)/2
            do 160 j=shlmin(1),shlmax(1)
               ij=ia+j
               q(ij)=fa(ij)-ka(ij)
 160        continue
 170     continue
c
c     ----- core-virtual -----
c
         do 190 i=shlmin(4),shlmax(4)
            ia=i*(i-1)/2
            do 180 j=shlmin(1),shlmax(1)
               ij=ia+j
               q(ij)=fa(ij)+fb(ij)
 180        continue
 190     continue
c
c     ----- active1-active2 -----
c
         i=shlmin(3)
         j=shlmin(2)
         ij=i*(i-1)/2+j
         q(ij)=fa(ij)-fb(ij)-ka(ij)+kb(ij)
c
c     ----- active1-virtual -----
c
         j=shlmin(2)
         do 200 i=shlmin(4),shlmax(4)
            ij=i*(i-1)/2+j
            q(ij)=fa(ij)+kb(ij)
 200     continue
c
c     ----- active2-virtual -----
c
         j=shlmin(3)
         do 210 i=shlmin(4),shlmax(4)
            ij=i*(i-1)/2+j
            q(ij)=fb(ij)+ka(ij)
 210     continue
c
      else if ((calc.eq.'open').and.pulay) then
c
c     ----- pulay's high-spin open-shell -----
c
         call rzero(q,nnp)
c
         do 250 ishell=1,nshell
            do 240 i=shlmin(ishell),shlmax(ishell)
               ia=i*(i-1)/2
               do 230 j=shlmin(ishell),i
                  ij=ia+j
                  q(ij)=fpo(ij)
 230           continue
 240        continue
 250     continue
c
c     ----- closed-open block -----
c
         do 270 i=shlmin(2),shlmax(2)
            ia=i*(i-1)/2
            do 260 j=shlmin(1),shlmax(1)
               ij=ia+j
               q(ij)=fpo(ij)+half*ko(ij)
 260        continue
 270     continue
c
c     ----- closed-virtual block -----
c
         do 290 i=shlmin(3),shlmax(3)
            ia=i*(i-1)/2
            do 280 j=shlmin(1),shlmax(1)
               ij=ia+j
               q(ij)=fpo(ij)
 280        continue
 290     continue
c
c     ----- open-virtual block -----
c
         do 310 i=shlmin(3),shlmax(3)
            ia=i*(i-1)/2
            do 300 j=shlmin(2),shlmax(2)
               ij=ia+j
               q(ij)=fpo(ij)-half*ko(ij)
 300        continue
 310     continue
c
      else if ((calc.eq.'closed').and.pulay) then
c
c     ----- fock matrix is what we want, so make it -----
c
         do 333 i=1,nnp
            q(i)=f(i)
 333     continue
c
      else if (calc.eq.'general') then
c
c     ----- general case using f's, a's and b's -----
c
         call rzero(q,nnp)
c
c ---- diagonal block
c
c..bhl
c.. virtual-virtual block
c
               do 421 i=shlmin(nshell),shlmax(nshell)
                  ia=i*(i-1)/2
                  do 411 j=shlmin(nshell),i
                     ij=ia+j
                           q(ij)=0.5d+00*h(ij)
                     do 401 lshell=1,nshell-1
                        aa=fcoef(lshell)
                        bb=-fcoef(lshell)*0.5d+00
                        q(ij)=q(ij)+
     $                        aa*jmat(ij,lshell)+
     $                        bb*kmat(ij,lshell)
 401                 continue
 411              continue
 421           continue
c
         do 931 ishell=1,nshell-1
               fact=damp(ishell)*fcoef(1)/fcoef(ishell)
               do 921 i=shlmin(ishell),shlmax(ishell)
                  ia=i*(i-1)/2
                  do 911 j=shlmin(ishell),i
                     ij=ia+j
                     q(ij)=fcoef(ishell)*h(ij)
                     do 901 lshell=1,nshell-1
                        q(ij)=q(ij)+
     $                   alpha(ishell,lshell)*jmat(ij,lshell)+
     $                    beta(ishell,lshell)*kmat(ij,lshell)
 901                 continue
                     q(ij)=q(ij)*fact
 911              continue
 921           continue
 931        continue
c..bhl
         do 440 ishell=2,nshell
            do 430 jshell=1,ishell-1
               jatyp=jshell*(jshell-1)/2
               do 420 i=shlmin(ishell),shlmax(ishell)
                  ia=i*(i-1)/2
                  do 410 j=shlmin(jshell),shlmax(jshell)
                     ij=ia+j
                     q(ij)=(fcoef(jshell)-fcoef(ishell))*h(ij)
                     do 400 lshell=1,nshell-1
                        q(ij)=q(ij)+
     $                       (alpha(jshell,lshell)-alpha(ishell,lshell))
     $                       *jmat(ij,lshell)+
     $                       (beta(jshell,lshell)-beta(ishell,lshell))*
     $                       kmat(ij,lshell)
 400                 continue
 410              continue
 420           continue
 430        continue
 440     continue
c
      end if
c
c     ----- form 'capital q' from 'small q' -----
c
      if(calc.ne.'general') then
         if(.not.pulay) then
            ij=0
            do 91 i=1,nbf
               do 90 j=1,i-1
                  ij=ij+1
                  q(ij)=q(ij)*lambda
 90            continue
               ij=ij+1
               if ((diagnl.eq.0).or.(calc.eq.'closed')) then
                  q(ij)=i
               end if
 91         continue
         end if
      end if
c
c
      return
      end
