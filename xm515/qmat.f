*deck @(#)qmat.f	1.1  4/25/95
      subroutine qmat(jmat,kmat,h,nnp,ncoul,nexch,q,pulay,page,
     $                calc,shlmin,shlmax,nshell,nbf,energy,
     $                fcoef,alpha,beta,damp)
c
c***begin prologue     qmat.f
c***date written       850601  
c***revision date      11/6/94      
c   10 june, 1993      rlm at lanl
c      modifying earlier version for kohn-sham
c***keywords           pseudo-canonical, scf, newton-raphson
c***author             saxe, paul(lanl) 
c***source             @(#)qmat.f	1.1   4/25/95
c***purpose            to form 'fock' matrices 
c***description
c     
c   to form the 'q', eqs. 66, 84, or 111 and capital Q, eqs 3,4,24
c   of page and mciver, or the fock matrices of hamilton and pulay.
c***references
c                      m. page and j. w. mciver, jr., jcp 79, 4985 (1983).
c                      t.p. hamilton and p.pulay, jcp 84, 5728(1986).
c
c***routines called
c
c***end prologue       qmat.f
      implicit none
c     --- input variables -----
      character*(*) calc
      logical pulay,page
      integer nnp,nbf,ncoul,nexch,nshell
      real*8 energy
c     --- input arrays (unmodified) ---
      integer shlmin(nshell),shlmax(nshell)
      integer fcoef(nshell),alpha(nshell,nshell),beta(nshell,nshell)
      real*8 jmat(nnp,nshell),kmat(nnp,nshell),h(nnp)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 q(nnp)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      logical called
      integer i,j,k,ishell,jshell,lshell
      integer ii,jj,ij,ia,ja
      integer inp,iout
      real*8 elast,lambda,llast,dedl,denom,numer
      real*8 two,half,zero,one,fact,damp(*)
      real*8 f,fc,fo,fpo,ko,aa,bb,fdiag,fov,fco
c
      parameter (two=2.0d+00, one=1.0d+00, half=0.5d+00, zero=0.0d+00)
c
      common /io/     inp,iout
c
      data called/.false./
      data lambda /1.0d+00/
      save called,lambda
      save elast,dedl,llast
c
c     --- statement functions ---
c
c     --- closed shell fock matrices ---
      f(ij)=h(ij)+two*jmat(ij,1)-kmat(ij,1)
c     --- mike page's open-shell cases ---
      fc(ij)=f(ij)+jmat(ij,2)-half*kmat(ij,2)
      fo(ij)=half*(f(ij)+jmat(ij,2)-kmat(ij,2))
c$$$c     --- pulay's open-shell ---
c$$$      fpo(ij)=f(ij)+half*(two*jmat(ij,2)-kmat(ij,2))
c$$$      ko(ij)=kmat(ij,2)
c
c     -- modified pulay's open-shell to work with alpha and beta k's --
c
      fdiag(ij)=h(ij)+two*jmat(ij,1)+jmat(ij,2)
     $     -.5*(kmat(ij,1)+kmat(ij,2))
      fco(ij)=h(ij)+two*jmat(ij,1)+jmat(ij,2)-kmat(ij,2)
      fov(ij)=h(ij)+two*jmat(ij,1)+jmat(ij,2)-kmat(ij,1)
c
 1000 format (1x,'qmat:negative hessian element:',2i4,f15.9)
c
c
      call rzero(q,nnp)
      if ((calc.eq.'closed').and.page) then
c        --- page's closed-shell scf ---
c
c        --- calculate lambda ---
         if (called) then
c            why the half?
            lambda=half*llast**2*dedl/(elast-energy+llast*dedl)
            if (lambda.lt.zero) lambda=one
         end if
c        calculate d(energy) by d(lambda)
         dedl=zero
         do 20 j=shlmin(2),shlmax(2)
            ja=j*(j-1)/2
            jj=ja+j
            do 10 i=shlmin(1),shlmax(1)
               ij=ja+i
               ii=(i+1)*i/2
               denom=f(jj)-f(ii)
               if (denom.lt.zero) then
                  write (iout,1000) j,i,denom
               end if
               dedl=dedl-4*f(ij)**2/(denom)
   10       continue
   20    continue
         llast=lambda
         elast=energy
c
c        if doing page-mciver, we are in the pseudo-canonical
c        representation. the closed-closed and virtual-virtual
c        blocks have been diagonalized. form the off-diagonal
c        matrix elements and set diagonal elements to orbital indices.
         do 40 i=shlmin(1),shlmax(1)
            ii=(i+1)*i/2
            do 30 j=shlmin(2),shlmax(2)
               ij=j*(j-1)/2+i
               jj=(j+1)*j/2
               numer=f(jj)-f(ii)
               if (numer.lt.zero) then
                  write (iout,1000) j,i,numer
               end if
               q(ij)=f(ij)*(lambda*numer/float(j-i))
   30       continue
   40    continue
         do 60 k=1,nshell
            do 50 i=shlmin(k),shlmax(k)
               ii=(i+1)*i/2
               q(ii)=float(i)
   50       continue
   60    continue
c
c
      else if ((calc.eq.'closed').and.pulay) then
c        --- pulay/roothan closed shell scf ---
c
c           fock matrix is what we want, so make it 
         do 100 i=1,nnp
            q(i)=f(i)
 100     continue
c
c
      else if ((calc.eq.'open').and.page) then
c        --- pages's high-spin open-shell scf ---
c
         if (called) then
            lambda=half*llast**2*dedl/(elast-energy+llast*dedl)
            if (lambda.lt.zero) lambda=one
         end if
         llast=lambda
         elast=energy
c
c        --- calculate small-q matrix and derivative of energy wrt lambda ---
         dedl=zero
c        --- closed-open block ---
         do 160 j=shlmin(2),shlmax(2)
            jj=(j+1)*j/2
            ja=j*(j-1)/2
            do 150 i=shlmin(1),shlmax(1)
               ii=(i+1)*i/2
               ij=ja+i
               denom=fc(jj)-fo(jj)-fc(ii)+fo(ii)
               if (denom.lt.zero) then
                  write (iout,1000) j,i,denom
               end if
               q(ij)=lambda*float(j-i)*(fc(ij)-fo(ij))/(denom)
               dedl=dedl+(fc(ij)-fo(ij))**2/(denom)
  150       continue
  160    continue
c        --- closed-virtual block ---
         do 190 j=shlmin(3),shlmax(3)
            jj=(j+1)*j/2
            ja=j*(j-1)/2
            do 170 i=shlmin(1),shlmax(1)
               ii=(i+1)*i/2
               ij=ja+i
               denom=fc(jj)-fc(ii)
               if (denom.lt.zero) then
                  write (iout,1000) j,i,denom
               end if
               q(ij)=lambda*float(j-i)*fc(ij)/(denom)
               dedl=dedl+fc(ij)**2/(denom)
  170        continue
c           --- open-virtual block ---
            do 180 i=shlmin(2),shlmax(2)
               ii=(i+1)*i/2
               ij=ja+i
               denom=fo(jj)-fo(ii)
               if (denom.lt.zero) then
                  write (iout,1000) j,i,denom
               end if
               q(ij)=lambda*float(j-i)*fo(ij)/(denom)
               dedl=dedl+fo(ij)**2/(denom)
  180        continue
  190     continue
         dedl=-4*dedl
      else if ((calc.eq.'open').and.pulay) then
c        --- pulay's high-spin open-shell ---
c
c        --- diagonal blocks ---
         do 250 ishell=1,nshell
            do 240 i=shlmin(ishell),shlmax(ishell)
               ia=i*(i-1)/2
               do 230 j=shlmin(ishell),i
                  ij=ia+j
                  q(ij)=fdiag(ij)
 230           continue
 240        continue
 250     continue
c
c        --- closed-open block ---
         do 270 i=shlmin(2),shlmax(2)
            ia=i*(i-1)/2
            do 260 j=shlmin(1),shlmax(1)
               ij=ia+j
               q(ij)=fco(ij)
 260        continue
 270     continue
c
c        --- closed-virtual block ---
         do 290 i=shlmin(3),shlmax(3)
            ia=i*(i-1)/2
            do 280 j=shlmin(1),shlmax(1)
               ij=ia+j
               q(ij)=fdiag(ij)
 280        continue
 290     continue
c
c        --- open-virtual block ---
         do 310 i=shlmin(3),shlmax(3)
            ia=i*(i-1)/2
            do 300 j=shlmin(2),shlmax(2)
               ij=ia+j
               q(ij)=fov(ij)
 300        continue
 310     continue
c
c
      else if (calc.eq.'general') then
c
c        --- general case using f's, a's and b's ---
c
c        --- virtual blocks ---
         do 430 i=shlmin(nshell),shlmax(nshell)
            ia=i*(i-1)/2
            do 420 j=shlmin(nshell),i
               ij=ia+j
               q(ij)=half*h(ij)
               do 410 lshell=1,nshell-1
                  aa=fcoef(lshell)
                  bb=-fcoef(lshell)*half
                  q(ij)=q(ij)+
     $                  aa*jmat(ij,lshell)+
     $                  bb*kmat(ij,lshell)
 410           continue
 420        continue
 430     continue
c
c        ---  closed-open blocks ---
         do 480 ishell=1,nshell-1
            fact=damp(ishell)*fcoef(1)/fcoef(ishell)
            do 470 i=shlmin(ishell),shlmax(ishell)
               ia=i*(i-1)/2
               do 460 j=shlmin(ishell),i
                  ij=ia+j
                  q(ij)=fcoef(ishell)*h(ij)
                  do 450 lshell=1,nshell-1
                     q(ij)=q(ij)+
     $                     alpha(ishell,lshell)*jmat(ij,lshell)+
     $                     beta(ishell,lshell)*kmat(ij,lshell)
 450              continue
                  q(ij)=q(ij)*fact
 460           continue
 470        continue
 480     continue
c
c        --- finishing touches ---
         do 540 ishell=2,nshell
            do 530 jshell=1,ishell-1
               do 520 i=shlmin(ishell),shlmax(ishell)
                  ia=i*(i-1)/2
                  do 510 j=shlmin(jshell),shlmax(jshell)
                     ij=ia+j
                     q(ij)=(fcoef(jshell)-fcoef(ishell))*h(ij)
                     do 500 lshell=1,nshell-1
                        q(ij)=q(ij)+
     $                       (alpha(jshell,lshell)-alpha(ishell,lshell))
     $                       *jmat(ij,lshell)+
     $                       (beta(jshell,lshell)-beta(ishell,lshell))*
     $                       kmat(ij,lshell)
 500                 continue
 510              continue
 520           continue
 530        continue
 540     continue
      end if
c
c
      called=.true.
c
c
      return
      end
