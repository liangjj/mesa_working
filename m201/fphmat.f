*deck @(#)fphmat.f	5.1  11/6/94
      subroutine fphmat(nvar,fpcycl,dump,vname,pool0,pool1,d1var,d2var,
     $                  xi,h,yold,d1vold,hdotd1,d1doth,mone,mtwo)
c***begin prologue     fphmat.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)fphmat.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       fphmat.f
      implicit none
c     --- input variables -----
      integer nvar,fpcycl
      logical dump
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
      real*8 pool0(nvar),d2var(nvar)
      real*8 xi(nvar),h(nvar,nvar),yold(nvar),d1vold(nvar)
      real*8 hdotd1(nvar),d1doth(nvar),mone(nvar,nvar)
      real*8 mtwo(nvar,nvar)
c     --- output arrays ---
      real*8 pool1(nvar),d1var(nvar)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j
      real*8 zero,one,sqrtd2,zdotd1,hbar
c
      parameter (zero=0.0d+00,one=1.0d+00)
c
      common/io/inp,iout
c
 1000 format(5x,'at step ',i6,' zdotd1=',e20.10,' and hbar=',e20.10)
 1001 format(5x,'fletcher-powell h-matrix:')
c
c     --- we have completed a complete pass over the variables.
c         now form the fletcher-powell approximation to the inverse
c         hessian.
c
c         if this is the first pass, set the h-matrix to the identity.
c         this ensures we begin by going down the path of steepest
c         descent.
      if(fpcycl.eq.0) then
         call rzero(h,nvar*nvar)
         do 20 i=1,nvar
            sqrtd2=sqrt(d2var(i))
            pool1(i)=pool1(i)*sqrtd2
            d1var(i)=d1var(i)/sqrtd2
            h(i,i)=one
   20    continue
      else
c        --- transform variables to xi-space, and move current points and
c            derivatives into yold and d1vold.
         do 30 i=1,nvar
             sqrtd2=sqrt(d2var(i))
             pool1(i)=pool1(i)*sqrtd2
             yold(i)=yold(i)*sqrtd2
             yold(i)=pool1(i)-  yold(i)
             d1var(i)=d1var(i)/sqrtd2
             d1vold(i)=d1vold(i)/sqrtd2
             d1vold(i)=d1var(i)-d1vold(i)
   30    continue
c        --- evaluate d1doth,hdotd1,zdotd1,and hbar.
         call rzero(hdotd1,nvar)
         call rzero(d1doth,nvar)
         do 40 i=1,nvar
             do 40 j=1,nvar
                    hdotd1(i)=hdotd1(i)+(h(i,j)*d1vold(j))
                    d1doth(i)=d1doth(i)+(h(i,j)*d1vold(j))
   40    continue
         zdotd1=zero
         hbar  =zero
         do 50 i=1,nvar
             zdotd1=zdotd1+(yold(i)*d1vold(i))
             hbar  =hbar  +(d1vold(i)*hdotd1(i))
   50    continue
c        --- now, put all the information together and evaluate the
c            corrective matrices mone and mtwo and combine them with
c            the old h-matrix to get h(fpcycl+1).
         do 60 i=1,nvar
             do 60 j=1,nvar
                 mone(i,j)=(yold(i)*yold(j))/zdotd1
                 mtwo(i,j)=(-hdotd1(i)*d1doth(j))/hbar
   60            h(i,j)=h(i,j)+mone(i,j)+mtwo(i,j)
      endif
c
c     --- print the h-matrix.
      if(dump) then
         write(iout,1000) fpcycl,zdotd1,hbar
         write(iout,1001)
         call matprt(h,nvar,nvar,nvar,nvar,1,1,vname,vname,1,h,.false.)
      endif
c
c     --- now evaluate the new variables and set up for a scaling
c         sequence.
      call rzero(xi,nvar)
      do 70 i=1,nvar
         do 65 j=1,nvar
            xi(i)=xi(i)-(h(i,j)*d1var(j))
   65    continue
   70 continue
c
c     --- set up pool0 in xi-space and transform to x-space.
      call vwxs(pool0,pool1,xi,one,1,nvar)
      do 80 i=1,nvar
         pool0(i)=pool0(i)/sqrt(d2var(i))
         pool1(i)=pool1(i)/sqrt(d2var(i))
         d1var(i)=d1var(i)*sqrt(d2var(i))
         if(fpcycl.gt.0) then
            yold(i)=yold(i)/sqrt(d2var(i))
            d1vold(i)=d1vold(i)*sqrt(d2var(i))
         endif
   80 continue
c
c
      return
      end
