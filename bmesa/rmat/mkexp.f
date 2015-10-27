*deck mkexp
c***begin prologue     mkexp
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           bessel, basis
c***author             schneider, barry (nsf)
c***source             
c***purpose            generate exponential basis sets
c***description        
c***references       
c
c***routines called
c***end prologue       mkexp
      subroutine mkexp(f,df,ddf,point,l,nbf,npts,chn,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 f, df, ddf,  point, fpkey, alpha
      character*3 itoc
      character*16 cpass
      character*80 title
      character*80 card
      logical prnt
      dimension f(npts,nbf), df(npts,nbf), ddf(npts,nbf), point(npts)
      call posinp('$expbfn-chn-'//itoc(chn),cpass)
      do 10 i=1,nbf
         read(inp,1) card
         power=intkey(card,'power',0,' ')
         alpha=fpkey(card,'exponent',1.d0,' ')
         pval=power+l+1
         if (pval.eq.0) then
             do 20 j=1,npts
                f(j,i)=exp(-alpha*point(j))
                df(j,i)=-alfa*f(j,i)
                ddf(j,i)=alfa*alfa*f(j,i)
 20          continue
         elseif(pval.eq.1) then
             do 30 j=1,npts
                f(j,i)=point(j)*exp(-alpha*point(j))
                df(j,i)=f(j,i)/point(j)-alfa*f(j,i)   
                ddf(j,i)=-2.d0*alfa*f(j,i)/point(j) + alpha*alpha*f(j,i)
 30          continue
         else
             do 40 j=1,npts
                f(j,i)=(point(j)**pval)*exp(-alpha*point(j))
                df(j,i)=pval*f(j,i)/point(j)-alfa*f(j,i)
                ddf(j,i)=pval*(pval-1)*f(j,i)/(point(j)*point(j)) -
     1                   2.d0*alfa*pval*f(j,i)/point(j)-alfa*alfs*f(j,i)
 40          continue
         endif
 10   continue
      if (prnt) then
          title='basis functions'
          call prntrm(title,f,npts,nbf,npts,nbf,iout)
          title='derivative of basis functions'
          call prntrm(title,df,npts,nbf,npts,nbf,iout)
          title='second derivative of basis functions'
          call prntrm(title,ddf,npts,nbf,npts,nbf,iout)
      endif    
      return
 1    format(a80)
      end
