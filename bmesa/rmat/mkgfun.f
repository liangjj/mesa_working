*deck mkgfun
c***begin prologue     mkgfun
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           gaussian, basis
c***author             schneider, barry (nsf)
c***source             
c***purpose            generate basis sets from gaussians
c***description        
c***references       
c
c***routines called
c***end prologue       mkgfun
      subroutine mkgfun(f,df,ddf,point,alpha,n,nbf,npts,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 f, df, ddf, point, alpha, fac
      character*80 title
      logical prnt
      dimension f(npts,nbf), df(npts,nbf), ddf(npts,nbf)
      dimension point(npts), alpha(nbf), n(nbf)
      do 10 bf=1,nbf
         do 20 pt=1,npts
            fac=exp(-alpha(bf)*point(pt)*point(pt))
            f(pt,bf)=fac
            df(pt,bf)=fac
            ddf(pt,bf)=fac
 20      continue   
         if (n(bf).gt.1) then
             do 30 pt=1,npts
                pp2=point(pt)**(n(bf)+2)
                pp1=point(pt)**(n(bf)+1)
                pp=point(pt)**n(bf)
                ppm1=point(pt)**(n(bf)-1)
                ppm2=point(pt)**(n(bf)-2)
                f(pt,bf)=pp*f(pt,bf)
                df(pt,bf)=df(pt,bf)*( n(bf)*ppm1 - 2.d0*alpha(bf)*pp1 )
                ddf(pt,bf)=ddf(pt,bf) * ( n(bf)*(n(bf)-1)*ppm2 -
     1                                    2.d0*alpha(bf)*
     2                                    (n(bf)+n(bf)+1)*pp +
     3                                    4.d0*alpha(bf)*alpha(bf)*pp2 )
 30          continue   
         endif
         if (n(bf).eq.1) then
             do 40 pt=1,npts
                f(pt,bf)=point(pt)*f(pt,bf)
                df(pt,bf)=df(pt,bf)*( 1.d0 - 
     1                           2.d0*alpha(bf)*point(pt)*point(pt) )
                ddf(pt,bf)=ddf(pt,bf) * ( -6.d0*alpha(bf)*point(pt) +
     1                                     4.d0*alpha(bf)*alpha(bf)*
     2                                     point(pt)*point(pt)*
     3                                              point(pt) )
 40          continue   
         endif
         if (n(bf).eq.0) then
             do 50 pt=1,npts
                df(pt,bf)=df(pt,bf)*( -2.d0*alpha(bf)*point(pt) )
                ddf(pt,bf)=ddf(pt,bf) * ( -2.d0*alpha(bf) +
     1                                     4.d0*alpha(bf)*alpha(bf)*
     2                                     point(pt)*point(pt) )
 50          continue   
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
      end
