*deck mkcfun
c***begin prologue     mkcfun
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           bessel, basis
c***author             schneider, barry (nsf)
c***source             
c***purpose            generate basis sets from orthogonal functions
c***description        
c***references       
c
c***routines called
c***end prologue       mkcfun
      subroutine mkcfun(f,df,ddf,x,a,b,left,right,l,nbf,npts,npoly,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 f, df, ddf, x, a, b, left, right 
      character*80 title
      character*2 itoc
      logical prnt
      dimension f(npts,nbf), df(npts,nbf), ddf(npts,nbf)
      dimension a(0:npoly+1), b(0:npoly+1), x(npts)
c     read in the recursion coefficients
      call iosys('read real "polynomial a coefficients for l='//
     1           itoc(l)//'" from lamdat',npoly+1,a,0,' ')
      call iosys('read real "polynomial b coefficients for l='//
     1           itoc(l)//'" from lamdat',npoly+1,b,0,' ')
      call iosys('read integer "order of leading left polynomial '//
     1           'for l='//itoc(l)//'" from lamdat',1,
     2            pleft,0,' ')
      call iosys('read integer "order of leading right polynomial '
     1           //'for l='//itoc(l)//'" from lamdat',1,
     2           pright,0,' ')
c     calculate the polynomials needed
      call gpoly(f,df,ddf,x,a,b,left,right,pleft,pright,nbf,npts,
     1           .false.)     
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
