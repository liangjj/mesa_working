*deck mkcfun
c***begin prologue     mkcfun
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           bessel, basis
c***author             schneider, barry (nsf)
c***source             
c***purpose            generate basis sets from ricatti-bessel functions
c***description        
c***references       
c
c***routines called
c***end prologue       mkcfun
      subroutine mkcfun(f,df,ddf,point,root,x,j,jp,y,yp,wron,
     1                 scr,l,nbf,npts,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 f, df, ddf,  point, root, x, j, jp, y, yp, wron, tmp
      character*80 title
      logical prnt
      dimension f(npts,nbf), df(npts,nbf), ddf(npts,nbf)
      dimension root(nbf), point(npts)
      dimension x(npts), j(npts,0:*), jp(npts,0:*), y(npts,0:*)
      dimension yp(npts,0:*), scr(*), wron(*)
      tmp=l*30.d0
      ltop=l+sqrt(tmp)
      ltop=max(ltop,l)
      do 10 iroot=1,nbf
         call vscale(x,point,root(iroot),npts)
         call rcbes(x,j,jp,y,yp,wron,scr,npts,l,ltop,
     1                            'derivatives',.false.)
         call copy(j(1,l),f(1,iroot),npts)
         call copy(jp(1,l),df(1,iroot),npts)
         call secder(j(1,l),ddf(1,iroot),x,l,npts)
         call sscal(npts,root(iroot),df(1,iroot),1)
         call sscal(npts,root(iroot)*root(iroot),ddf(1,iroot),1)           
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
