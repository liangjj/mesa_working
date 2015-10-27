*deck @(#)extcou.f	1.1 9/8/91
c***begin prologue     m7000
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7000, link 7000, spline
c***author             schneider, b. (nsf)
c***source             m7004
c***purpose            get values of coulomb functions at arbitrary
c***                   points.
c***description        series and spline fits are used to get values
c***                   of regular and irregular coulomb functions
c***                   at arbitrary points. the spline fits are used
c***                   for all values of the regular function and for
c***                   those values of the irregular function where
c***                   the independent variable is above rswtch. for
c***                   smaller values there is a loss of accuracy due
c***                   to the non-polynomic behavior of the function
c***                   at small rho. for those values a series expansion
c***                   is employed which makes use of the regular function.
c***
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m7000
      subroutine extcou(f,df,g,dg,a,b,a0,b0,c0,all,x,xinv,break,c,
     1                  ind,lval,energy,n,ntrms,order,nbreak,
     2                  toskp,filnm,nl)
      implicit integer (a-z)
      real*8 f, df, g, dg, a, b, a0, b0, c0, x, xinv, break, c, all
      real*8 dum, energy, rswtch, cl, dl, c0sq, pre, qlpl
      character*2 itoc
      character*4 rowlab, collab(6)
      character *16 fptoc, enchr
      character *18 label
      character *(*) filnm
      dimension f(n), df(n), g(n), dg(n), x(n), xinv(n)
      dimension all(n,5), lval(nl), c(order,nbreak), break(nbreak+1)
      dimension ind(n), a(*), b(*), a0(*), b0(*), c0(*), dum(4)
      common /io/ inp, iout
      enchr=fptoc(energy)
      call iosys('read real "switching r" from atomci',1,rswtch,0,' ')
      call fndswt(x,rswtch,n,iswtch)
      ibeg=iswtch+1
      ndel=n-iswtch
      write(iout,*)
      write(iout,*) '     testing spline routine'
      write(iout,1) energy
      collab(1)='r'
      collab(2)='f'
      collab(3)='df'
      collab(4)='g'
      collab(5)='dg'
      collab(6)='wron'
      call fndbrk(x,break,ind,n,nbreak)
      do 200 l=1,nl
         label=enchr//'-'//itoc(lval(l))
         if (energy.ge.0.d0) then
             filnm='"'//label//'-a"' 
             call iosys ('read real '//filnm//' from atomci',ntrms,
     1                    a,0,' ')
             filnm='"'//label//'-b"' 
             call iosys ('read real '//filnm//' from atomci',ntrms,
     1                    b,0,' ')
             filnm='"'//label//'-cl"' 
             call iosys ('read real '//filnm//' from atomci',1,cl,0,' ')
             filnm='"'//label//'-factors"'
             call iosys ('read real '//filnm//' from atomci',4,
     1                    dum,0,' ')
             dl=dum(1)
             c0sq=dum(2)
             pre=dum(3)
             qlpl=dum(4)
             write(iout,2) dl, c0sq, pre, qlpl
             filnm='"'//label//'-c reg"'
             call iosys ('read real '//filnm//' from atomci',
     1                    nbreak*order,c,0,' ')
             call ppval(x(ibeg),f(ibeg),break,c,ndel,nbreak,order,
     1                  ind(ibeg),0)
             call ppval(x(ibeg),df(ibeg),break,c,ndel,nbreak,order,
     1                  ind(ibeg),1)
             call regexp(f,df,x,a,lval(l),cl,iswtch,ntrms,.false.)
             filnm='"'//label//'-c ireg"' 
             call iosys ('read real '//filnm//' from atomci',
     1                    nbreak*order,c,0,' ')
             call ppval(x(ibeg),g(ibeg),break,c,ndel,nbreak,order,
     1                  ind(ibeg),0)
             call ppval(x(ibeg),dg(ibeg),break,c,ndel,nbreak,order,
     1                  ind(ibeg),1)
             call iregxp(g,dg,f,df,x,b,lval(l),dl,pre,qlpl,iswtch,
     1                   ntrms,.false.)
         else
             filnm='"'//label//'-a0"' 
             call iosys ('read real '//filnm//' from atomci',ntrms,
     1                    a0,0,' ')             
             filnm='"'//label//'-b0"' 
             call iosys ('read real '//filnm//' from atomci',ntrms,
     1                    b0,0,' ')             
             filnm='"'//label//'-c0"' 
             call iosys ('read real '//filnm//' from atomci',ntrms,
     1                    c0,0,' ')             
             filnm='"'//label//'-factors"'
             call iosys ('read real '//filnm//' from atomci',
     1                    1,pre,0,' ')
             filnm='"'//label//'-c reg"'
             call iosys ('read real '//filnm//' from atomci',
     1                    nbreak*order,c,0,' ')
             call ppval(x(ibeg),f(ibeg),break,c,ndel,nbreak,order,
     1                  ind(ibeg),0)
             call ppval(x(ibeg),df(ibeg),break,c,ndel,nbreak,order,
     1                  ind(ibeg),1)
             call regexn(f,df,x,xinv,a0,lval(l),iswtch,ntrms,.false.)
             filnm='"'//label//'-c ireg"' 
             call iosys ('read real '//filnm//' from atomci',
     1                    nbreak*order,c,0,' ')
             call ppval(x(ibeg),g(ibeg),break,c,ndel,nbreak,order,
     1                  ind(ibeg),0)
             call ppval(x(ibeg),dg(ibeg),break,c,ndel,nbreak,order,
     1                  ind(ibeg),1)
             call iregxn(g,dg,a0,b0,c0,pre,x,xinv,lval(l),iswtch,
     1                   ntrms,.false.)
         endif
         do 300 pt=1,n
            all(pt,6)=f(pt)*dg(pt)-g(pt)*df(pt)
  300    continue
         call matpre(all,n,6,n,6,0,1,rowlab,collab,0,dum,.false.)
  200 continue
      return
    1 format(/,15x,'energy = ',e15.8)
    2 format(/,5x,'dl = ',e15.8,1x,'c0sq = ',e15.8,1x,'pre = ',e15.8,/,
     1            20x,'qlpl = ',e15.8)
      end

























