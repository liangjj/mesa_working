      subroutine stuff0(n,arg,fvec,n1,n2,iarg1,iarg2)
      parameter (mxbuf=1000)
      implicit real*8 (a-h,o-z)
      dimension arg(mxbuf),fvec(1:mxbuf,0:6),iarg1(mxbuf),iarg2(mxbuf)
      common/store/str0(350),str1(350),str2(350),str3(350),str4(350),
     1str5(350),str6(350),str7(350),str8(350),str9(350),str10(350)
c
c  interpolating taylor series for f integrals
c  vectorizable algorithm
c
      do 1000 iscat = 1,n1
      iii=iarg1(iscat)
      t=arg(iii)
      x= 10.e0*(t+0.05e0)
      iit=x
      ti=iit
      it=iit+1
      delt=t-0.1e0*ti
      delt2=0.5e0*delt
      delt3=-.33333333e0*delt
      delt4=-0.25e0*delt
      ft0=str0(it)
      ft1=str1(it)
      ft2=str2(it)
      ft3=str3(it)
      ft4=str4(it)
 1000 fvec(iii,0)=ft0+delt*(-ft1+delt2*(ft2+delt3*(ft3+delt4* ft4)))
      do 1001 iscat=1,n2
      iii=iarg2(iscat)
      t=arg(iii)
      xd=1.e0/t
      fvec(iii,0) =.88622692e0* sqrt(xd)
 1001 continue
      return
      end
