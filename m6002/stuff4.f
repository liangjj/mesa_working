*deck @(#)stuff4.f	1.1 9/7/91
c***begin prologue     stuff4
c***date written                (yymmdd)
c***revision date      890411   (yymmdd)
c***keywords           stuff4, link 6003
c***authors            unknown
c***                   
c***source             m6002
c***purpose            interpolating taylor series for f functions
c***references       
c
c***routines called    
c***end prologue       stuff4
      subroutine stuff4(n,nlt0,nge0,arg,pointr,fvec)
      parameter (maxr=350)
      implicit real *8 (a-h,o-z)
      integer pointr
      dimension arg(n), fvec(1:n,0:6), pointr(n,2)
      common/store/str0(maxr),str1(maxr),str2(maxr),str3(maxr),
     1             str4(maxr),str5(maxr),str6(maxr),str7(maxr),
     2             str8(maxr),str9(maxr),str10(maxr)
c
c  interpolating taylor series for f integrals
c  vectorizable algorithm
c
c
      do 1000 i = 1,nlt0
         iii=pointr(i,1)
         t=arg(iii)
         x= 10.d0*(t+0.05d0)
         iit=x
         ti=iit
         it=iit+1
         delt=t-0.1d0*ti
         delt2=0.5d0*delt
         delt3=-.33333333d0*delt
         delt4=-0.25d0*delt
         ft0=str0(it)
         ft1=str1(it)
         ft2=str2(it)
         ft3=str3(it)
         ft4=str4(it)
         ft5=str5(it)
         ft6=str6(it)
         ft7=str7(it)
         ft8=str8(it)
         fvec(iii,4)=ft4+delt*(-ft5+delt2*(ft6+delt3*(ft7+delt4* ft8)))
         fvec(iii,3)=ft3+delt*(-ft4+delt2*(ft5+delt3*(ft6+delt4* ft7)))
         fvec(iii,2)=ft2+delt*(-ft3+delt2*(ft4+delt3*(ft5+delt4* ft6)))
         fvec(iii,1)=ft1+delt*(-ft2+delt2*(ft3+delt3*(ft4+delt4* ft5)))
         fvec(iii,0)=ft0+delt*(-ft1+delt2*(ft2+delt3*(ft3+delt4* ft4)))
1000  continue
c
      do 2000 i = 1,nge0
         iii=pointr(i,2)
         t=arg(iii)
         xd=1.d0/t
         fvec(iii,0) = .88622692d0* sqrt(xd)
         fvec(iii,1) = 0.5d0*xd*fvec(iii,0)
         fvec(iii,2) = 1.5d0*xd*fvec(iii,1)
         fvec(iii,3) = 2.5d0*xd*fvec(iii,2)
         fvec(iii,4) = 3.5d0*xd*fvec(iii,3)
2000  continue
      return
      end
