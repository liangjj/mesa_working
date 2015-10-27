      subroutine vints(l1,l2,m1,m2,n1,n2,icn,jcn,eta1,eta2,grid,narg
     x ,valint,den,soo,istate,maxp,maxprim1)
      parameter (mxbuf=1000)
      implicit real*8 (a-h,o-z)
c
c  computes matrix elements of 1/abs(r1-r2) for points
c on a grid of r2 vector values for a quadrature
c
c l1 and l2 are powers of x for functions 1 and 2
c m1 and m2 are powers of y for functions 1 and 2
c n1 and n2 are powers of z for functions 1 and 2
c icn and jcn are the number of primitives in functions 1 and 2
c eta1 eta2 refer to functions 1 and 2 and contain:
c    x y and z coordinates in locations 1 to 3, orbital exponents in
c    location 4, and contraction coefficients(x normalization) in location 5
c     first index of eta1 and eta2 is the index of the primitive
c grid(i,1 to 3) contains narg r2 cartesian coordinates
c den is the density matrix element which goes with this
c basis function pair
c
      dimension l1(maxprim1),l2(maxprim1),m1(maxprim1)
      dimension m2(maxprim1),n1(maxprim1),n2(maxprim1)
      dimension valint(mxbuf,maxp), eta1(maxprim1,5)
      dimension eta2(maxprim1,5), g(mxbuf,7,3)
      dimension grid(mxbuf,3),pcx(mxbuf),pcy(mxbuf),pcz(mxbuf)
      dimension fvec(mxbuf,7),iarg1(mxbuf),iarg2(mxbuf)
      dimension den(maxp),soo(maxp),arg(mxbuf)
c
      dimension vrint(mxbuf)
c
      common /nmbrs /  pi, piterm, pitern, acrcy, scale, icanon
c
      do 635 ii=1,icn
         a=eta1(ii,4)
         do 1635 jj=1,jcn
c..   bhl
            do 6635 ipi=1,mxbuf
               vrint(ipi)=0.0
 6635       continue
c..bhl 
            mx=l1(ii)+l2(jj)+1
            my=m1(ii)+m2(jj)+1
            mz=n1(ii)+n2(jj)+1
            b=eta2(jj,4)
            t1=a+b
            t=1.e0/t1
            p1=(a*eta1(ii,1)+b*eta2(jj,1))*t
            p2=(a*eta1(ii,2)+b*eta2(jj,2))*t
            p3=(a*eta1(ii,3)+b*eta2(jj,3))*t
            ab1=eta1(ii,1)-eta2(jj,1)
            ab2=eta1(ii,2)-eta2(jj,2)
            ab3=eta1(ii,3)-eta2(jj,3)
            distab=ab1*ab1+ab2*ab2+ab3*ab3
            do 1 i=1,istate
               soo(i)=den(i)*2.e0*pi*t* 
     &              exp(-a*b*distab*t)*eta1(ii,5)*eta2(jj,5)
c               write(66,"('den(i)=',f10.4,'   ai=',f10.4,'   aj=',f10.4,
c     &              '   ci=',f10.4,'   cj=',f10.4,'   soo=',f10.4)")
c     &              den(i),
c     &              eta1(ii,4),eta2(jj,4),
c     &              eta1(ii,5),eta2(jj,5),soo(i)/pi/2/t
 1          continue
c     
c  grid contains a block of 500 or fewer grid points
c
c compute the arguments of the f functions
            do 690 ig=1,narg
               pcx(ig)=p1-grid(ig,1)
               pcy(ig)=p2-grid(ig,2)
               pcz(ig)=p3-grid(ig,3)
               pcsq=pcx(ig)*pcx(ig)+pcy(ig)*pcy(ig)+pcz(ig)*pcz(ig)
               arg(ig)=t1*pcsq
 690        continue
c
c preprocess arg vector to determine appropriate way to compute integrals
c
            i1=0
            i2=0
            do 1000 i=1,narg
               switch=arg(i)-34.9e0
               if(switch.le.0.0)then
                  i1=i1+1
                  iarg1(i1)=i
               else
                  i2=i2+1
                  iarg2(i2)=i
               endif
 1000       continue
            narg1=i1
            narg2=i2
            if((narg1+narg2).ne.narg)then
               write(6,*)' vints sort inconsistency'
               stop
            endif
c
c compute the f functions
c
            maxtyp=mx+my+mz-2
c     write(66,895) maxtyp
 895        format(' fint type',i3)
            go to (691,692,693,694,695,696,697) maxtyp
 691        call stuff0(narg,arg,fvec,narg1,narg2,iarg1,iarg2)
            go to 698
 692        call stuff1(narg,arg,fvec,narg1,narg2,iarg1,iarg2)
            go to 698
 693        call stuff2(narg,arg,fvec,narg1,narg2,iarg1,iarg2)
            go to 698
 694        call stuff3(narg,arg,fvec,narg1,narg2,iarg1,iarg2)
            go to 698
 695        call stuff4(narg,arg,fvec,narg1,narg2,iarg1,iarg2)
            go to 698
 696        call stuff5(narg,arg,fvec,narg1,narg2,iarg1,iarg2)
            go to 698
 697        call stuff6(narg,arg,fvec,narg1,narg2,iarg1,iarg2)
 698        continue
c
            pax=p1-eta1(ii,1)
            pbx=p1-eta2(jj,1)
            call  gfunct(l1(ii),l2(jj),pax,pbx,pcx,t,g,1,narg)
            pay=p2-eta1(ii,2)
            pby=p2-eta2(jj,2)
            call  gfunct(m1(ii),m2(jj),pay,pby,pcy,t,g,2,narg)
            paz=p3-eta1(ii,3)
            pbz=p3-eta2(jj,3)
            call  gfunct(n1(ii),n2(jj),paz,pbz,pcz,t,g,3,narg)
            do 506 ix=1,mx
               do 507 jy=1,my
                  do 508 kz=1,mz
                     mxyz=ix+jy+kz-2
c..bhl
c     do 508 is=1,istate
c..bhl
                     do 5081 ig=1,narg
c..bhl
c     508 valint(ig,is)=valint(ig,is)+g(ig,ix,1)*g(ig,jy,2)*g(ig,kz,3)
c     x *fvec(ig,mxyz)*soo(is)
c..   bhl
                        vrint(ig)=vrint(ig)+
     X                       g(ig,ix,1)*g(ig,jy,2)*g(ig,kz,3)
     x                       *fvec(ig,mxyz)
c            if(l1(ii)+m1(ii)+n1(ii).eq.2.and.
c     &                       l2(jj)+m2(jj)+m2(jj).eq.2.and.
c     &                       eta1(ii,3)*eta2(jj,3) .lt. 0)then
c               write(66,*)"ix,iy,iz,jx,jy,jz,xv,yy,zv",
c     &              ix,iy,iz,jx,jy,jz,xv,yy,zv 
c               write(66,*)"g(ig,ix,1),g(ig,jy,2),g(ig,kz,3),
c     &              fvec(ig,mxyz)",
c     &              g(ig,ix,1),g(ig,jy,2),g(ig,kz,3),
c     &              fvec(ig,mxyz)
c            endif
 5081                continue
 508              continue
 507           continue
 506        continue
c     close contraction loops
c..   bhl
            do 55081 is=1,istate
               do 5508 ig=1,narg
                  valint(ig,is)=valint(ig,is)+vrint(ig)*soo(is)
c                  write(66,"('vint=', f12.8)")valint(ig,is)
 5508          continue
55081       continue
c..   bhl
c
 1635    continue
 635  continue
      do 55091 is=1,istate
         do 5509 ig=1,narg
c            write(66,*) is,ig,valint(ig,is)
 5509    continue
55091 continue
      return
      end
