      subroutine inputs(nbf,lnew,mnew,nnew,etanew,nfirst,nlast,
     x nnuc,cent,maxnbf,mxcen,maxprim)
c
c read basis function input from file made by m950 in mesa
c
c
c  the way this works is that there are nbf contracted functions
c
c   for the ith contracted function nfirst(i) and nlast(i) are
c  the first and last locations for the relevant quantities in
c   the arrays lnew(ii), mnew(ii), nnew(ii), eta(ii, jj=1,5)
c   for that contracted function
c
c lnew, mnew, nnew, are the powers of x y and z and eta(ii,jj=1,5)
c contains center, exponent and contraction coefficient
c
      implicit real*8 (a-h,o-z)
      dimension lnew(maxprim),mnew(maxprim),nnew(maxprim)
      dimension etanew(maxprim,5)
      dimension nfirst(maxnbf),nlast(maxnbf)
      dimension cent(4,mxcen)
      intp=5
c
      read(8) nnuc
      if (nnuc .gt. mxcen) then
         write(6,*) 'Error nnuc > mxcen'
         stop 'bad nnuc'
      end if
      do 10 i=1,nnuc
      read(8) (cent(jjj,i),jjj=1,4)
  10  continue
      read(8) nbf
      if (nbf .gt. maxnbf) then
         write(6,*) 'Error too many functions in xustatnuc',
     $         nbf, maxnbf
         stop
      end if
      read(8) (nfirst(i),nlast(i),i=1,nbf)
      do 11 i=1,nbf
      ilow=nfirst(i)
      ihi=nlast(i)
      if (ihi .gt. maxprim) then
         write(6,*) 'Error too many primitives in xustatnuc',
     $         ihn, maxprim
         stop
      end if

      do 12 j=ilow,ihi
      read(8) lnew(j),mnew(j),nnew(j),(etanew(j,jj),jj=1,5)
c avoid 0 to the power 0 later in gausval
      etanew(j,1)=etanew(j,1)+1.e-20
      etanew(j,2)=etanew(j,2)+1.e-20
      etanew(j,3)=etanew(j,3)+1.e-20
  12  continue
  11  continue
c
c  write out basis function info
c
      write(66,200)
200   format(///,' basis functions from mesa interface')
      write(66,201)(inuc,(cent(jjj,inuc),jjj=1,4),inuc=1,nnuc)
201   format(' center ',i3,3x,4f10.5)
      do 500 i=1,nbf
      istrt=nfirst(i)
      il=nlast(i)
      write(66,101) i, istrt,il
101   format(//,' basis fcn ',i3,' nfirst and nlast are ',2i4)
      do 400 ii=istrt,il
      write(66,102) lnew(ii),mnew(ii),nnew(ii),(etanew(ii,jj),jj=1,5)
102   format(1x,3i3,5f10.5)
400   continue
500   continue
      return
 900  format(10a8)
 901  format( 1h1///10x,10a8)
 902  format(20i3)
 903  format(20a4)
 904  format(a8,2x,4f12.0)
 905  format(i5,5x,a8,3f10.6,f10.1)
 906  format(3i3,5x,20a6)
 907   format(a8,2x,a8,2x,f10.0)
 908   format(2i5,5x,2a8,2f20.5)
 909  format(5e20.10)
 910  format(2i5,5x,a3)
 911  format(10e12.5)
 998  format(9htype name,a8,31hdoes not match basis fntn. type,a8)
 999  format(11hcenter name,a8,31hdoes not match basis fntn. name,a8)
      end
