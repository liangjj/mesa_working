      subroutine gaussval(ibas,maxpts,narg,grid,icn,eta1,l1,m1,n1,v,vpp,
     1maxnbf,maxprim1)
      implicit real*8 (a-h,o-z)
      dimension grid(maxpts,4),eta1(maxprim1,5),l1(maxprim1)
      dimension m1(maxprim1),n1(maxprim1)
      dimension v(maxpts,maxnbf),vpp(maxpts,maxnbf)
      do 500 ii=1,icn
      do 400 k=1,narg
      term = eta1(ii,5)*(grid(k,1)-eta1(ii,1))**l1(ii)
     x *(grid(k,2)-eta1(ii,2))**m1(ii)
     x *(grid(k,3)-eta1(ii,3))**n1(ii)
     x *exp(-eta1(ii,4)*((grid(k,1)-eta1(ii,1))**2+
     x (grid(k,2)-eta1(ii,2))**2 + (grid(k,3)-eta1(ii,3))**2))
      v(k,ibas) = v(k,ibas) + term
      tx = (grid(k,1) - eta1(ii,1))**2
      ty = (grid(k,2) - eta1(ii,2))**2
      tz = (grid(k,3) - eta1(ii,3))**2
      vpp(k,ibas) = vpp(k,ibas) - 0.5*term*(
     $l1(ii)*(l1(ii)-1)/tx-2*eta1(ii,4)*(2*l1(ii)+1)+4.*eta1(ii,4)**2*tx
     $ +
     $m1(ii)*(m1(ii)-1)/ty-2*eta1(ii,4)*(2*m1(ii)+1)+4.*eta1(ii,4)**2*ty
     $ +
     $n1(ii)*(n1(ii)-1)/tz-2*eta1(ii,4)*(2*n1(ii)+1)+4.*eta1(ii,4)**2*tz
     $ )
400   continue
500   continue
      return
      end
