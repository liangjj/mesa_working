*deck %W%  %G%
      subroutine davsdr(cr,icr,nwint,buf,
     $     thrv,thrn,thrs,
     $     n,nmax,nroot,nsml,nlook,msml)
c
c***begin prologue     davsdr
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       davsdr
c
      implicit real*8(a-h,o-z)
      integer wpadti
c
      dimension cr(2),icr(2),buf(2),thrv(20)
c
      common /io/ inp,iout
c
c...
c...  equivalence(cr,icr)
c...  integer icr(nwint)
c...
      mdim=max(nsml,nmax*msml+msml)
      itr=(mdim*(mdim+1))/2
      isqr=mdim*mdim
      iblk=nmax*n*msml+msml*n
c
      ib=1
      ic=ib+iblk
      ish=ic+iblk+msml*n
      idiag=ish+itr
      itemp=idiag+n
      ieval=itemp+itr
      ievec=ieval+mdim
      ipt=wpadti(ievec+isqr)
      iscr=iadtwp(ipt+nsml)
      need=wpadti(iscr+5*n)
c
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(need,icr,ncore,'davidson driver',0)
c
      call davson(cr(ic),cr(ib),cr(ish),
     $     cr(idiag),cr(itemp),cr(ieval),
     $     cr(ievec),buf,icr(ipt),
     $     cr(iscr),thrv,thrn,thrs,
     $     n,nmax,nroot,nsml,nlook,msml)
c
      return
      end
