*deck @(#)rdham.f	5.1  11/6/94
      subroutine rdham(bptr,nbins,xb,ib,xm,binsiz,binsz2,
     $                 tblks,nblks,lsize,blksiz,mxblk,iu,iprt,
     $                 diag,mdim,energy,filtyp)
      implicit integer(a-z)
      real*8 xb(binsiz),xm(*),diag(mdim),energy
      integer ib(binsz2),bptr(nbins)
      character*(*) filtyp
      common /io/ inp,iout
c
c
c
      ii=1
      jx=0
      nx=1
c
c     ----- read the blocked hamiltonian into core and save diagonals -----
c
       call rzero(xm,blksiz)
       ptr=bptr(1)
       if(ptr.ne.0) then
  1       continue
          call getham(iu,xb,wptoin(binsiz),ptr,filtyp)
          ptr=ptr+wptoin(binsiz)
          call getham(iu,ib,binsz2,ptr,filtyp)
          lbin=ib(binsiz+1)
          call bscattr(lbin,xm,ib,xb)
          ptr=ib(binsz2)
          if(ptr.ne.0) go to 1
          nz=1
       else
          nz=0
       end if
c
       ic=0
       ip=1
       do 2 i=1,mxblk
          leni=min(lsize,mdim-ip+1)
          jp=1
          do 3 j=1,i
             lenj=min(lsize,mdim-jp+1)
             if(jx.eq.nblks) then
                jx=0
                nx=nx+1
                call rzero(xm(ii),blksiz)
                ptr=bptr(nx)
                if(ptr.ne.0) then
  4                continue
                   call getham(iu,xb,wptoin(binsiz),ptr,filtyp)
                   ptr=ptr+wptoin(binsiz)
                   call getham(iu,ib,binsz2,ptr,filtyp)
                   lbin=ib(binsiz+1)
                   call bscattr(lbin,xm(ii),ib,xb)
                   ptr=ib(binsz2)
                   if(ptr.ne.0) go to 4
                   nz=1
                else
                   nz=0
                end if
             end if
c
             jx=jx+1
c
             if(nz.ne.0) then
                if(i.eq.j) then
                   kx=1
                   do 5 k=1,leni
                      ic=ic+1
                      diag(ic)=xm(ii-1+kx)-energy
                      kx=kx+lsize+1
   5               continue
                end if
c
                if(iprt.gt.0) then
                   write(iout,*)' '
                   write(iout,*)' hamiltonian block i j ',i,j
                   call matout(xm(ii),lsize,lsize,lsize,lsize,iout)
                end if
             end if
c
             ii=ii+lsize*lsize
             jp=jp+lsize
   3       continue
           ip=ip+lsize
   2   continue
c
      return
      end
