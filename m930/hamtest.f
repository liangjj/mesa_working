*deck  @(#)test.f	5.1 11/6/94
      subroutine hamtest(bptr,nbins,xb,ib,xm,binsiz,
     $                   binsz2,tblks,nblks,lsize,blksiz,
     $                   mxblk,iu,iprt,vec,svec,tvec,mdim,npvec,filtyp)
      implicit integer(a-z)
      character*(*) filtyp
      real*8 xb(binsiz),xm(*),vec(mdim,npvec),sdot,
     $       svec(mdim,npvec),tvec(lsize,npvec),energy
      integer ib(binsz2),bptr(nbins)
      common /io/ inp,iout
c
c
      ii=1
      jx=0
      nx=1
c
      call rzero(svec,mdim*npvec)
      call rzero(xm,blksiz)
c
      ptr=bptr(1)
      if (ptr.ne.0) then
    1    continue
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
      endif
c
      ip=1
      do 2 i=1,mxblk
         leni=min(lsize,mdim-ip+1)
         jp=1
         do 3 j=1,i
            lenj=min(lsize,mdim-jp+1)
            if(jx.eq.nblks) then
               ii=1
               jx=0
               nx=nx+1
               call rzero(xm,blksiz)
               ptr=bptr(nx)
               if (ptr.ne.0) then
  4               continue
                  call getham(iu,xb,wptoin(binsiz),ptr,filtyp)
                  ptr=ptr+wptoin(binsiz)
                  call getham(iu,ib,binsz2,ptr,filtyp)
                  lbin=ib(binsiz+1)
                  call bscattr(lbin,xm,ib,xb)
                  ptr=ib(binsz2)
                  if(ptr.ne.0) go to 4
                  nz=1
               else
                  nz=0 
               end if
            endif
c
            jx=jx+1
c
            if (nz.ne.0) then
               call rzero(tvec,lsize*npvec)
               call mxma(xm(ii),1,lsize,vec(jp,1),1,mdim,tvec,1,lsize,
     $                   leni,lenj,npvec)
               do 11 ia=1,npvec
                  do 12 ja=1,leni
                     svec(ip-1+ja,ia)=svec(ip-1+ja,ia)+tvec(ja,ia)
   12             continue
   11          continue
c
               if(i.ne.j) then
                  call rzero(tvec,lsize*npvec)
                  call mxma(xm(ii),lsize,1,vec(ip,1),1,mdim,tvec,1,
     $                      lsize,lenj,leni,npvec)
                  do 13 ia=1,npvec
                     do 14 ja=1,lenj
                        svec(jp-1+ja,ia)=svec(jp-1+ja,ia)+tvec(ja,ia)
   14                continue
   13             continue
               end if
            endif
c
            if(iprt.gt.0) then
               write(iout,*)' '
               write(iout,*)' hamiltonian block i j ',i,j
               call matout(xm(ii),lsize,lsize,lsize,lsize,6)
            end if
c 
            ii=ii+lsize*lsize
            jp=jp+lsize
    3    continue
         ip=ip+lsize
    2 continue
c
      energy=sdot(mdim,vec,1,svec,1)
c
      write(iout,201) energy
 201  format(5x,' hamiltonian test completed',/,
     #            10x,'energy = ',f20.10)
c
      return
      end
