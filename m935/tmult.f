*deck @(#)tmult.f	5.1  11/6/94
      subroutine tmult(bptr,nbins,xb,ib,xm,binsiz,
     $                 binsz2,tblks,nblks,lsize,blksiz,
     $                 mxblk,iu,iprt,vec,svec,tvec,
     $                 mdim,npvec,incore,energy,filtyp)
      implicit integer(a-z)
      real*8 xb(binsiz),xm(*),vec(mdim,npvec),
     $       svec(mdim,npvec),tvec(lsize,npvec),energy
      integer ib(binsz2),bptr(nbins)
      character*(*) filtyp
c
c
      ii=1
      jx=0
      nx=1
c
      call rzero(svec,mdim*npvec)
c
      if(incore.eq.0) then
c
c        -- out-of-core matrix multiplication
c
         call rzero(xm,blksiz)
         ptr=bptr(1)
         if (ptr.ne.0) then      
  1         continue
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
  4                  continue
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
               end if
c
               jx=jx+1
c
               if (nz.ne.0) then
                  call rzero(tvec,lsize*npvec)
                  call mxma(xm(ii),1,lsize,vec(jp,1),1,mdim,
     $                      tvec,1,lsize,leni,lenj,npvec)
                  do 11 ia=1,npvec
                     do 12 ja=1,leni
                        svec(ip-1+ja,ia)=svec(ip-1+ja,ia)+
     $                                   tvec(ja,ia)
   12                continue
   11             continue
                  if(i.ne.j) then
                     call rzero(tvec,lsize*npvec)
                     call mxma(xm(ii),lsize,1,vec(ip,1),1,mdim,
     $                         tvec,1,lsize,lenj,leni,npvec)
                     do 13 ia=1,npvec
                        do 14 ja=1,lenj
                           svec(jp-1+ja,ia)=svec(jp-1+ja,ia)+
     $                                      tvec(ja,ia)
   14                   continue
   13                continue
                  end if
               endif
c
c
               ii=ii+lsize*lsize
               jp=jp+lsize
    3       continue
            ip=ip+lsize
    2    continue
c
      else
c
c        -- incore matrix multiplication
c
         ip=1
         do 102 i=1,mxblk
            leni=min(lsize,mdim-ip+1)
            jp=1
            do 103 j=1,i
               lenj=min(lsize,mdim-jp+1)
c
c
               call rzero(tvec,lsize*npvec)
               call mxma(xm(ii),1,lsize,vec(jp,1),1,mdim,
     $                   tvec,1,lsize,leni,lenj,npvec)
               do 111 ia=1,npvec
                  do 112 ja=1,leni
                     svec(ip-1+ja,ia)=svec(ip-1+ja,ia)+tvec(ja,ia)
  112             continue
  111          continue
c
               if(i.ne.j) then
                  call rzero(tvec,lsize*npvec)
                  call mxma(xm(ii),lsize,1,vec(ip,1),1,mdim,t
     $                      vec,1,lsize,lenj,leni,npvec)
                  do 113 ia=1,npvec
                     do 114 ja=1,lenj
                        svec(jp-1+ja,ia)=svec(jp-1+ja,ia)+tvec(ja,ia)
  114                continue
  113             continue
               end if
c
c
               ii=ii+lsize*lsize
               jp=jp+lsize
  103       continue
            ip=ip+lsize
  102    continue
c
      end if
c
      do 201 i=1,npvec
         call saxpy(mdim,-energy,vec(1,i),1,svec(1,i),1)
 201  continue
c
      return
      end
