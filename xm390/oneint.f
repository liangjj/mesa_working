*deck @(#)oneint.f	2.1  10/10/91
 
      subroutine oneint(x,n1int,lbli,stvi,lenbuf,nnp,ijpt,symoff,nsym,
     #                  itape,bfsym,bfnum)
c
      implicit integer (a-z)
c
      real*8 x(n1int),stvi(lenbuf)
      integer lbli(lenbuf),ijpt(nnp),symoff(nsym)
      integer bfsym(lenbuf),bfnum(lenbuf)
c
      common /io/ inp,iout
c
      ioff(i,j)=i*(i-1)/2+j
c
c
      call rzero(x,n1int)
 
       read(14) kk
       write(iout,26)
   26 format(/,10x,'next one-electron integral set')
c
  200 continue
c **** comment out pauls reads
c        call sread(itape,lbli,lenbuf)
c        ilsti=lbli(1)
c        nbuf=lbli(2)
 
c **** read in from ijkl tape
 
      read(14) nbuf,ilsti,lbli,stvi
       if(nbuf.le.0) go to 110
 
         do 101 ii=1,nbuf
c ***** comment out pauls indexing and unpack ijkl
 
c           jsm=rshift(lbli(ii+2),8)
c           ior=rshift(jsm,3)
c           ism=rshift(ior,8) + 1
c           ior=and(ior,255)
c           jsm=and(jsm,7) + 1
c           jor=and(lbli(ii+2),255)
            j1=shiftr(lbli(ii),40).a.(.n.mask(54))
            i1=shiftr(lbli(ii),50)
            ism=bfsym(i1)
            jsm=bfsym(j1)
            ior=bfnum(i1)
            jor=bfnum(j1)
            i=symoff(ism)+ior
            j=symoff(jsm)+jor
            ij=ioff(i,j)
            x(ijpt(ij))=stvi(ii)
cps          write(iout,25) i,j,ism,jsm,x(ijpt(ij))
cps 25      format(5x,' i=',i3,' j=',i3,' isym=',i3,' jsym=',i3,/,5x,
cps     #' the int value=',e13.6)
  101    continue
  110 if(ilsti.eq.0) go to 200
c
      return
      end
