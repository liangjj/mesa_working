*deck pac2.f
c***begin prologue     pac2
c***date written       000710   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           one-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            non zero hamiltonian elements and indices for 
c***                   two dimensional hamiltonian. 
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       pac2
      subroutine pac2(ph,index,pham,px,py,pv,len,lenx,leny,nonz,
     1                nonzx,nonzy,n,nx,ny,ngot,prn,sym,region)
      implicit integer (a-z)
      real*8 h, hx, hy, v, ham, scr
      logical incore, prn
      character*(*) sym, region
      dimension prn(*), ngot(2), index(*)
      common/io/inp, iout
      pointer (ph,h(1))
      pointer (ph,ih(1))
      pointer (px,hx(1))
      pointer (px,ihx(1))
      pointer (py,hy(1))
      pointer (py,ihy(1))
      pointer (pv,v(1))
      pointer (pham,ham(1))
      pointer (pscr,scr(1))
c
      write(iout,1) sym
c
c     set the array pointers for the two one dimensional arrays
c
      hbufx=1
      bufx=wpadti(hbufx+lenx)
      diagx=iadtwp(bufx+2*lenx)
      if(sym.eq.'unsymmetric') then
         hbufy=1
         bufy=wpadti(hbufy+leny)
         diagy=iadtwp(bufy+2*leny)
      endif
c
c     get the memory for the buffered two dimensional hamiltonian
c
      hbuf=1
      buf=wpadti(hbuf+len)
      diag=iadtwp(buf+2*len)                
      words=diag+n
      if(sym.eq.'symmetric') then
         one=words
         words=one+nx
      endif      
      need=wpadti(words)
      call memory(need,ph,ngot(1),'hambuf',0)
      incore=.true.
      if (len.le.n) then
          incore=.false.
          call iosys('create integer "hamiltonian buffers" on ham',
     1               -1,0,0,' ')
      endif
      call hamfl2(h(hbuf),ih(buf),h(diag),v,index,len,nonz, 
     1            hx(hbufx),ihx(bufx),hx(diagx),lenx,nonzx,
     2            hy(hbufy),ihy(bufy),hy(diagy),leny,nonzy,
     3            nx,ny,n,incore)
      write(iout,2) nonz
      if(prn(1)) then
         call rdham(h(hbuf),ih(buf),h(diag),len,nonz,incore,n)
      endif
      hamil=1
      eig=hamil+n*n
      need=wpadti(eig+n)
      call memory(need,pham,ngot(2),'ham',0)
      need=wpadti(1+5*n)
      call memory(need,pscr,words,'scr',0)
      call diagh(ham(hamil),ham(eig),scr,h(hbuf),ih(buf),h(diag),
     1           len,nonz,incore,n,region)
      call memory(-words,pscr,idum,'scr',idum)
      return
 1    format(/,1x,'constructing hamiltonian for a ',a16,
     1            ' hamiltonian')
 2    format(/,1x,'number of non-zero matrix elements = ',i5)
      end       










