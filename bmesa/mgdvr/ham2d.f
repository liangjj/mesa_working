*deck ham2d.f
c***begin prologue     ham2d
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
c***end prologue       ham2d
      subroutine ham2d(ph,index,pv,px,py,pham,len,nonz,n,nx,ny,
     1                 ngot,prn,sym,region)
      implicit integer (a-z)
      real*8 h, hx, hy, v, ham, scr, ngot(2)
      logical incore, prn
      character*(*) sym, region
      dimension prn(*), index(*)
      common/io/inp, iout
      pointer (ph,h(1))
      pointer (ph,ih(1))
      pointer (pv,v(1))
      pointer (px,hx(1))
      pointer (py,hy(1))
      pointer (pham,ham(1))
      pointer (pscr,scr(1))
c
      write(iout,1) sym
c
c     get the memory for the buffered two dimensional hamiltonian
c
      hbuf=1
      buf=wpadti(hbuf+len)
      diag=iadtwp(buf+2*len)                
      words=diag+n
      need=wpadti(words)
      call memory(need,ph,ngot(1),'hambuf',0)
      incore=.true.
      if (len.le.n) then
          incore=.false.
          call iosys('create integer "hamiltonian buffers" on ham',
     1               -1,0,0,' ')
      endif
      if(sym.eq.'unsymmetric') then
         call conh2(h(hbuf),ih(buf),h(diag),v,index,len,nonz,
     1              hx,hy,nx,ny,n,incore)
         size=n
      elseif(sym.eq.'symmetric') then      
         call symh2(h(hbuf),ih(buf),h(diag),v,index,len,nonz, 
     1              hx,nx,n,incore)
         size=nx*(nx+1)/2
      endif
      write(iout,2) nonz
      if(prn(1)) then
         write(iout,*) 'calling rdham'
         call rdham(h(hbuf),ih(buf),h(diag),len,nonz,incore,size)
      endif
      hamil=1
      eig=hamil+n*n
      need=wpadti(eig+n)
      call memory(need,pham,ngot(2),'ham',0)
      need=wpadti(1+5*n)
      call memory(need,pscr,words,'scr',0)
      call diagh(ham(hamil),ham(eig),scr,h(hbuf),ih(buf),h(diag),
     1           len,nonz,incore,size,region)
      call memory(-words,pscr,idum,'scr',idum)
      return
 1    format(/,1x,'constructing hamiltonian for a ',a16,
     1            ' hamiltonian')
 2    format(/,1x,'number of non-zero matrix elements = ',i5)
      end       










