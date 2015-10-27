*deck drvh2.f
c***begin prologue     drvh2
c***date written       000710   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           one-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            hamiltonian elements  for two dimensional hamiltonian. 
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       drvh2
      subroutine drvh2(pham,phone,phamil,n,len,nonz,ngot,
     1                 pack,drctv,type,ondisk,sym,dim,prn)
      implicit integer (a-z)
#ifdef DECPOINTER
      integer*8 pham, phone, phamil
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pham, phone, phamil
#endif SGIPOINTER
      real*8 h, hx, hy, v, hamil, scr
      logical incore, pack, ondisk, prn
      character*(*) sym, drctv, type
      character*8 key
      dimension pham(dim), phone(dim)
      dimension len(dim+1), nonz(dim+1), n(dim+1)
      dimension ngot(*), wrds(5)
      dimension prn(*) 
      dimension hbuf(4), buf(4), diag(4)
      common/io/inp, iout
      pointer (phx,hx(1))
      pointer (phx,ihx(1))
      pointer (phy,hy(1))
      pointer (phy,ihy(1))
      pointer (ph,h(1))
      pointer (ph,ih(1))
      pointer (pvt,v(1))
      pointer (phamil,hamil(1))
      pointer (pind,index(1))
      pointer (pscr,scr(1))
c
      write(iout,1) sym
      need=wptoin(n(dim+1))
      call getmem(need,pvt,wrds(1),'vint',0)
      call rzero(v,n(dim+1))
      if(type.ne.'h0') then
         key='$vpert'
         phx=pham(1)
         phy=pham(2)
         vx=1+n(1)*n(1)
         vy=1+n(2)*n(2)
         x=vx+n(1)+n(1)*n(1)+3
         y=vy+n(2)+n(2)*n(2)+3                      
         call vpert(v,hx(x),hy(y),scr,hx(vx),hy(vy),scr,n,
     1              dim,key,wrds(1),.true.,prn(1))     
      endif
      call getmem(n(1)*n(2),pind,wrds(2),'index',0)
      call setp2(index,n,sym,prn(2))
c
c        get the memory for the buffered two dimensional hamiltonian
c
      hbuf(dim+1)=1
      buf(dim+1)=wpadti(hbuf(dim+1)+len(dim+1))
      diag(dim+1)=iadtwp(buf(dim+1)+2*len(dim+1))                
      words=diag(dim+1)+n(dim+1)
      if(sym.eq.'symmetric') then
         one=words
         words=one+n(1)
      endif      
      need=wpadti(words)
      call getmem(need,ph,wrds(3),'hambuf',0)
      incore=.true.
      if (len(dim+1).le.n(dim+1).or.ondisk) then
          incore=.false.
          call iosys('create integer "hamiltonian buffers" on ham',
     1               -1,0,0,' ')
      endif
      if(incore) then
         write(iout,2)
      else
         write(iout,3)
      endif	  	 
c
c     set pointers
c
      if(pack) then
         phx=phone(1)            
         phy=phone(2)
c
c     set the array pointers for the two one dimensional arrays
c
         hbuf(1)=1
         buf(1)=wpadti(hbuf(1)+len(1))
         diag(1)=iadtwp(buf(1)+2*len(1))
         if(sym.eq.'unsymmetric') then
            hbuf(2)=1
            buf(2)=wpadti(hbuf(2)+len(2))
            diag(2)=iadtwp(buf(2)+2*len(2))
         endif
         call h2pac(h(hbuf(dim+1)),ih(buf(dim+1)),h(diag(dim+1)),
     1              v,index,hx(hbuf(1)),ihx(buf(1)),hx(diag(1)),
     2              hy(hbuf(2)),ihy(buf(2)),hy(diag(2)),dim,
     3              len,nonz,n,incore)
         size=n(dim+1)
      else
         phx=pham(1)
         phy=pham(2)
         if(sym.eq.'unsymmetric') then
            call ham2d(h(hbuf(dim+1)),ih(buf(dim+1)),h(diag(dim+1)),
     1                 v,index,hx,hy,dim,len(dim+1),nonz(dim+1),n,
     2                 incore)
            size=n(dim+1)
         elseif(sym.eq.'symmetric') then
            call sym2d(h(hbuf(dim+1)),ih(buf(dim+1)),h(diag(dim+1)),
     1                 v,index,hx,dim,len(dim+1),nonz(dim+1),n,
     2                 incore)
            size=n(1)*(n(1)+1)/2
         endif          
      endif
c
c     release memory for interaction potential and 2d-packing index.
c
      call getmem(-wrds(1),pvt,idum,'vint',idum)
      call getmem(-wrds(2),pind,idum,'index',idum)
      write(iout,4) nonz(dim+1)
      if(prn(3)) then
         call rdham(h(hbuf(dim+1)),ih(buf(dim+1)),h(diag(dim+1)),
     1              len(dim+1),nonz(dim+1),incore,size)
      endif
c
      if(drctv.eq.'form') then
         htot=1
         need=wpadti(htot+size*size)
         call getmem(need,phamil,words,'hamil',0)
      elseif(drctv.eq.'diagonalize') then
         htot=1
         eig=htot+size*size
         need=wpadti(eig+size)
         call getmem(need,phamil,words,'hamil',0)
         scrwds=wptoin(5*size)
         call getmem(scrwds,pscr,got,'scr',0)
      else
         call getmem(-wrds(3),ph,idum,'hambuf',idum)
         return
      endif
      call diagh(hamil(htot),hamil(eig),scr,h(hbuf(dim+1)),
     1           ih(buf(dim+1)),h(diag(dim+1)),len(dim+1),
     2           nonz(dim+1),incore,drctv,size)
      if(drctv.eq.'diagonalize') then
         call getmem(-got,pscr,idum,'scr',idum)
      endif
      call getmem(-wrds(3),ph,idum,'hambuf',idum)
      return
 1    format(/,1x,'constructing hamiltonian for a ',a16,
     1            ' hamiltonian')
 2    format(/,1x,'hamiltonian stored in memory')
 3    format(/,1x,'hamiltonian stored on disk')
 4    format(/,1x,'number of non-zero matrix elements = ',i10)
      end       










