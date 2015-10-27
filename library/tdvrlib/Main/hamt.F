*deck hamt.f
c***begin prologue     hamt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           1-d hamiltonian
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate one body time matrix.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       hamt
      subroutine hamt(phamt,pvt,key,typ,n,ngot,reuse,prn)
      implicit integer (a-z)
#ifdef DECPOINTER
      integer*8 phamt, pvt
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 phamt, pvt
#endif SGIPOINTER
      real*8 ht, vt, hbar
      character*(*) key, typ 
      character*80 loctyp
      character*80 title, cpass, units, chrkey
      character*240 card
      logical prn, dollar, reuse
      common/io/inp, iout
      save loctyp, units
      data hbar/1.054571596d-34/                                  
      pointer(phamt,ht(1))
      pointer(pvt,vt(1))
      need=wptoin(n)
      call getmem(need,pvt,ngot,'v0t',0)
      h=1
      t0=h+n*n
      t=t0+1
      call rzero(vt,n) 
      if(reuse) then
         typ=loctyp
      else
         if(dollar(key,card,cpass,inp)) then
            typ=chrkey(card,'potential','none',' ')
            loctyp=typ
            units=chrkey(card,'units','atomic-units',' ') 
         endif
      endif
      call vt0(card,ht(t),vt,hbar,units,n,typ)
      call vscale(ht(h),ht(h),hbar,n)
      if(prn) then
         title='time-dependent one body potential'
         call prntrm(title,vt,n,1,n,1,iout)
         title='time-dependent one body hamiltonian'
         call prntrm(title,ht(h),n,n,n,n,iout)
      endif
      return         
      end       




