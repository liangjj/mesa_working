*deck %W% %G%
      subroutine genacint(iatom,jatom,imax,jmax,nprimi,nprimj,
     $     alpha,ainv,inarry,
     $     outarry,tmparry,tmpsiz,t,nv,mmax,pi252,lval,
     $     mycart,yn0,amaxm,mindx,mini,minj,nx,ny,nz)
c***begin prologue     %M%
c***date written       930208 (yymmdd)
c***revision date      %G%
c
c***keywords
c***author             russo, thomas (lanl)
c***source             %W% %G%
c***purpose
c     Generate two-center, two-electron integrals of form (A,C)
c     using the (a,0) integrals previously generated
c***description
c    uses basic obara-saika recurrences, taking advantage of the two-center
c    special case and using translational invariance to build the (a,c+1)
c    recurrence
c     *NOTE*  THIS ROUTINE WILL ONLY WORK WHEN imax>= jmax!  Twoel.f has
c     been written to assure this.  If you change that, you need to change
c     this!  The reason for the restriction is that the integrals required
c     when imax>=jmax follow a simpler rule than when imax<jmax.  This
c     restriction should have no effect on the rest of the code, as we
c     only store the lower triangle of the matrix of integrals anyway, and
c     putel will take care of putting things in the right place.  
c***references
c     obara and saika, j.chem.phys., 84,3963 (1986)
c     
c***routines called
c
c***end prologue     %M%
      implicit integer(a-z)
      logical debug
      real*8 half,one,two
      parameter (debug=.false.)
      parameter (half=0.5d0,one=1.0d0,two=2.0d0)
      common /io/ inp,iout
c
c inarry has the (a,0) integrals, tmparry is where we'll build the 
c (a-l',l') guys.
c
      real*8 inarry(nv,*),outarry(nv,*),tmparry(nv,tmpsiz,0:2),t(nv),
     $     jnk
      integer nx(*),ny(*),nz(*)
      integer oc(3),oc2(3),oa1(3),oa2(3)
      integer lval(0:amaxm,mindx,3),yn0(0:amaxm),mycart(0:amaxm)
      real*8 alpha(nv,2),ainv(nv,5),sclri
c
c statement function to return offset of the (lx,ly,lz) triple within my
c ordering.  In using it to find (lx,ly,lz|cx,cy,cz) one should note that
c for each (lx,ly,lz) there are mycart(cx+cy+cz) guys, so to find the offset 
c of (lx,ly,lz|cx,cy,cz) from (l,0,0|c,0,0) you need
c     intloc(lx,ly,lz)*mycart(c)+intloc(cx,cy,cz)
c which will be 0 if trying to find (l,0,0|c,0,0)
c
      intloc(lx,ly,lz)=((ly+lz)*(ly+lz+1))/2+lz

      last=0
      this=1
      next=2
c
c load up temp array with (a,0)'s, only need to load evens if a+c=even,
c odds if a+c is odd.
c (either start at 1 or 2, i.e. (0,0) or (1,0)
c
      a0off=1+mod(mmax,2)
      toff=1
      if (jmax .ne. 0) then
         lmin=mod(mmax,2)
      else
         a0off=1
         do 666 l=0,mmax-1
            a0off=a0off+mycart(l)
 666     continue 
         lmin=mmax
      endif
      do 10 l=lmin,mmax,2
         call vmove(tmparry(1,toff,this),inarry(1,a0off),nv*mycart(l))
         if (debug) then
            write(iout,*) "check: the (",l,",0) ints are:"
            write(iout,*) (tmparry(1,toff+i,this),i=0,mycart(l)-1)
            write(iout,*)"------------"
         endif
         a0off=a0off+mycart(l)+mycart(l+1)
         toff=toff+mycart(l)
 10   continue 
c
c form  the alphaA/alphaC -- ainv(1,1)=1/alphaA,ainv(1,2)=1/alphaC,
c ainv(1,3)becomes -alphaA/alphaC, other elements unused here.
c
      if (debug) then
         write(iout,*) "Now forming the ainv(1,3)"
         write(iout,*) "Input (alpha(.,1) and ainv(.,2)):"
         write(iout,*) (alpha(i,1),i=1,nv)
         write(iout,*) (ainv(i,2),i=1,nv)
      endif
      call vmul(ainv(1,3),alpha(1,1),ainv(1,2),nv)
      jnk=-1.0
      call smul(ainv(1,3),ainv(1,3),jnk,nv)
      if (debug) then
         write(iout,*) "Output (ainv(.,3)=-alpha(.,1)*ainv(.,2))"
         write(iout,*) (ainv(i,3),i=1,nv)
      endif
c
c------
c     form (la,lc+1) from previously available data
c
c     at each step, tamin is the minimum la in the previous generation,
c     tamax is the maximum, namin=minimum we need to generate, namax=max to gen
c tam1off=offset  of (la-1,lc), tap1off is offset of (la+1,lc)
c laoff=offset (la,lc-1), naoff is where to stick the new ones.  All offsets
c are relative to the beginning of the array, element 1.
c
      tamin=lmin
      tamax=mmax
c 
c generate each of the (la,lc+1) that we'll need for the next step
c start by looping over each triple with momentum la, then form each triple
c with lc+1 and that bra.
c
      do 20 lc=0,jmax-1
         if (debug) then
            write(iout,*)"last is ",last," This is ",this, " Next=",
     $           next
         endif
         namin=tamin+1
         namax=tamax-1
         tam1off=0
         tap1off=mycart(tamin)*mycart(lc)
         if (tamin .gt. 0 .and. lc .gt. 0) then
            laoff=mycart(tamin-1)*mycart(lc-1)
         else
            laoff=0
         endif
         naoff=0
         do 30 la=namin,namax,2
            do 40 thebra=1,mycart(la)
               ax=lval(la,thebra,1)
               ay=lval(la,thebra,2)
               az=lval(la,thebra,3)
               ail=intloc(ax,ay,az)
               bumpaux=0
               if (ax .eq. 0) bumpaux=1
               if (ay .eq. 0) bumpaux=2
               if (az .eq. 0) bumpaux=3
               do 50 theket=1,mycart(lc+1)
                  cx=lval(lc+1,theket,1)
                  cy=lval(lc+1,theket,2)
                  cz=lval(lc+1,theket,3)
c
c now, how what component of what integral should we bump up to get this?
c welp, if any of the cx,cy,cz is 1, do that one, coz the recurrence is
c cheaper.  If not, check to see if any of the ax,ay,az is 0 and bump that
c one for the same reason.  If there is no component on either side which
c is pickable by this reasoning, do any nonzero one, and there should
c be one of those!
c 
                  bumpup=0
                  if (cx .eq. 1) bumpup=1
                  if (cy .eq. 1) bumpup=2
                  if (cz .eq. 1) bumpup=3
                  if (bumpup .eq. 0) then
                     if (bumpaux .eq. 1 .and. cx .ne. 0) bumpup=1
                     if (bumpaux .eq. 2 .and. cy .ne. 0) bumpup=2
                     if (bumpaux .eq. 3 .and. cz .ne. 0) bumpup=3
                     if ( bumpup .eq. 0) then
                        if (cx .ne. 0) then
                           bumpup=1
                        else if (cy .ne. 0) then 
                           bumpup=2
                        else if (cz .ne. 0) then
                           bumpup=3
                        else
                           call lnkerr('m352:genacint---this cannot'//
     $                          ' be !')
                        endif
                     endif
                  endif
c
c now we know that we need to recur on the bumpup component.  Let's do it.
c we'll need the (a-1_bumpup,c), (a+1_bumpup,c) and (a,c-1_bumpup)
c remember that the cx,cy,cz are already c+1_bumpup!
c
                  oc(1)=cx
                  oc(2)=cy
                  oc(3)=cz
                  oc(bumpup)=oc(bumpup)-1
                  cil=intloc(oc(1),oc(2),oc(3))
                  oc2(1)=cx
                  oc2(2)=cy
                  oc2(3)=cz
                  oc2(bumpup)=oc2(bumpup)-2
                  c2il=intloc(oc2(1),oc2(2),oc2(3))
                  oa1(1)=ax
                  oa2(1)=ax
                  oa1(2)=ay
                  oa2(2)=ay
                  oa1(3)=az
                  oa2(3)=az
                  oa1(bumpup)=oa1(bumpup)-1
                  oa2(bumpup)=oa2(bumpup)+1
c
c so now, we know where the (la-1,lc), (la+1,lc) and (la,lc-1) start
c (tam1off,tap1off,laoff) and we can get offsets relative to those...
c
                  tam1ind=1+tam1off+intloc(oa1(1),oa1(2),oa1(3))*
     $                 mycart(lc)+cil
                  tap1ind=1+tap1off+intloc(oa2(1),oa2(2),oa2(3))*
     $                 mycart(lc)+cil
                  laind=1+laoff+ail*mycart(lc-1)+
     $                 c2il
                  if (debug) then
                     write(iout,*)"Forming (",ax,",",ay,",",az,"|",
     $                    cx,",",cy,",",cz,") in ",naoff+theket,"from"
                     if (oa1(bumpup) .ge. 0) then
                     write(iout,*)" (",oa1(1),",",oa1(2),",",oa1(3)
     $                    ,"|",oc(1),",",oc(2),",",oc(3),"), at ",
     $                    tam1ind
                     write(iout,*)(tmparry(i,tam1ind,this),i=1,nv)
                     endif
                     write(iout,*)"  (",oa2(1),",",oa2(2),",",oa2(3)
     $                    ,"|",oc(1),",",oc(2),",",oc(3),"), at ",
     $                    tap1ind
                     write(iout,*)(tmparry(i,tap1ind,this),i=1,nv)
                     if (oc(bumpup).gt.0) then
                     write(iout,*)"  (",ax,",",ay,",",az,"|",
     $                    oc2(1),",",oc2(2),",",oc2(3),") at ",laind
                     write(iout,*)(tmparry(i,laind,last),i=1,nv)
                     endif
                  endif
c
c (a,c+1_i)=-alphaA/alphaC(a+1_i,c)+N_i(a)/2alphaC (a-1_i,c)+
c     N_i(c)/2alphaC (a,c-1_i)
c
                  call vmul(tmparry(1,naoff+theket,next),
     $                 ainv(1,3),tmparry(1,tap1ind,this),nv)
                  if (debug) then
                     write(iout,*) "-alphaA/alphaC(A+1_i,c)="
                     write(iout,*) (tmparry(i,naoff+theket,next),
     $                    i=1,nv)
                  endif
                  if ( oa1(bumpup).ge. 0 ) then
c
c then a_bumpup, not easily accessed without making another stupid array,
c     is greater than or equal to 1 (coz oa1(bumpup)=a_bumpup-1)
c so form (N_i(a)/2alphaC)(a-1_i,c)
                     sclri=float(oa1(bumpup)+1.0)*half
                     if (debug) then
                        write(iout,*)"Sclri=",sclri
                     endif
                     call smul(t,ainv(1,2),sclri,nv)
                     if (debug) then
                        write(iout,*)"sclri/alphaC="
                        write(iout,*)(t(i),i=1,nv)
                     endif
                     call vwxy(tmparry(1,naoff+theket,next),
     $                    tmparry(1,naoff+theket,next),t,
     $                    tmparry(1,tam1ind,this),+1,nv)
                     if (debug) then
                        write(iout,*) "after adding in N_i(a)/2alphac",
     $                       " (a-1_i,c):"
                        write(iout,*) (tmparry(i,naoff+theket,next),
     $                    i=1,nv)
                     endif
                  endif
                  if (oc(bumpup).gt.0) then
                     sclri=float(oc(bumpup))*half
                     call smul(t,ainv(1,2),sclri,nv)
                     call vwxy(tmparry(1,naoff+theket,next),
     $                    tmparry(1,naoff+theket,next),t,
     $                    tmparry(1,laind,last),+1,nv)
                     if (debug) then
                        write(iout,*) "after adding in N_i(c)/2alphac",
     $                       " (a,c-1_i):"
                        write(iout,*) (tmparry(i,naoff+theket,next),
     $                    i=1,nv)
                     endif
                  endif
                  if (debug) then
                     write(iout,*)"And da winna iz (",ax,",",ay,",",az,
     $                    "|",cx,",",cy,",",cz,")="
                     write(iout,*)(tmparry(i,naoff+theket,next),i=1,nv)
                  endif
 50            continue 
c     get ready for next bra
               naoff=naoff+mycart(lc+1)
 40         continue 
c get ready for next value of la
            tam1off=tam1off+mycart(la-1)*mycart(lc)
            tap1off=tap1off+mycart(la+1)*mycart(lc)
            if (lc .gt. 0) then
               laoff=laoff+mycart(la)*mycart(lc-1)
            endif
 30      continue 
c 
c get ready for next value of lc
c
         temp=last
         last=this
         this=next
         next=temp
         tamin=namin
         tamax=namax
 20   continue 
c
c      Hah.  Now copy over the end results into outarry.
c     This involves figuring out what the nx,ny, and nz want us to be and
c     pulling out the appropriate piece from tmparry(....,this)
c     Off1 is offset of to first integral with the given bra relative to
c     beginning of array.
c
      theint=1
      do 60 thebra=1,mycart(imax)
         ax=nx(mini+thebra-1)
         ay=ny(mini+thebra-1)
         az=nz(mini+thebra-1)
         off1=intloc(ax,ay,az)*mycart(jmax)
         do 70 theket=1,mycart(jmax)
            cx=nx(minj+theket-1)
            cy=ny(minj+theket-1)
            cz=nz(minj+theket-1)
            off2=off1+intloc(cx,cy,cz)
            if (debug) then
               write(iout,*) "We're putting (",ax,",",ay,",",az,
     $              "|",cx,",",cy,",",cz,"), which is in ",off2+1,
     $              " into ",theint
            endif
            call vmove(outarry(1,theint),tmparry(1,off2+1,this),nv)
            theint=theint+1
 70      continue 
 60   continue 

      return
      end
