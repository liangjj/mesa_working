        program mcomelp
c
c       ===================================================
c       purpose: this program computes complete elliptic 
c                integrals k(k) and e(k) using subroutine 
c                comelp
c       input  : k  --- modulus k ( 0 ó k ó 1 )
c       output : ck --- k(k)
c                ce --- e(k)
c       example:
c                  k         k(k)          e(k)
c                ---------------------------------
c                 .00      1.570796      1.570796
c                 .25      1.596242      1.545957
c                 .50      1.685750      1.467462
c                 .75      1.910990      1.318472
c                1.00       ì            1.000000
c       ===================================================
c
        double precision hk,ck,ce
        write(*,*)'please enter the modulus k '
        read(*,*) hk
        write(*,*)'    k         k(k)          e(k)'
        write(*,*)'  ---------------------------------'
        call comelp(hk,ck,ce)
        if (hk.ne.1.0) write(*,10) hk,ck,ce
        if (hk.eq.1.0) write(*,20) hk,ce
10      format(2x,f5.2,2f14.6)
20      format(2x,f5.2,7x,'ì',6x,f14.6)
        end


        subroutine comelp(hk,ck,ce)
c
c       ==================================================
c       purpose: compute complete elliptic integrals k(k)
c                and e(k)
c       input  : k  --- modulus k ( 0 ó k ó 1 )
c       output : ck --- k(k)
c                ce --- e(k)
c       ==================================================
c
        implicit double precision (a-h,o-z)
        pk=1.0d0-hk*hk
        if (hk.eq.1.0) then
           ck=1.0d+300
           ce=1.0d0
        else
           ak=(((.01451196212d0*pk+.03742563713d0)*pk
     &        +.03590092383d0)*pk+.09666344259d0)*pk+
     &        1.38629436112d0
           bk=(((.00441787012d0*pk+.03328355346d0)*pk+
     &        .06880248576d0)*pk+.12498593597d0)*pk+.5d0
           ck=ak-bk*dlog(pk)
           ae=(((.01736506451d0*pk+.04757383546d0)*pk+
     &        .0626060122d0)*pk+.44325141463d0)*pk+1.0d0
           be=(((.00526449639d0*pk+.04069697526d0)*pk+
     &        .09200180037d0)*pk+.2499836831d0)*pk
           ce=ae-be*dlog(pk)
        endif
        return
        end
