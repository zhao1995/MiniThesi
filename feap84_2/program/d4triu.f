c$Id:$
      subroutine d4triu( al, au, ad, jp, neqs, neqt )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'dlog10' to 'log10'                       17/11/2006
c       2. Add prints during 'echo' mode                    21/07/2007
c       3. Add prflg                                        16/12/2008
c       4. Add 'setups.h' and check rank for prints         24/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Triangular decomposition of a sym./unsym. matrix
c               stored in profile form.  Programmed for super-scalar
c               processors.

c      Inputs:
c         al(*) - Unreduced rows    of lower triangular part.
c         au(*) - Unreduced columns of upper triangular part.
c         ad(*) - Unreduced diagonals of array.
c         jp(*) - Pointers to bottom of row/columns of 'al','au'.
c         neqs  - Number of symmetric equations (1 or more).
c         neqt  - Number of total equations to be solved.
c                 (N.B. Last "neqt - neqs" equations are unsymmetric)

c      Outputs:
c         al(*) - Reduced rows    of lower triangular part.
c         au(*) - Reduced columns of upper triangular part.
c         ad(*) - Reduced diagonals of array.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'counts.h'
      include  'iofile.h'
      include  'fdata.h'
      include  'print.h'
      include  'setups.h'

      logical   prflg
      integer   i, ie, ic, j, j1, j2, je, jf, jq, js, jt, jp(*), jh(2)
      integer   n0, n1, n2, neq, neqs, neqt, nprti, nerrs, ifig
      real*4    etime, tt, tary(2)
      real*8    dd,at2,tol, dot, dsred, dured, dfig,dimx,dimn
      real*8    al(*), au(*), ad(*)

      save

      data      tol /0.5d-07/
      data      nerrs /0/

c     Reduce 1 equation if neq odd

      n0  = 0
      n1  = 0
      n2  = 0
      dfig = 1.0d0
      dimx = 0.0d0
      dimn = 0.0d0
      do j = 1,neqt
        dimn = max(dimn,abs(ad(j)))
      end do ! j
      prflg = ior.lt.0 .and. pfr .and. rank.eq.0

      neq   = max(neqs,1)
      nprti = max(neqt/10,1001)

      if(ad(1).ne.0.0d0) then
        dimx  = max(dimx,abs(ad(1)))
        dimn  = min(dimn,abs(ad(1)))
        ad(1) = 1.d0/ad(1)
      endif
      js = 2
      if(mod(neq,2).eq.0) then
        if(jp(2).eq.1) then
          at2   = au(1)
          au(1) = au(1)*ad(1)
          dd    = ad(2)
          ad(2) = ad(2) - at2*au(1)

c         Count errors

          if(abs(ad(2)).lt.tol*abs(dd))    n0 = n0 + 1
          if(dd.lt.0.d0.and.ad(2).gt.0.d0) n1 = n1 + 1
          if(dd.gt.0.d0.and.ad(2).lt.0.d0) n1 = n1 + 1

        endif
        if(ad(2).eq.0.d0) then
          if(jp(2) .gt.0) n2 = n2 + 1
        else
          dimx  = max(dimx,abs(ad(2)))
          dimn  = min(dimn,abs(ad(2)))
          dfig  = max(dfig,abs(dd/ad(2)))
          ad(2) = 1.d0/ad(2)
        endif
        js = 3
      endif

c     Set a threshold to initiate reduction

      jt = 5

c     Reduce remaining symmetric equations in pairs

      if(prflg .and. neq.ge.1000) then
        write(*,3000)
      endif

      do j = js, neq, 2

        if(prflg .and. (mod(j,nprti).eq.0.or.mod(j+1,nprti).eq.0)) then
          tt = etime(tary)
          write(*,3001) j,tary
        endif

        jh(1) = jp(j  ) - jp(j-1)
        jh(2) = jp(j+1) - jp(j  ) - 1
        ic    = min(jh(1),jh(2))

c       Reduce short column cases

        if(ic.lt.jt) then

          do i = 0,1
            j1      = jp(j-1+i) + 1
            call dlfac( au(j1), au, jp(j-jh(i+1)), jh(i+1)+i)
            dd      = ad(j+i)
            ad(j+i) = ad(j+i) - dsred(au(j1),ad(j-jh(i+1)),jh(i+1)+i)

c           Count errors

            if(abs(ad(j+i)) .lt.    tol*abs(dd))    n0 = n0 + 1
            if(dd.lt.0.d0   .and.  ad(j+i).gt.0.d0) n1 = n1 + 1
            if(dd.gt.0.d0   .and.  ad(j+i).lt.0.d0) n1 = n1 + 1
            if(ad(j+i).eq.0.d0 .and. jh(i+1).gt.-i) n2 = n2 + 1

            if(ad(j+i).ne.0.0d0) then
              dimx  = max(dimx,abs(ad(j+i)))
              dimn  = min(dimn,abs(ad(j+i)))
              dfig  = max(dfig,abs(dd/ad(j+i)))
              ad(j+i) = 1.d0/ad(j+i)
            endif
          end do ! i

        else

          ie    = j - ic
          ie    = min(j, ie + mod( j+ie,2 ))

c         Reduce column j above ie

          call dlfac( au(jp(j-1)+1), au, jp(j-jh(1)), ie-j+jh(1) )

c         Reduce column j+1 above ie

          call dlfac( au(jp(j  )+1), au, jp(j-jh(2)), ie-j+jh(2) )

c         Reduce rows associated with column pairs

          call d4blks( au , jp , jh, ie, j)

c         Reduce diagonal of column j

          dd      = ad(j)
          ad(j  ) = ad(j  ) - dsred(au(jp(j-1)+1),ad(j-jh(1)),jh(1))

c         Count errors

          if(abs(ad(j)) .lt. tol*abs(dd))    n0 = n0 + 1
          if(dd.lt.0.d0 .and. ad(j).gt.0.d0) n1 = n1 + 1
          if(dd.gt.0.d0 .and. ad(j).lt.0.d0) n1 = n1 + 1
          if(ad(j).eq.0.d0 .and. jh(1).gt.0) n2 = n2 + 1

          if(ad(j  ).ne.0.0d0) then
            dimx    = max(dimx,abs(ad(j)))
            dimn    = min(dimn,abs(ad(j)))
            dfig    = max(dfig,abs(dd/ad(j)))
            ad(j  ) = 1.d0/ad(j  )
          endif

c         Reduce last term in column j+1

          au(jp(j+1)) = au(jp(j+1))
     &                - dot(au(jp(j)+1-ic),au(jp(j+1)-ic),ic)

c         Reduce diagonal of column j+1

          dd      = ad(j+1)
          ad(j+1) = ad(j+1) - dsred(au(jp(j)+1),ad(j-jh(2)),jh(2)+1)

c         count errors

          if(abs(ad(j+1)) .lt. tol*abs(dd))    n0 = n0 + 1
          if(dd.lt.0.d0 .and. ad(j+1).gt.0.d0) n1 = n1 + 1
          if(dd.gt.0.d0 .and. ad(j+1).lt.0.d0) n1 = n1 + 1
          if(ad(j+1).eq.0.d0 .and. jh(2).ge.0) n2 = n2 + 1

          if(ad(j+1).ne.0.0d0) then
            dimx    = max(dimx,abs(ad(j+1)))
            dimn    = min(dimn,abs(ad(j+1)))
            dfig    = max(dfig,abs(dd/ad(j+1)))
            ad(j+1) = 1.d0/ad(j+1)
          endif
        endif

      end do ! j

c     Reduce unsymmetric equations in pairs

      je = jp(neq)
      jf = 1 - je
      do j = max(js,neq+1), neqt-mod(neqt-max(js-1,neq),2), 2

        if(prflg .and. (mod(j,nprti).eq.0.or.mod(j+1,nprti).eq.0)) then
          tt = etime(tary)
          write(*,3001) j,tary
        endif

        jh(1) = jp(j  ) - jp(j-1)
        jh(2) = jp(j+1) - jp(j  ) - 1
        ic    = min(jh(1),jh(2))

c       Reduce short column cases

        if(ic.lt.jt) then

          jq = j - neq - 1
          do i = 0,1
            j1      = jp(j-1+i) + 1
            call dufac(au(j1),al,au, jp(j-jh(i+1)), jh(i+1)+i, je,jq+i)
            call dlfac(al(j1-je),au, jp(j-jh(i+1)), jh(i+1)+i )
            dd      = ad(j+i)
            ad(j+i) = ad(j+i) - dured(al(j1-je),au(j1),ad(j-jh(i+1)),
     &                              jh(i+1)+i)

c           Count errors

            if(abs(ad(j+i)) .lt.    tol*abs(dd))    n0 = n0 + 1
            if(dd.lt.0.d0   .and.  ad(j+i).gt.0.d0) n1 = n1 + 1
            if(dd.gt.0.d0   .and.  ad(j+i).lt.0.d0) n1 = n1 + 1
            if(ad(j+i).eq.0.d0 .and. jh(i+1).gt.-i) n2 = n2 + 1

            if(ad(j+i).ne.0.0d0) then
              dimx  = max(dimx,abs(ad(j+i)))
              dimn  = min(dimn,abs(ad(j+i)))
              dfig  = max(dfig,abs(dd/ad(j+i)))
              ad(j+i) = 1.d0/ad(j+i)
            endif
          end do ! j

        else

          ie = j - ic
          ie = min(j, ie + mod( j+ie,2 ))
          jq = max(0, ie - neq - 1)

c         Reduce row/column j above ie

          call dlfac(al(jp(j-1)+jf),au,   jp(j-jh(1)),ie-j+jh(1) )
          call dufac(au(jp(j-1)+1 ),al,au,jp(j-jh(1)),ie-j+jh(1),je,jq)

c         Reduce row/column j+1 above ie

          call dlfac(al(jp(j)+jf),au,   jp(j-jh(2)),ie-j+jh(2) )
          call dufac(au(jp(j)+1 ),al,au,jp(j-jh(2)),ie-j+jh(2),je,jq)

c         Reduce rows associated with column pairs

          call d4blku( al , au , jp , jh, ie, j, je, neq)

c         Reduce diagonal of column j

          dd      = ad(j)
          ad(j  ) = ad(j  ) - dured(al(jp(j-1)+jf),au(jp(j-1)+1),
     &                              ad(j-jh(1)),jh(1))

c         Count errors

          if(abs(ad(j)) .lt. tol*abs(dd))    n0 = n0 + 1
          if(dd.lt.0.d0 .and. ad(j).gt.0.d0) n1 = n1 + 1
          if(dd.gt.0.d0 .and. ad(j).lt.0.d0) n1 = n1 + 1
          if(ad(j).eq.0.d0 .and. jh(1).gt.0) n2 = n2 + 1

          if(ad(j  ).ne.0.0d0) then
            dimx    = max(dimx,abs(ad(j)))
            dimn    = min(dimn,abs(ad(j)))
            dfig    = max(dfig,abs(dd/ad(j)))
            ad(j  ) = 1.d0/ad(j  )
          endif

c         Reduce last term in row/column j+1

          al(jp(j+1)-je) = al(jp(j+1)-je)
     &                - dot(au(jp(j)+1-ic),al(jp(j+1)-ic-je),ic)
          au(jp(j+1)   ) = au(jp(j+1)   )
     &                - dot(al(jp(j)+1-ic-je),au(jp(j+1)-ic),ic)

c         Reduce diagonal of column j+1

          dd      = ad(j+1)
          ad(j+1) = ad(j+1) - dured(al(jp(j )+jf),au(jp(j  )+1),
     &                              ad(j-jh(2)),jh(2)+1)

c         Count errors

          if(abs(ad(j+1)) .lt. tol*abs(dd))    n0 = n0 + 1
          if(dd.lt.0.d0 .and. ad(j+1).gt.0.d0) n1 = n1 + 1
          if(dd.gt.0.d0 .and. ad(j+1).lt.0.d0) n1 = n1 + 1
          if(ad(j+1).eq.0.d0 .and. jh(2).ge.0) n2 = n2 + 1

          if(ad(j+1).ne.0.0d0) then
            dimx    = max(dimx,abs(ad(j+1)))
            dimn    = min(dimn,abs(ad(j+1)))
            dfig    = max(dfig,abs(dd/ad(j+1)))
            ad(j+1) = 1.d0/ad(j+1)
          endif
        endif

      end do ! j

c     Reduce last equation if (neqt - neqs) odd

      if( mod(neqt-max(neq,js-1),2).ne.0 ) then
        j  = neqt
        j1 = jp(j-1) + 1
        j2 = jp(j) - jp(j-1)
        jq = j - neq - 1
        call dufac( au(j1)   , al,au,jp(j-j2), j2, je, jq )
        call dlfac( al(j1-je), au,   jp(j-j2), j2 )
        dd    = ad(j)
        ad(j) = ad(j) - dured(al(j1-je),au(j1),ad(j-j2),j2)

c       Count errors

        if(abs(ad(j)).lt.tol*abs(dd))    n0 = n0 + 1
        if(dd.lt.0.d0.and.ad(j).gt.0.d0) n1 = n1 + 1
        if(dd.gt.0.d0.and.ad(j).lt.0.d0) n1 = n1 + 1
        if(ad(j).eq.0.d0 .and. j2.gt.0)  n2 = n2 + 1

        if(ad(j).ne.0.0d0) then
          dimx    = max(dimx,abs(ad(j)))
          dimn    = min(dimn,abs(ad(j)))
          dfig    = max(dfig,abs(dd/ad(j)))
          ad(j) = 1.d0/ad(j)
        endif
      endif

c     Output error summary and conditioning information

      if(n0+n1+n2.gt.0) then
        nerrs = nerrs + 1
      endif

      if(nerrs.lt.40) then
        if(n0.gt.0) write(iow,2000) n0,nstep,niter
        if(n1.gt.0) write(iow,2001) n1,nstep,niter
        if(n2.gt.0) write(iow,2002) n2,nstep,niter
        if(ior.lt.0 .or. echo) then
          if(n0.gt.0) write(*,2000) n0,nstep,niter
          if(n1.gt.0) write(*,2001) n1,nstep,niter
          if(n2.gt.0) write(*,2002) n2,nstep,niter
        endif
      endif

c     Output conditioning check

      if(pfr .and. rank.eq.0) then
        dd = 0.0d0
        if(dimn.ne.0.0d0) dd = dimx/dimn
        ifig = nint(log10(dfig) + 0.6d0)
        write(iow,2003) dimx,dimn,dd,ifig
        if(ior.lt.0) then
          write(*,2003) dimx,dimn,dd,ifig
        endif
      endif

c     Formats

2000  format(' *D4TRI WARNING* Lost at least 7 digits in',
     & ' reducing diagonals of',i10,' equations.'/
     &   16x,' Step:',i7,' Iteration:',i6)

2001  format(' *D4TRI WARNING* Sign of diagonal changed when',
     & ' reducing',i10,' equations.'/16x,' Step:',i7,' Iteration:',i6)

2002  format(' *D4TRI WARNING* Reduced diagonal is zero in',
     & i10,' equations.'/16x,' Step:',i7,' Iteration:',i6)

2003  format(' Condition check: D-max',e11.4,'; D-min',e11.4,
     & '; Ratio',e11.4/' Maximum no. diagonal digits lost:',i3)

3000  format('   Solution Status')
3001  format('   --> Equation = ',i10,' Time: CPU = ',f12.2,
     &       ' System = ',f12.2)

      end
