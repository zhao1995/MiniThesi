c$Id:$
      subroutine datri(al,au,ad,jp, neqs, neqt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'dlog10' to 'log10'                       17/11/2006
c       2. Allocate both SPTAN and IPTAN for sparse         01/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Triangular decomposition of matrix stored in profile or
c               sparse form.
c               N.B. Sparse form for symmetric only.

c         Equations '   1  ' to 'neqs' are symmetric.
c         Equations 'neqs+1' to 'neqt' are unsymmetric.

c         Use:
c          a.) All equations are unsymmetric       : neqs = 1 (or 0)
c              N.B.  The top 1 x 1 submatrix is always symmetric.
c                    Both 'al' and 'au' must be provided.

c          b.) All equations are symmetric         : neqs = neqt
c              N.B.  In this case array 'al' is not used.

c          c.) First 'neqs' equations are symmetric: 1 < neqs < neqt
c              N.B.  Storage of 'al' for unsymmetric equations only.

c      Inputs:
c         al(*)  - Unfactored lower triangular terms
c         au(*)  - Unfactored upper triangular terms
c         ad(*)  - Unfactored diagonal terms
c         jp(*)  - Pointer to row/column ends in 'al' and 'au'.
c         neqs   - Number of symmetric equations
c         neqt   - Number of equations in A.

c      Outputs:
c         al(*)  - Factored lower triangular terms
c         au(*)  - Factored upper triangular terms
c         ad(*)  - Factored diagonal terms
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'compas.h'
      include  'complx.h'
      include  'counts.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'print.h'
      include  'setups.h'
      include  'ssolve.h'
      include  'comblk.h'

      logical   setval, palloc
      integer   espo, neqs, neqt, neq , n0, n1, n2, n3, jp(*)
      integer   i, id, ie, ih, is, idh, ifig, j, jd, je, jh, jr, jrh
      real*4    etime, tt, tary(2)
      real*8    tol, dd, daval, dfig, dimx, dimn, dot, dsred, dured
      real*8    al(*),au(*),ad(*)

      save

c     N.B.  tol should be set to approximate half-word precision.

      data      tol/0.5d-07/

c     Factor for complex equations

      if(cplxfl) then
        is = jp(neqt) + 1
        id = neqt     + 1
        call cdatri(al,al(is),au,au(is),ad,ad(id),jp,neqs,neqt)
        return
      endif

c     Use Sparse solver with minimum degree ordering (symmetric only)

      if(ittyp.eq.-2) then

        esp    = 3*neqt + 4*mr(np(93)+2*neqt) + 1
        setval = palloc(111, 'TEMP1' ,2*neqt,  1) ! temp variables
        if(np(67).ne.0) then
          setval = palloc( 67, 'SPTAN' ,   0  ,  2) ! Delete present one
          setval = palloc(307, 'IPTAN' ,   0  ,  1) ! Delete present one
        endif
        setval = palloc( 67, 'SPTAN' , esp  ,  2) ! Add new one at end
        setval = palloc(307, 'IPTAN' , esp  ,  1) ! Add new one at end
        espo = esp

c       With minimum degree

        if(domd .and. ddom) then
          call odrv(neqt,mr(np(93)+neqt),mr(np(94)),ad,mr(np(47)),
     &              mr(np(48)),esp,mr(np(307)),2,je)
          if (je.ne.0) call oderr(je)
          ddom = .false.
        endif

c       Symbolic triangular decomposition

        call sdrv(neqt,mr(np(47)),mr(np(48)),mr(np(93)+neqt),mr(np(94)),
     &            ad,ad,ad,mr(np(307)),hr(np(67)),esp,4,je)

c       Minimal dimension for working arrays

        if(esp.gt.espo) then
          setval = palloc( 67, 'SPTAN' , esp, 2)
          setval = palloc(307, 'IPTAN' , esp, 1)
        endif
        setval = palloc(111, 'TEMP1' ,   0, 1)

c       Numerical triangular decomposition

        call sdrv(neqt,mr(np(47)),mr(np(48)),mr(np(93)+neqt),mr(np(94)),
     &            ad,ad,ad,mr(np(307)),hr(np(67)),esp,6,je)
        if (je.ne.0) call sderr(je,neqt)
        return

c     Use two column elimination

      elseif(soltyp.eq.2) then
        call d4triu(al,au,ad,jp, neqs, neqt)
        return
      endif

c     Set error counters and initial values for conditioning check

      n0   = 0
      n1   = 0
      n2   = 0
      n3   = 0
      dfig = 1.0d0
      dimx = 0.0d0
      dimn = 0.0d0
      do j = 1,neqt
        dimn = max(dimn,abs(ad(j)))
      end do ! j

c     Loop through columns to perform triangular decomposition

      neq = max(1,neqs)

c     Reduce symmetric equations

      jd  = 1
      do j = 1,neq
        jr = jd + 1
        jd = jp(j)
        jh = jd - jr
        if(jh.gt.0) then
          is = j - jh
          ie = j - 1

c         if diagonal is zero compute a norm for singularity test

          if(ad(j).eq.0.0d0) call datest(au(jr),jh,daval)

          do i = is,min(ie,neq)
            jr = jr + 1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              jrh = jr - ih
              idh = id - ih + 1
              au(jr) = au(jr) - dot( au(jrh), au(idh), ih )
            endif
          end do ! i
        endif

c       Reduce diagonal

        dd = ad(j)
        if(jh.ge.0) then
          jr    = jd - jh
          jrh   = j - jh - 1
          ad(j) = ad(j) - dsred(au(jr),ad(jrh),jh+1)

c         Count errors

          if(abs(ad(j)).lt.tol*abs(dd))      n0 = n0 + 1
          if(dd.lt.0.0d0.and.ad(j).gt.0.0d0) n1 = n1 + 1
          if(dd.gt.0.0d0.and.ad(j).lt.0.0d0) n1 = n1 + 1
          if(ad(j) .eq.  0.0d0)              n2 = n2 + 1

c         Complete rank test for a zero diagonal case

          if(dd.eq.0.0d0.and.jh.gt.0) then
            if(abs(ad(j)).lt.tol*daval)      n3 = n3 + 1
          endif
        endif

c       Store reciprocal of diagonal

        if(ad(j).ne.0.0d0) then
          dimx  = max(dimx,abs(ad(j)))
          dimn  = min(dimn,abs(ad(j)))
          dfig  = max(dfig,abs(dd/ad(j)))
          ad(j) = 1.0d0/ad(j)
        endif

c       Output interactive equation processing

        if(ior.lt.0 .and. mod(j,1000).eq.0 .and. prnt) then
          tt = etime(tary)
          write(*,3000) j,tary
        endif
      end do ! j

c     Reduce unsymmetric equations

      je  = jp(neq)
      jd  = je
      do j = neq+1,neqt
        jr = jd + 1
        jd = jp(j)
        jh = jd - jr
        if(jh.gt.0) then
          is = j - jh
          ie = j - 1

c         If diagonal is zero compute a norm for singularity test

          if(ad(j).eq.0.0d0) call datest(au(jr),jh,daval)

          do i = is,min(ie,neq)
            jr = jr + 1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              jrh = jr - ih
              idh = id - ih + 1
              au(jr)    = au(jr   ) - dot( au(jrh   ), au(idh), ih )
              al(jr-je) = al(jr-je) - dot( al(jrh-je), au(idh), ih )
            endif
          end do ! i

          do i = max(neq+1,is),ie
            jr = jr + 1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              jrh = jr - ih
              idh = id - ih + 1
              au(jr)    = au(jr)    - dot( au(jrh), al(idh-je), ih )
              al(jr-je) = al(jr-je) - dot( al(jrh-je), au(idh), ih )
            endif
          end do ! i
        endif

c       Reduce diagonal

        dd = ad(j)
        if(jh.ge.0) then
          jr    = jd - jh
          jrh   = j - jh - 1
          ad(j) = ad(j) - dured(al(jr-je),au(jr),ad(jrh),jh+1)

c         Count errors

          if(abs(ad(j)).lt.tol*abs(dd))      n0 = n0 + 1
          if(dd.lt.0.0d0.and.ad(j).gt.0.0d0) n1 = n1 + 1
          if(dd.gt.0.0d0.and.ad(j).lt.0.0d0) n1 = n1 + 1
          if(ad(j) .eq.  0.0d0)              n2 = n2 + 1

c         Complete rank test for a zero diagonal case

          if(dd.eq.0.0d0.and.jh.gt.0) then
            if(abs(ad(j)).lt.tol*daval)      n3 = n3 + 1
          endif
        endif

c       Store reciprocal of diagonal

        if(ad(j).ne.0.0d0) then
          dimx  = max(dimx,abs(ad(j)))
          dimn  = min(dimn,abs(ad(j)))
          dfig  = max(dfig,abs(dd/ad(j)))
          ad(j) = 1.0d0/ad(j)
        endif

c       Output interactive equation processing

        if(ior.lt.0 .and. mod(j,1000).eq.0 .and. prnt) then
          tt = etime(tary)
          write(*,3000) j,tary
        endif
      end do ! j

c     Print error summary and conditioning information

      dd = 0.0d0
      if(dimn.ne.0.0d0) dd = dimx/dimn
      if(dfig.gt.0.0d0) then
        ifig = nint(log10(dfig) + 0.6d0)
      else
        ifig = -1
      endif
      if(prnt) then
        if(n0.gt.0) write(iow,2000) n0,nstep,niter
        if(n1.gt.0) write(iow,2001) n1,nstep,niter
        if(n2.gt.0) write(iow,2002) n2,nstep,niter
        if(n3.gt.0) write(iow,2003) n3,nstep,niter
        if(pfr) write(iow,2004) dimx,dimn,dd,ifig
        if(ior.lt.0) then
          if(n0.gt.0) write(*,2000) n0,nstep,niter
          if(n1.gt.0) write(*,2001) n1,nstep,niter
          if(n2.gt.0) write(*,2002) n2,nstep,niter
          if(n3.gt.0) write(*,2003) n3,nstep,niter
          if(pfr) write(*,2004) dimx,dimn,dd,ifig
        endif
      endif
          if(n1.gt.0) write(*,2001) n1,nstep,niter

c     Formats

2000  format(' *WARNING* Lost at least 7 digits in',
     &       ' reducing diagonals of',i10,' equations.'/
     &       ' Step:',i7,' Iteration:',i6)

2001  format(' *WARNING* sign of diagonal changed when reducing',
     &    i10,' equations.'/' Step:',i7,' Iteration:',i6)

2002  format(' *WARNING* reduced diagonal is zero in',
     &    i10,' equations.'/' Step:',i7,' Iteration:',i6)

2003  format(' *WARNING* rank failure for zero unreduced diagonal in',
     &    i10,' equations.'/' Step:',i7,' Iteration:',i6)

2004  format(' Condition check: D-max',e11.4,'; D-min',e11.4,
     & '; Ratio',e11.4/' Maximum no. diagonal digits lost:',i3)

3000  format(' ---> Equation = ',i10,' , Time: CPU = ',f12.2,
     &       ' , System = ',f12.2)

      end
