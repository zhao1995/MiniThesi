c$Id:$
      subroutine cdatri(alr,ali,aur,aui,adr,adi,jp,neqs,neqt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'dlog10' to 'log10'                       17/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Triangular decomposition of a matrix stored in profile
c               form for solution of algebraic equations stored in
c               profile form

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
c         alr(*) - Real      part of lower terms in unfactored array
c         ali(*) - Imaginary part of lower terms in unfactored array
c         aur(*) - Real      part of upper terms in unfactored array
c         aui(*) - Imaginary part of upper terms in unfactored array
c         adr(*) - Real      part of diagonal terms in unfactored array
c         adi(*) - Imaginary part of diagonal terms in unfactored array
c         jp(*)  - Pointer to end of columns/rows in A array
c         neqs   - Number of symmetric equations
c         neqt   - Number of total     equations

c      Outputs:
c         alr(*) - Real      part of lower terms in factored array
c         ali(*) - Imaginary part of lower terms in factored array
c         aur(*) - Real      part of upper terms in factored array
c         aui(*) - Imaginary part of upper terms in factored array
c         adr(*) - Real      part of diagonal terms in factored array
c         adi(*) - Imaginary part of diagonal terms in factored array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'fdata.h'
      include  'iofile.h'

      integer   neqs, neqt, neq , n0, n1, jp(*)
      integer   i, id, ie, ih, is, ifig, j, jd, je, jh, jr, jrh
      integer   ijh(2), n
      real*4    etime,tary(2),tt
      real*8    dfig, dimx,dimn, dd, del0, delj
      real*8    srr, sri, sir, sii
      real*8    alr(*),aur(*),adr(*)
      real*8    ali(*),aui(*),adi(*)

      save

c     Set error counters and initial values for checks

      n0   = 0
      n1   = 0
      dfig = 1.d0
      dimx = 0.d0
      dimn = 0.d0
      do j = 1,neqt
        dimn = max(dimn,abs(adr(j)),abs(adi(j)))
      end do ! j

c     Loop through columns to perform triangular decomposition

      neq = max(1,neqs)

c     Symetric part

      jd = 1
      do j = 1,neq
        jr = jd + 1
        jd = jp(j)
        jh = jd - jr
        if(jh.gt.0) then
          is = j - jh
          ie = j - 1

c         Construct upper triangular matrix

          do i = is,min(ie,neq)
            jr = jr + 1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              ijh(1) = jr - ih
              ijh(2) = id - ih + 1
              srr    = 0.0d0
              sri    = 0.0d0
              sir    = 0.0d0
              sii    = 0.0d0
              do n = 0,ih-1
                srr = srr + aur(ijh(1)+n)*aur(ijh(2)+n)
                sri = sri + aur(ijh(1)+n)*aui(ijh(2)+n)
                sir = sir + aui(ijh(1)+n)*aur(ijh(2)+n)
                sii = sii + aui(ijh(1)+n)*aui(ijh(2)+n)
              end do ! n
              aur(jr) = aur(jr) - srr + sii
              aui(jr) = aui(jr) - sri - sir
            endif
          end do ! i
        endif

c       Reduce diagonal

        dd = adr(j)
        del0 = adr(j)*adr(j) + adi(j)*adi(j)
        if(jh.ge.0) then
          jrh = j - jh - 1
          call cdsred(aur(jd-jh),aui(jd-jh),
     &                adr(jrh),adi(jrh),jh+1, adr(j),adi(j))
        endif

c       Count errors

        delj  = adr(j)*adr(j) + adi(j)*adi(j)

        if(delj .lt. 0.5d-7*del0) n0 = n0 + 1
        if(delj .eq. 0.0d0) then
          if(jh.ge.0) n1 = n1 + 1

c       Store reciprocal of diagonal

        else
          dimx  = max(dimx,abs(adr(j)))
          dimn  = min(dimn,abs(adr(j)))
          dfig  = max(dfig,del0/delj)
          adr(j) =   adr(j)/delj
          adi(j) = - adi(j)/delj
        endif

c       Output interactive equation processing

        if(ior.lt.0 .and. mod(j,1000).eq.0 ) then
          tt = etime(tary)
          write(*,3000) j,tary
        endif
      end do ! j

c     Unsymmetric part

      je  = jp(neq)
      jd  = je
      do j = neq+1,neqt
        jr = jd + 1
        jd = jp(j)
        jh = jd - jr
        if(jh.gt.0) then
          is = j - jh
          ie = j - 1

c         Construct upper and lower triangular matrix

          do i = is,min(ie,neq)
            jr = jr + 1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              ijh(1) = jr - ih
              ijh(2) = id - ih + 1
              srr    = 0.0d0
              sri    = 0.0d0
              sir    = 0.0d0
              sii    = 0.0d0
              do n = 0,ih-1
                srr = srr + aur(ijh(1)+n)*aur(ijh(2)+n)
                sri = sri + aur(ijh(1)+n)*aui(ijh(2)+n)
                sir = sir + aui(ijh(1)+n)*aur(ijh(2)+n)
                sii = sii + aui(ijh(1)+n)*aui(ijh(2)+n)
              end do ! n
              aur(jr) = aur(jr) - srr + sii
              aui(jr) = aui(jr) - sri - sir

              ijh(1) = jr - ih - je
              srr    = 0.0d0
              sri    = 0.0d0
              sir    = 0.0d0
              sii    = 0.0d0
              do n = 0,ih-1
                srr = srr + alr(ijh(1)+n)*aur(ijh(2)+n)
                sri = sri + alr(ijh(1)+n)*aui(ijh(2)+n)
                sir = sir + ali(ijh(1)+n)*aur(ijh(2)+n)
                sii = sii + ali(ijh(1)+n)*aui(ijh(2)+n)
              end do ! n
              alr(jr-je) = alr(jr-je) - srr + sii
              ali(jr-je) = ali(jr-je) - sri - sir
            endif
          end do ! i

          do i = max(neq+1,is),ie
            jr = jr + 1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              ijh(1) = jr - ih
              ijh(2) = id - ih + 1
              srr    = 0.0d0
              sri    = 0.0d0
              sir    = 0.0d0
              sii    = 0.0d0
              do n = 0,ih-1
                srr = srr + aur(ijh(1)+n)*alr(ijh(2)+n)
                sri = sri + aur(ijh(1)+n)*ali(ijh(2)+n)
                sir = sir + aui(ijh(1)+n)*alr(ijh(2)+n)
                sii = sii + aui(ijh(1)+n)*ali(ijh(2)+n)
              end do ! n
              aur(jr) = aur(jr) - srr + sii
              aui(jr) = aui(jr) - sri - sir

              ijh(1) = jr - ih - je
              srr    = 0.0d0
              sri    = 0.0d0
              sir    = 0.0d0
              sii    = 0.0d0
              do n = 0,ih-1
                srr = srr + alr(ijh(1)+n)*aur(ijh(2)+n)
                sri = sri + alr(ijh(1)+n)*aui(ijh(2)+n)
                sir = sir + ali(ijh(1)+n)*aur(ijh(2)+n)
                sii = sii + ali(ijh(1)+n)*aui(ijh(2)+n)
              end do ! n
              alr(jr-je) = alr(jr-je) - srr + sii
              ali(jr-je) = ali(jr-je) - sri - sir
            endif
          end do ! i
        endif

c       Reduce diagonal

        dd = adr(j)
        del0 = adr(j)*adr(j) + adi(j)*adi(j)
        if(jh.ge.0) then
          jrh = j - jh - 1
          call cdured(alr(jd-jh),ali(jd-jh),aur(jd-jh),aui(jd-jh),
     &                adr(jrh),adi(jrh),jh+1, adr(j),adi(j))
        endif


c       Count errors

        delj  = adr(j)*adr(j) + adi(j)*adi(j)

        if(delj .lt. 0.5d-7*del0) n0 = n0 + 1
        if(delj .eq. 0.0d0) then

          if(jh.ge.0) n1 = n1 + 1

c       Store reciprocal of diagonal

        else
          dimx  = max(dimx,abs(adr(j)))
          dimn  = min(dimn,abs(adr(j)))
          dfig  = max(dfig,del0/delj)
          adr(j) =   adr(j)/delj
          adi(j) = - adi(j)/delj
        endif

c       Output interactive equation processing

        if(ior.lt.0 .and. mod(j,1000).eq.0 ) then
          write(*,3000) j,tary
        endif
      end do ! j

c     Output error summary and check information

      dd = 0.d0
      if(dimn.ne.0.d0) dd = dimx/dimn
      ifig = nint(log10(dfig) + 0.6d0)
      if(n0.gt.0) write(iow,2000) n0
      if(n1.gt.0) write(iow,2001) n1
      if(pfr) write(iow,2002) dimx,dimn,dd,ifig
      if(ior.lt.0) then
        if(n0.gt.0) write(*,2000) n0
        if(n1.gt.0) write(*,2001) n1
        if(pfr) write(*,2002) dimx,dimn,dd,ifig
      endif

c     Formats

2000  format(' *WARNING* loss of at least 7 digits in',
     1 ' reducing diagonals of',i5,' equations.')

2001  format(' *WARNING* reduced diagonal is zero in',
     1 i5,' equations.')

2002  format(' Condition check: D-max',e11.4,'; D-min',e11.4,
     1 '; Ratio',e11.4/' Maximum no. diagonal digits lost:',i3)

3000  format(' ---> Equation = ',i5:,' , Time: CPU = ',f12.2,
     1       ' , System = ',f12.2)

      end
