c$Id:$
      subroutine cgapt2(xi0,gn0, xm,cxs,nn, mel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Keep point within slave facet in 'slavt2' fn.    18/07/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor             6 March 2003            1.0

c      Purpose: Compute intersection segments

c      Inputs :

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      logical    noconv
      integer    mel, m,nc,i
      real*8     xi0,gn0, dxi,dgn, denom, tol
      real*8     shpm(2,4),xm(2,*),cxs(2),cxm(2),dxm(2),nn(2)
      real*8     ff(2)

      data       tol / 1.0d-12 /

      nc     = 0
      noconv = .true.
      do while (noconv)
        call shp1dn(xi0-1.d0,shpm,mel)
        do i = 1,2
          cxm(i) = 0.0d0
          dxm(i) = 0.0d0
          do m = 1,mel
            cxm(i) = cxm(i) + shpm(2,m)*xm(i,m)
            dxm(i) = dxm(i) + shpm(1,m)*xm(i,m)
          end do ! m
          ff(i) = cxs(i) + nn(i)*gn0 - cxm(i)
        end do ! i
        denom = nn(1)*dxm(2) - nn(2)*dxm(1)
        dxi   = ( nn(1)*ff(2) -  nn(2)*ff(1))/denom
        dgn   = (dxm(1)*ff(2) - dxm(2)*ff(1))/denom
        xi0   = xi0 + dxi
        gn0   = gn0 + dgn
        if(abs(dxi) .lt. tol) then
          noconv = .false.
        endif

        nc = nc + 1
        if(nc.gt.10) then
          noconv = .false.
          write(*,*) ' No convergence in 10 iterations ',xi0,dxi
        endif
      end do ! while

      end

      logical function cmast2(sg,x,xs,xm,ix2,ndm,nel,me,xi0,gn0,xjac,nn)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Find coordinate, gap, and master element number

c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'c_geom.h'

      integer  mel,nel,ndm, i,j,me,mn, mecl
      real*8   sg,xjac,denom,xi0,gn0, tol, ximin
      integer  ix2(dnope2,*)
      real*8   cxs(2),dxs(2),nn(2),x(ndm,*),xs(2,*)
      real*8   shp(2,4),dxm(2),ff(2),xm(2,*)

      data     tol / 1.d-8 /

      call shp1dn(sg,shp,nel)
      do j = 1,2
        cxs(j) = 0.0d0
        dxs(j) = 0.0d0
        do i = 1,nel
          cxs(j) = cxs(j) + shp(2,i)*xs(j,i)
          dxs(j) = dxs(j) + shp(1,i)*xs(j,i)
        end do ! i
      end do ! j
      xjac  =  sqrt(dxs(1)**2 + dxs(2)**2)
      nn(1) =  dxs(2)/xjac
      nn(2) = -dxs(1)/xjac

c     Search for correct master segment

      ximin  = 1.0d20
      mecl   = 1
      cmast2 = .false.
      do me = 1,neps2
        do j = 1,2
          xm(j,1) = x(j,ix2(1,me))
          xm(j,2) = x(j,ix2(2,me))
          ff(j)   = cxs(j) - xm(j,1)
          dxm(j)  = 0.5d0*(xm(j,2) - xm(j,1))
        end do ! j
        denom = nn(1)*dxm(2) - nn(2)*dxm(1)
        xi0   = (nn(1)*ff(2) - nn(2)*ff(1))/denom
        gn0   = (dxm(1)*ff(2) - dxm(2)*ff(1))/denom
        if(nope2.gt.2) then
          mel = 0
          do mn = 1,nope2
            if(ix2(mn,me).gt.0) then
              mel       = mel + 1
              xm(1,mel) = x(1,ix2(mn,me))
              xm(2,mel) = x(2,ix2(mn,me))
            endif
          end do ! mn
          call cgapt2(xi0,gn0, xm,cxs,nn, mel)
        endif
        if(xi0.gt.-tol .and. xi0.le.2.0d0+tol) then
          xi0    =  xi0 - 1.0d0
          cmast2 = .true.
          return
        else
          if(abs(xi0-1.0d0).lt.ximin) then
            ximin = abs(xi0-1.0d0)
            mecl  = me
          endif
        endif
      end do ! me

      me = mecl

      end

      logical function slavt2( cxm,xs,nel,xic )

      implicit   none

      logical    noconv
      integer    nel, i,j,count
      real*8     cxm(2), xc(2), xs(2,*)
      real*8     xic,dxic, gnn,dgnn, delta, tol
      real*8     rc(2), shp(3,5), tp(2), tc(2), tt(2), nn(2)

      data       tol / 1.d-9 /

      xic    = 0.0d0
      gnn    = 0.0d0
      noconv = .true.
      count  = 0
      do while (noconv)
        call shp1dt(xic,shp,nel)

        do i = 1,2
          xc(i) = 0.0d0
          tc(i) = 0.0d0
          tt(i) = 0.0d0
          do j = 1,nel
            tc(i) = tc(i) + shp(3,j)*xs(i,j)
            xc(i) = xc(i) + shp(2,j)*xs(i,j)
            tt(i) = tt(i) + shp(1,j)*xs(i,j)
          end do ! j
        end do ! i
        nn(1) =  tt(2)
        nn(2) = -tt(1)
        tp(1) =  tt(1) + gnn*tc(2)
        tp(2) =  tt(2) - gnn*tc(1)

        rc(1) = cxm(1) - xc(1) - nn(1)*gnn
        rc(2) = cxm(2) - xc(2) - nn(2)*gnn
        delta = 1.d0/(tt(1)*tp(1) + tt(2)*tp(2))
        dgnn  = (tp(2)*rc(1) - tp(1)*rc(2))*delta
        dxic  = (tt(1)*rc(1) + tt(2)*rc(2))*delta

        gnn   = gnn + dgnn
        xic   = xic + dxic

        if(abs(dxic).lt.tol) then
          noconv = .false.
        endif
        count = count + 1
        if(count.gt.25) then
          noconv = .false.
        endif
      end do ! while

c     Keep within Slave facet

      if(xic.lt.-1.d0) then
        xic = -1.d0
      endif
      if(xic.gt. 1.d0) then
        xic =  1.d0
      endif

      slavt2 = .not.noconv

      end

      logical function cmasp2(sg,x,xs,xm,ix2,ndm,nel,me,xi0,gn0,xjac,nn)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Find coordinate, gap, and master element number

c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'c_geom.h'

      integer  mel,nel,ndm, i,j,me,mn
      real*8   sg,xjac,denom,xi0,gn0, tol
      integer  ix2(dnope2,*)
      real*8   cxs(2),dxs(2),nn(2),x(ndm,*),xs(2,*)
      real*8   shp(2,4),dxm(2),ff(2),xm(2,*)

      data     tol / 1.d-8 /

      call shp1dn(sg,shp,nel)
      do j = 1,2
        cxs(j) = 0.0d0
        dxs(j) = 0.0d0
        do i = 1,nel
          cxs(j) = cxs(j) + shp(2,i)*xs(j,i)
          dxs(j) = dxs(j) + shp(1,i)*xs(j,i)
        end do ! i
      end do ! j
      xjac  =  sqrt(dxs(1)**2 + dxs(2)**2)
      nn(1) =  dxs(2)/xjac
      nn(2) = -dxs(1)/xjac

c     Set master segment geometry

      cmasp2 = .false.
      do j = 1,2
        xm(j,1) = x(j,ix2(1,me))
        xm(j,2) = x(j,ix2(2,me))
        ff(j)   = cxs(j) - xm(j,1)
        dxm(j)  = 0.5d0*(xm(j,2) - xm(j,1))
      end do ! j
      denom = nn(1)*dxm(2) - nn(2)*dxm(1)
      xi0   = (nn(1)*ff(2) - nn(2)*ff(1))/denom
      gn0   = (dxm(1)*ff(2) - dxm(2)*ff(1))/denom
      if(nope2.gt.2) then
        mel = 0
        do mn = 1,nope2
          if(ix2(mn,me).gt.0) then
            mel       = mel + 1
            xm(1,mel) = x(1,ix2(mn,me))
            xm(2,mel) = x(2,ix2(mn,me))
          endif
        end do ! mn
        call cgapt2(xi0,gn0, xm,cxs,nn, mel)
      endif
      if(xi0.gt.-tol .and. xi0.le.2.0d0+tol) then
        cmasp2 = .true.
      endif
      xi0 =  xi0 - 1.0d0

      end
