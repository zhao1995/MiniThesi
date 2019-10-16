c$Id:$
      subroutine crigplt(cs0,ch2,cn,xs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set plot value for rigid/deformable surface node

c     Inputs:
c        cs0(nr0,*) - Surface prameter table
c        ch2(*)     - History variables: Augmented or Lagrange force
c        cn         - Surface penalty
c        xs(3)      - Coordinate

c     Outputs:
c        cpl(1)     - Plot value
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_comnd.h'
      include   'c_keyh.h'
      include   'c_pair.h'
      include   'prld1.h'
      include   'prlod.h'
      include   'ptdat4.h'
      include   'sdata.h'
      include   'trdata.h'

      integer    i,np, nsopt, dn
      real*8     cs0(nr0,n0c1:*),ch2(*),xs(3)
      real*8     xx(3), gn, rr,un, fn, cn

c     Extract surface type and parameters

      nsopt = nint(cs0(2,0))
      dn    = nint(cs0(3,0))
      rr    = cs0(4,0)
      un    = cs0(5,0)
      np    = nint(cs0(6,0))
      gn    = 0.0d0

      if(np.gt.0) then
        un = un*prldv(np)
      else
        un = un*prop
      endif

c     Cylindrical or Spherical surface

      if(nsopt.le.2) then

        do i = 1,ndm
          xx(i) = tr(1,i)*(xs(1) - xr(1))
     &          + tr(2,i)*(xs(2) - xr(2))
     &          + tr(3,i)*(xs(3) - xr(3))
        end do ! i
        if(nsopt.eq.1) then
          xx(3) = 0.0d0
        endif

c       Compute gap function

        if(dn.gt.0) then
          gn = sqrt(xx(1)**2 + xx(2)**2 + xx(3)**2) - rr - un
        else
          gn = rr + un - sqrt(xx(1)**2 + xx(2)**2 + xx(3)**2)
        endif

c     Cartesian surface

      elseif(nsopt.eq.3) then

        if(dn.gt.0) then
          gn = xs(dn) - rr - un
        else
          gn = rr + un - xs(-dn)
        endif

      endif ! nsopt

c     Set contact force

      if(ifaugm.le.1) then
        fn = cn*gn
      else
        fn = cn*gn + ch2(p1(151))
      endif

c     Lagrange multiplier

      if(ifsolm.eq.2) then
        fn = fn + ch2(p1(21))
      endif

c     Test on contact force and accumlate contact forces

      if(fn.lt.0.0d0) then
        cpl(1) = cpl(1) + abs(fn)
      endif ! g

      end
