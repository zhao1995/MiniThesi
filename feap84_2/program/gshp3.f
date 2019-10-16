c$Id:$
      subroutine gshp3(xi, ixl, shp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Shape functions for 3-d mesh generation by block.
c               8 to 27 nodes.  No derivatives computed.

c      Inputs:
c         xi(3)   - Natural coordinates for point
c         ixl(*)  - List of active nodes

c      Outputs:
c         shp(*)  - Shape functions for point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, j, nxi(3,27),ixl(27)
      real*8    f1,f2,f3, xi(3),shp(27)

      save

      data      nxi/-1,-1,-1,   1,-1,-1,   1, 1,-1, -1, 1,-1,
     &              -1,-1, 1,   1,-1, 1,   1, 1, 1, -1, 1, 1,
     &              -1,-1, 0,   1,-1, 0,   1, 1, 0, -1, 1, 0,
     &               0,-1,-1,   1, 0,-1,   0, 1,-1, -1, 0,-1,
     &               0, 0,-1,
     &               0,-1, 1,   1, 0, 1,   0, 1, 1, -1, 0, 1,
     &               0, 0, 1,
     &               0,-1, 0,   1, 0, 0,   0, 1, 0, -1, 0, 0,
     &               0, 0, 0/

c     Generate first functions

      do i = 1,4
        f1 = 0.25*(1.+nxi(1,i)*xi(1))*(1.+nxi(2,i)*xi(2))
        f2 = 1.-xi(3)
        f3 = 1.+xi(3)
        shp(i  ) = 0.50*f1*f2
        shp(i+4) = 0.50*f1*f3
        shp(i+8) =      f1*f2*f3
      end do ! i

c     Form relative shape functions - level 1

      f1      = 0.25*(1.-xi(1)**2)
      if(ixl(13).ne.0)
     & shp(13) = f1*(1.+nxi(2,13)*xi(2))*(1.+nxi(3,13)*xi(3))
      if(ixl(15).ne.0)
     & shp(15) = f1*(1.+nxi(2,15)*xi(2))*(1.+nxi(3,15)*xi(3))
      if(ixl(18).ne.0)
     & shp(18) = f1*(1.+nxi(2,18)*xi(2))*(1.+nxi(3,18)*xi(3))
      if(ixl(20).ne.0)
     & shp(20) = f1*(1.+nxi(2,20)*xi(2))*(1.+nxi(3,20)*xi(3))
      f1      = 0.25*(1.-xi(2)**2)
      if(ixl(14).ne.0)
     & shp(14) = f1*(1.+nxi(1,14)*xi(1))*(1.+nxi(3,14)*xi(3))
      if(ixl(16).ne.0)
     & shp(16) = f1*(1.+nxi(1,16)*xi(1))*(1.+nxi(3,16)*xi(3))
      if(ixl(19).ne.0)
     & shp(19) = f1*(1.+nxi(1,19)*xi(1))*(1.+nxi(3,19)*xi(3))
      if(ixl(21).ne.0)
     & shp(21) = f1*(1.+nxi(1,21)*xi(1))*(1.+nxi(3,21)*xi(3))

c     Form relative shape functions - level 2

      f1      = 1.-xi(1)**2
      f2      = 1.-xi(2)**2
      f3      = 1.-xi(3)**2

      if(ixl(17).ne.0)
     & shp(17) = 0.50*f1*f2*(1.+nxi(3,17)*xi(3))
      if(ixl(22).ne.0)
     & shp(22) = 0.50*f1*f2*(1.+nxi(3,22)*xi(3))
      if(ixl(23).ne.0)
     & shp(23) = 0.50*f3*f1*(1.+nxi(2,23)*xi(2))
      if(ixl(24).ne.0)
     & shp(24) = 0.50*f2*f3*(1.+nxi(1,24)*xi(1))
      if(ixl(25).ne.0)
     & shp(25) = 0.50*f3*f1*(1.+nxi(2,25)*xi(2))
      if(ixl(26).ne.0)
     & shp(26) = 0.50*f2*f3*(1.+nxi(1,26)*xi(1))

c     Form relative shape functions - level 3

c     Convert to absolute shape functions

c     (1) Modify for center node

      if(ixl(27).ne.0) then
        shp(27) = f1*f2*f3
        f1 = 0.125d0*shp(27)
        f2 = 0.250d0*shp(27)
        do i = 1,4
          shp(i  ) = shp(i  ) - f1
          shp(i+4) = shp(i+4) - f1
          if(ixl(i+8 ).ne.0) shp(i+8 ) = shp(i+8 ) - f2
          if(ixl(i+18).ne.0) shp(i+18) = shp(i+18) - f2
          if(ixl(  17).ne.0) shp(  17) = shp(  17) - f2
          if(ixl(  26).ne.0) shp(  26) = shp(  26) - f2
        end do ! i
      endif

c     (2) Modify for mid-face nodes

      if(ixl(17).ne.0) then
        f1 = 0.5d0*shp(17)
        if(ixl(13).ne.0) shp(13) = shp(13) - f1
        if(ixl(14).ne.0) shp(14) = shp(14) - f1
        if(ixl(15).ne.0) shp(15) = shp(15) - f1
        if(ixl(16).ne.0) shp(16) = shp(16) - f1
        f1 = 0.5d0*f1
        shp(1) = shp(1) - f1
        shp(2) = shp(2) - f1
        shp(3) = shp(3) - f1
        shp(4) = shp(4) - f1
      endif

      if(ixl(22).ne.0) then
        f1 = 0.5d0*shp(22)
        if(ixl(18).ne.0) shp(18) = shp(18) - f1
        if(ixl(19).ne.0) shp(19) = shp(19) - f1
        if(ixl(20).ne.0) shp(20) = shp(20) - f1
        if(ixl(21).ne.0) shp(21) = shp(21) - f1
        f1 = 0.5d0*f1
        shp(5) = shp(5) - f1
        shp(6) = shp(6) - f1
        shp(7) = shp(7) - f1
        shp(8) = shp(8) - f1
      endif

      do i = 1,4
        j  = (mod(i,4)) + 1
        if(ixl(22+i).ne.0) then
          f1 = 0.5d0*shp(22+i)
          if(ixl(i+8 ).ne.0) shp(i+8 ) = shp(i+8 ) - f1
          if(ixl(j+8 ).ne.0) shp(j+8 ) = shp(j+8 ) - f1
          if(ixl(i+12).ne.0) shp(i+12) = shp(i+12) - f1
          if(ixl(i+17).ne.0) shp(i+17) = shp(i+17) - f1
          f1 = 0.5d0*f1
          shp(i  ) = shp(i  ) - f1
          shp(i  ) = shp(i  ) - f1
          shp(i+4) = shp(i+4) - f1
          shp(j+4) = shp(j+4) - f1
        endif

c       (3) Modify for mid-edge nodes

        if(ixl( 8+i).ne.0) then
          f1 = 0.25d0*shp( 8+i)
          if(ixl(i+22).ne.0) shp(i+22) = shp(i+22) - f1
          shp(i  ) = shp(i  ) - f1
          shp(i+4) = shp(i+4) - f1
        endif

        if(ixl( 8+j).ne.0) then
          f1 = 0.25d0*shp( 8+j)
          if(ixl(i+22).ne.0) shp(i+22) = shp(i+22) - f1
          shp(j  ) = shp(j  ) - f1
          shp(j+4) = shp(j+4) - f1
        endif

        if(ixl(12+i).ne.0) then
          f1 = 0.5d0*shp(12+i)
          if(ixl(i+22).ne.0) shp(i+22) = shp(i+22) - f1
          shp(i  ) = shp(i  ) - f1
          shp(j  ) = shp(j  ) - f1
        endif

        if(ixl(17+i).ne.0) then
          f1 = 0.5d0*shp(17+i)
          if(ixl(i+22).ne.0) shp(i+22) = shp(i+22) - f1
          shp(i+4) = shp(i+4) - f1
          shp(j+4) = shp(j+4) - f1
        endif
      end do ! i

      end
