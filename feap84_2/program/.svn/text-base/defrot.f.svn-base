c$Id:$
      subroutine defrot(xln,f,id,numnp,ndf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'dabs' to 'abs'                           17/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set default vectors in rotation matrix, if not defined

c      Programmed by: M.S. Rifai/JCS for rotation matrix initialization

c      Inputs:
c         f(ndf,*)  - Input values defining boundary rotation
c         id(ndf,*) - Equation numbers for dof's
c         numnp     - Number of nodal points in mesh
c         ndf       - Number dof/node

c      Outputs:
c         xln(*,-)  - Rotation array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, numnp, ndf
      integer   id(ndf,*)
      real*8    t3t3, t1t3, t2t3, t1t1, t2t2, fac, tnrm, dot
      real*8    xln(3,3,6,numnp),f(ndf,*)

      save

      do i = 1 , numnp

c       If t3 is 0, set to old value

        t3t3 = dot(xln(1,3,1,i),xln(1,3,1,i),3)
        if (t3t3.lt.1.d-5) then
          xln(1,3,1,i) = xln(1,3,2,i)
          xln(2,3,1,i) = xln(2,3,2,i)
          xln(3,3,1,i) = xln(3,3,2,i)
        endif

c       Orthogonalize t1 wrt t

        t1t3 = dot(xln(1,1,1,i),xln(1,3,1,i),3)
        if (abs(t1t3).gt.1.d-5) then
          xln(1,1,1,i) = xln(1,1,1,i) - t1t3*xln(1,3,1,i)
          xln(2,1,1,i) = xln(2,1,1,i) - t1t3*xln(2,3,1,i)
          xln(3,1,1,i) = xln(3,1,1,i) - t1t3*xln(3,3,1,i)
        endif

c       Orthogonalize t2 wrt t

        t2t3 = dot(xln(1,2,1,i),xln(1,3,1,i),3)
        if (abs(t2t3).gt.1.d-5) then
          xln(1,2,1,i) = xln(1,2,1,i) - t2t3*xln(1,3,1,i)
          xln(2,2,1,i) = xln(2,2,1,i) - t2t3*xln(2,3,1,i)
          xln(3,2,1,i) = xln(3,2,1,i) - t2t3*xln(3,3,1,i)
        endif

c       If t1 is 0, set to t2*t3

        t1t1 = dot(xln(1,1,1,i),xln(1,1,1,i),3)
        if (t1t1.lt.1.d-5) then
          xln(1,1,1,i) = xln(2,2,1,i)*xln(3,3,1,i)
     &                 - xln(3,2,1,i)*xln(2,3,1,i)
          xln(2,1,1,i) = xln(3,2,1,i)*xln(1,3,1,i)
     &                 - xln(1,2,1,i)*xln(3,3,1,i)
          xln(3,1,1,i) = xln(1,2,1,i)*xln(2,3,1,i)
     &                 - xln(2,2,1,i)*xln(1,3,1,i)
        endif

c       If t2 is 0, set to t3*t1

        t2t2 = dot(xln(1,2,1,i),xln(1,2,1,i),3)
        if (t2t2.lt.1.d-5) then
          xln(1,2,1,i) = xln(2,3,1,i)*xln(3,1,1,i)
     &                 - xln(3,3,1,i)*xln(2,1,1,i)
          xln(2,2,1,i) = xln(3,3,1,i)*xln(1,1,1,i)
     &                 - xln(1,3,1,i)*xln(3,1,1,i)
          xln(3,2,1,i) = xln(1,3,1,i)*xln(2,1,1,i)
     &                 - xln(2,3,1,i)*xln(1,1,1,i)
        endif

c       If t1 or t2 is still 0 set to old value

        t1t1 = dot(xln(1,1,1,i),xln(1,1,1,i),3)
        t2t2 = dot(xln(1,2,1,i),xln(1,2,1,i),3)
        if (t1t1.lt.1.d-5 .or. t2t2.lt.1.d-5) then
          if (xln(3,3,1,i).gt.0.d0) then
            fac          = 1.d0 / ( 1.d0 + xln(3,3,1,i) )
            xln(1,1,1,i) = xln(3,3,1,i)+fac*xln(2,3,1,i)*xln(2,3,1,i)
            xln(1,2,1,i) =             -fac*xln(2,3,1,i)*xln(1,3,1,i)
            xln(2,1,1,i) =             -fac*xln(1,3,1,i)*xln(2,3,1,i)
            xln(2,2,1,i) = xln(3,3,1,i)+fac*xln(1,3,1,i)*xln(1,3,1,i)
            xln(3,1,1,i) =-xln(1,3,1,i)
            xln(3,2,1,i) =-xln(2,3,1,i)
          else
            fac          = 1.d0 / ( 1.d0 - xln(3,3,1,i) )
            xln(1,1,1,i) =-xln(3,3,1,i)+fac*xln(2,3,1,i)*xln(2,3,1,i)
            xln(1,2,1,i) =              fac*xln(2,3,1,i)*xln(1,3,1,i)
            xln(2,1,1,i) =             -fac*xln(1,3,1,i)*xln(2,3,1,i)
            xln(2,2,1,i) = xln(3,3,1,i)-fac*xln(1,3,1,i)*xln(1,3,1,i)
            xln(3,1,1,i) = xln(1,3,1,i)
            xln(3,2,1,i) =-xln(2,3,1,i)
          endif
        endif

c       Normalize t1

        t1t1 = dot(xln(1,1,1,i),xln(1,1,1,i),3)
        if (abs(t1t1-1.d0).gt.1.d-8) then
          tnrm         = sqrt(t1t1)
          xln(1,1,1,i) = xln(1,1,1,i) / tnrm
          xln(2,1,1,i) = xln(2,1,1,i) / tnrm
          xln(3,1,1,i) = xln(3,1,1,i) / tnrm
          if (id(4,i).le.0) f(4,i) = f(4,i) / tnrm
        endif

c       Normalize t2

        t2t2 = dot(xln(1,2,1,i),xln(1,2,1,i),3)
        if (abs(t2t2-1.d0).gt.1.d-8) then
          tnrm         = sqrt(t2t2)
          xln(1,2,1,i) = xln(1,2,1,i) / tnrm
          xln(2,2,1,i) = xln(2,2,1,i) / tnrm
          xln(3,2,1,i) = xln(3,2,1,i) / tnrm
          if (id(5,i).le.0) f(5,i) = f(5,i) / tnrm
        endif

c       Normalize t3

        t3t3 = dot(xln(1,3,1,i),xln(1,3,1,i),3)
        if (abs(t3t3-1.d0).gt.1.d-8) then
          tnrm         = sqrt(t3t3)
          xln(1,3,1,i) = xln(1,3,1,i) / tnrm
          xln(2,3,1,i) = xln(2,3,1,i) / tnrm
          xln(3,3,1,i) = xln(3,3,1,i) / tnrm
        endif
      end do ! i

      end
