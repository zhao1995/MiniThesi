c$Id:$
      subroutine dmat2d(d,psi,dmg,betag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct definition of qm(5,5) to qm(6,6)         01/09/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Rotation of material arrays from principal to
c                local element directions

c      Inputs:
c         d        - Array with material properties
c         psi      - Angle from y1-axis (local) to 1-axis (principal)

c      Outputs:
c         dmg(6,6) - Plane modulus matrix
c         betag(6) - global thermal stress/temperature

c-----[--.----+----.----+----.-----------------------------------------]
c      Variables used in subroutine

c         qm(6,6)  - Transformation matrix

c         dml(6,6) - Local (orthotropic ) modulus matrix
c         dmlqj(6) - intermediate matrix for triple product

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, j, k, imat, nsiz
      real*8    psi, si, co, s2, c2, cs, d(*)
      real*8    dml(6,6),betal(6), qm(6,6),dmlqj(6), dmg(6,6),betag(6)

      save

c     Assign material properties

      imat = nint(d(20))

c     No rotation

      if(psi.eq.0.0d0) then

        do i = 1,6
          betag(i) = 0.0d0
          do j = 1,6
            dmg(j,i) = 0.0d0
          end do ! j
        end do ! i

        if(imat.eq.8) then

          nsiz = nint(d(200))
          k    = 200
          do i = 1,nsiz
            do j = 1,i
              k        = k + 1
              dmg(i,j) = d(k)
              dmg(j,i) = d(k)
            end do ! j
          end do ! i

        else

          dmg(1,1) = d(21)
          dmg(2,2) = d(22)
          dmg(3,3) = d(23)

          dmg(1,2) = d(24)
          dmg(2,1) = dmg(1,2)

          dmg(2,3) = d(25)
          dmg(3,2) = dmg(2,3)

          dmg(3,1) = d(26)
          dmg(1,3) = dmg(3,1)

          dmg(4,4) = d(27)
          dmg(5,5) = d(28)
          dmg(6,6) = d(29)

          betag(1) = d(47)
          betag(2) = d(48)
          betag(3) = d(49)

        endif ! no rotation

c     Orthotropic material (rotations)

      else

c       Zero arrays

        do i = 1,6
          betal(i) = 0.0d0
          do j = 1,6
            dml(j,i) = 0.0d0
            qm(j,i)  = 0.0d0
          end do ! j
        end do ! i

c       Set up constants for transformation

        si = sin(psi)
        co = cos(psi)
        s2 = si*si
        c2 = co*co
        cs = co*si

c       Set up transformation matrix for plane problem

        qm(1,1) =  c2
        qm(1,2) =  s2
        qm(1,4) =  cs

        qm(2,1) =  s2
        qm(2,2) =  c2
        qm(2,4) = -cs

        qm(3,3) =  1.0d0

        qm(4,1) = -2.d0 * cs
        qm(4,2) =  2.d0 * cs
        qm(4,4) =  c2 - s2

        qm(5,5) =  co
        qm(5,6) = -si

        qm(6,5) =  si
        qm(6,6) =  co

        if(imat.eq.8) then

          nsiz = nint(d(200))
          k    = 200
          do i = 1,nsiz
            do j = 1,i
              k        = k + 1
              dml(i,j) = d(k)
              dml(j,i) = d(k)
            end do ! j
          end do ! i

c       Set up local (orthotropic) plane matrix

        else

          dml(1,1) = d(21)
          dml(2,2) = d(22)
          dml(3,3) = d(23)

          dml(1,2) = d(24)
          dml(2,1) = dml(1,2)

          dml(2,3) = d(25)
          dml(3,2) = dml(2,3)

          dml(3,1) = d(26)
          dml(1,3) = dml(3,1)

          dml(4,4) = d(27)
          dml(5,5) = d(28)
          dml(6,6) = d(29)

          betal(1) = d(47)
          betal(2) = d(48)
          betal(3) = d(49)

        endif

c       Convert plane local to global matrix

        do j = 1,6 ! {

c       Global thermal vector

          betag(j) = qm(1,j)*betal(1) + qm(2,j)*betal(2)
     &             + qm(3,j)*betal(3) + qm(4,j)*betal(4)
     &             + qm(5,j)*betal(5) + qm(6,j)*betal(6)

c         Global moduli array

          do i = 1,6
            dmlqj(i) = dml(i,1)*qm(1,j) + dml(i,2)*qm(2,j)
     &               + dml(i,3)*qm(3,j) + dml(i,4)*qm(4,j)
     &               + dml(i,5)*qm(5,j) + dml(i,6)*qm(6,j)
          end do ! i

          do i = 1,6 ! {
            dmg(i,j) = qm(1,i)*dmlqj(1) + qm(2,i)*dmlqj(2)
     &               + qm(3,i)*dmlqj(3) + qm(4,i)*dmlqj(4)
     &               + qm(5,i)*dmlqj(5) + qm(6,i)*dmlqj(6)
          end do ! i   }

        end do ! j   }

      endif

      end
