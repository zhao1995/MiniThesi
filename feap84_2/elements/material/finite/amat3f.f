c$Id:$
      subroutine amat3f(tautil,atilp,q, tau,aa)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Tangent Moduli: Constitutive Model in Principal Stretch

c      INPUT variables
c         tautil(3)  Principal deviatoric stresses
c         atilp(6,6) Deviatoric tangent matrix in principal basis
c         q(6,6)     Transformation matrix: [principal]-->[std]

c      OUTPUT variables
c         tau(6)     Kirchhoff stress tensor  (deviatoric)
c         aa(6,6)    Kirchhoff tangent moduli (deviatoric)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, j, l
      real*8    tautil(3), atilp(6,6), q(6,6), tau(6), aa(6,6),vi(6)

      save

c     Compute stress tensor

      do j = 1,6
        tau(j) = q(j,1)*tautil(1) + q(j,2)*tautil(2) + q(j,3)*tautil(3)
      end do ! j

c     Compute tangent matrix in standard basis

      do i = 1,6

c       Left multiplication: V = Q * AL

        do l = 1,3
          vi(l) = 0.0d0
          do j = 1,3
            vi(l) = vi(l) + q(i,j)*atilp(j,l)
          end do ! j
        end do ! l

        do l = 4,6
            vi(l) = q(i,l)*atilp(l,l)
        end do ! l

c       Right multiplication: AA = V * Q_t

        do j = 1,i
          aa(i,j) = 0.0d0
          do l = 1,6
            aa(i,j) = aa(i,j) + vi(l)*q(j,l)
          end do ! l

c         Fill in symmetric part

          aa(j,i) = aa(i,j)

        end do ! j

      end do ! i

      end
