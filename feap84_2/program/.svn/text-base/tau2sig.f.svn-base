c$Id:$
      subroutine tau2sig(tau,ctau,volmr,detf, sig,dd, ndm,ntau)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Convert Kirchhoff quatitites to Cauchy ones

c     Inputs:
c       tau(6)    - Kirchhoff stress
c       ctau(6,6) - Kirchhoff moduli
c       volmr     - Volume of RVE cell
c       detf      - Deformation gradient determinant
c       ndm       - Mesh spatial dimension

c     Outputs:
c       sig(6)    - Cauchy stress
c       dd(6,6)   - Spatial moduli
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit  none

      integer   ndm,ntau, a,b, i,j,k,ik
      integer   ismap(3,3)
      real*8    tau(*),ctau(ntau,*), volmr,detf,jinv, sig(*),dd(ntau,*)

      save

      data        ismap / 1,4,6,  4,2,5,  6,5,3 /

c     Transform Kirchhoff to Cauchy stress and spatial moduli

      do a = 1,ntau
        do b = 1,ntau
          dd(b,a)   = 0.0d0
        end do ! b
      end do ! a

      do k = 1,3
        do j = 1,3
          b  = ismap(j,k)
          do i = 1,3
            ik      = ismap(i,k)
            a       = ismap(i,j)
            dd(a,b) = dd(a,b) - tau(ik)
          end do ! k
        end do ! i
      end do ! j

c     Transform Kirchhoff stress to Cauchy stress

      do a = 4,6
        do b = 1,6
          dd(a,b) = dd(a,b)*0.5d0
          dd(b,a) = dd(b,a)*0.5d0
        end do ! b
      end do ! a

c     Divide by deformation gradient

      jinv = volmr/detf
      do b = 1,6
        sig(b) = tau(b)*jinv
        do a = 1,6
          dd(a,b) = (ctau(a,b) + dd(a,b))*jinv
        end do ! a
      end do ! b

c     Coupled case

      if(ntau.gt.6) then
        do a = 1,6
          dd(a,7) = ctau(a,7)*jinv
          dd(7,a) = ctau(7,a)*jinv
        end do ! a
        dd(7,7) = ctau(7,7)*jinv
      endif

c     For 2-d problems the stress has been averaged before

      if(ndm.eq.2) then
        sig(3) = tau(3)
      endif

      end
