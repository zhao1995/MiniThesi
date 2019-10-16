c$Id:$
      subroutine pblend3x(x1,nr1,ns1, x2,nr2,ns2, x3,nr3,ns3,
     &                    x4,nr4,ns4, x5,nr5,ns5, x6,nr6,ns6,
     &                    nty,x,ndm,iblend, nf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form coordinates for nodes in blended block

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer nr1,ns1, nr2,ns2, nr3,ns3, nr4,ns4, nr5,ns5, nr6,ns6
      integer nr,ns,nt,j,ndm,nf, iblend(*),nty(*)
      real*8  r,s,t, dr,ds,dt  , xl(3,8),sh(2,3),x(ndm,*)
      real*8  x1(3,0:nr1,0:ns1), x2(3,0:nr2,0:ns2), x3(3,0:nr3,0:ns3)
      real*8  x4(3,0:nr4,0:ns4), x5(3,0:nr5,0:ns5), x6(3,0:nr6,0:ns6)

      save

      dr = 2.d0/dble(iblend(3))
      ds = 2.d0/dble(iblend(1))
      dt = 2.d0/dble(iblend(2))

      do j = 1,ndm
        xl(j,1) = x5(j,  0,  0)
        xl(j,2) = x5(j,nr5,  0)
        xl(j,3) = x5(j,nr5,ns5)
        xl(j,4) = x5(j,  0,ns5)
        xl(j,5) = x6(j,  0,  0)
        xl(j,6) = x6(j,nr6,  0)
        xl(j,7) = x6(j,nr6,ns6)
        xl(j,8) = x6(j,  0,ns6)
      end do ! j

      nf = iblend(4) - 1
      t  = -1.d0
      do nt = 0,iblend(2)
        sh(1,3) = 0.5d0 - 0.5d0*t
        sh(2,3) = 0.5d0 + 0.5d0*t
        s = -1.d0
        do ns = 0,iblend(1)
          sh(1,2) = 0.5d0 - 0.5d0*s
          sh(2,2) = 0.5d0 + 0.5d0*s
          r = -1.d0
          do nr = 0,iblend(3)
            sh(1,1) = 0.5d0 - 0.5d0*r
            sh(2,1) = 0.5d0 + 0.5d0*r
            nf      = nf + 1
            nty(nf) = 0
            do j = 1,ndm
              x(j,nf) =(x1(j,ns,nt)*sh(1,1) + x2(j,ns,nt)*sh(2,1)
     &                + x3(j,nt,nr)*sh(1,2) + x4(j,nt,nr)*sh(2,2)
     &                + x5(j,nr,ns)*sh(1,3) + x6(j,nr,ns)*sh(2,3))*0.5d0
            end do ! j
            call xiso3d(r,s,t,xl,x(1,nf))
            r  = r  + dr
          end do ! nr
          s = s + ds
        end do ! ns
        t = t + dt
      end do ! nt

      end
