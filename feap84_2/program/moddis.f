c$Id:$
      subroutine moddis(xl,ul,ndm,ndf,nel,nen,rcg,rlam,rben)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Computes modal displacements at t_n+a
c              for stress outputs assuming
c              x_n+a = r_n+a + lam_n+a(X-R)
c                            + alph*lam_n+1*wf_n+1
c                            + (1-alph)*lam_n*wf_n
c              Added on 2/97 by Alecia J. Chen

c     Note:    This function is incomplete
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      include  'ddata.h'

      integer  ndm,ndf,nel,nen,rben, i,j,n,nrd
      real*8   xl(ndm,*),ul(ndf,nen,*),rcg(3,11,*),rlam(9,6,*), vh(3)
      real*8   lamn(3,3),lam1(3,3),y1(3),v1(3),yn(3),vn(3),bigy(3),alpn

      save

      if(nrk.eq.0) then
        nrd = 2
      else
        nrd = 7
      endif

      alpn = 1.d0 - theta(3)
      call quamat(rlam(1,1,rben),lamn)
      call quamat(rlam(1,3,rben),lam1)

      do n = 1,nel

c       Compute (x - r) at t_n, t_n+1, t_0

        do i = 1,ndm
          y1(i) = xl(i,n) + ul(i,n,1) + alpn*ul(i,n,2) - rcg(i,2,rben)
          yn(i) = y1(i) - ul(i,n,2) + rcg(i,3,rben)
          bigy(i)  = xl(i,n) - rcg(i,1,rben)
        end do ! i

        do i = 1,3
          vh(i) = 0.0d0
          v1(i) = 0.0d0
          vn(i) = 0.0d0
        end do ! i
        do i = 1,ndm
          do j = 1,ndm
            v1(i) = v1(i) + lam1(j,i)*y1(j)
            vn(i) = vn(i) + lamn(j,i)*yn(j)
            vh(i) = vh(i)+0.25d0*(lam1(j,i)+lamn(j,i))*(y1(j)+yn(j))
          end do ! j
        end do ! i

        do i = 1,ndm
          ul(i,n,1) = theta(3)*v1(i) + alpn*vn(i) - bigy(i)
        end do ! i

      end do ! n

      end
