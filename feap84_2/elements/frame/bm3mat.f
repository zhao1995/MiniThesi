c$Id:$
      subroutine bm3mat(d,hn,h1,nh,eps,sig,cc,isw)
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Purpose:    Computes stress resultants and elasticity matrix

c     Input:      d       : Material data
c                 eps     : Strain [e_x,g_xy, g_yz, k_y, k_z, the]
c                 hn(*)   : History variables at t_n
c                 nh      : Number of history variables/point
c                 isw     : Switch value from element

c     Output:     sig     : Stress resultants [N , M]
c                 cc      : Elasticity tensor
c                 h1(*)   : History variables at t_n+1

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      implicit   none

      integer    nh, isw,i,j
      real*8     d(*),hn(*),h1(*),eps(6),sig(6),cc(6,6),g

c     Resultant elastic model

      if(nint(d(100)).eq.0) then

c       Compute material moduli for resultants

        do i = 1,6
          do j = 1,6
            cc(j,i) = 0.0d0
          end do ! j
        end do ! i

        g       = 0.5d0*d(1)/(1.d0 + d(2))
        cc(1,1) = d(37)*d(32)*g
        cc(2,2) = d(38)*d(32)*g
        cc(3,3) = d(1)*d(32)          ! E A
        cc(4,4) = d(1)*d(33)          ! E I_y
        cc(4,5) = d(1)*d(35)          ! E I_yz
        cc(5,4) = d(1)*d(35)          ! E I_yz
        cc(5,5) = d(1)*d(34)          ! E I_z
        cc(6,6) = g*d(36)             ! J G

c       Compute stress resultants

        do i = 1,6
          sig(i) = cc(i,1)*eps(1) + cc(i,2)*eps(2)
     &           + cc(i,3)*eps(3) + cc(i,4)*eps(4)
     &           + cc(i,5)*eps(5) + cc(i,6)*eps(6)
        end do ! i

c     Integrated models

      else
        call bm3res(d,hn,h1,nh,eps,sig,cc,isw)
      endif

      end
