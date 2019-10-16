c$Id:$
      subroutine viscoe(d,ta,eps,en,qi,ntm,sig,dd,dr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'pconstant.h'                 14/11/2006
c       2. Take alpha from d(47) instead of d(3)            10/07/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Linear viscoelastic model (shear only)

c      Inputs:
c         d(*)    - Material parameters
c         ta      - Temperature
c         eps(*)  - Strains at t_n+1
c         en(*)   - Strains at t_n
c         qi(*)   - Viscoelastic strain
c         ntm     - Number of terms

c      Outputs:
c         sig(*)  - Stresses
c         dd(6,6) - Viscoelastic moduli
c         dr(6,6) - Rayleigh moduli
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'
      include  'sdata.h'
      include  'tdata.h'

      integer   i,j, n,nv,ntm
      real*8    G,Gg,K,Kg,Kth, gfac,exp_n,mu_0,mu_n,dq_n,dtau, theta
      real*8    alpha,ta, d(*),eps(*),en(*),qi(ntm,*),ee(6)
      real*8    sig(*),dd(6,6),dr(6,6), hvisc

      save

c     Set elastic parameters for G (mu) and lambda

      G     = d(1)/(2.d0*(1.d0 + d(2)))
      K     = d(1)/(3.d0*(1.d0 - 2.d0*d(2)))
      alpha = d(47)

c     Compute volumetric strain and deviatoric components

      theta = (eps(1) + eps(2) + eps(3))*one3
      do i = 1,3
        ee(i  ) = eps(i) - theta
      end do ! i
      do i = 4,ntm
        ee(i) = eps(i)*0.5d0
      end do ! i

c     Set properties for integrating the q_i terms

      do i = 1,ntm
        sig(i) = 0.0d0
      end do ! i
      mu_0 = 0.0d0
      gfac = 0.0d0

      nv   = nint(d(57))
      do n = 1,nv
        mu_n  = d(2*n+49)
        dtau  = dt/d(2*n+50)
        exp_n = exp(-dtau)

        dq_n = mu_n * hvisc(dtau,exp_n)
        gfac = gfac + dq_n
        mu_0 = mu_0 + mu_n

c       Update history and compute viscoelastic deviatoric stress

        do i = 1,ntm
          qi(i,n) = exp_n*qi(i,n) + dq_n*(ee(i) - en(i))
          sig(i)  = sig(i) + qi(i,n)
        end do ! i
      end do ! n

c     Finish updates and save the strains

      mu_0 = 1.d0 - mu_0
      gfac = gfac + mu_0
      do i = 1,ntm
        sig(i) = 2.d0*G*(mu_0*ee(i) + sig(i))
        en(i)  = ee(i)
      end do ! i

c     Add elastic bulk term

      Kth = K*(theta*3.0d0 - alpha*ta)
      do i = 1,3
        sig(i) = sig(i) + Kth
      end do ! i

c     Set tangent parameters

      Gg = G*gfac
      Kg = K - two3*Gg
      K  = K - two3*G
      do j =1,3
        do i = 1,3
          dd(i,j) = Kg
          dr(i,j) = K
        end do ! i
        dd(j,j) = dd(j,j) + 2.d0*Gg
        dr(j,j) = dr(j,j) + 2.d0*G
      end do ! i

      do i = 4,ntm
        dd(i,i) = Gg
        dr(i,i) = G
      end do ! i

      end
