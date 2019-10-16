c$Id:$
      subroutine cvisco(d,eps,sig,dr,di)


c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Form complex moduli and real stress parts

c     Inputs:
c       d(*)    - Material parameters
c       eps(*)  - Real strain

c     Outputs:
c       sig(*)  - Real Stress
c       dr(*,*) - Real moduli
c       di(*,*) - Imaginary moduli
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'pconstant.h'
      include   'sdata.h'
      include   'tdata.h'

      integer    i,j, n
      real*8     G,Gr,Gi,K,Kr,Ki, mu_0,mu_n,mu_r,mu_i, om_n,de_n,theta
      real*8     d(*),eps(*),sig(*),ee(6),dr(6,6),di(6,6),one(6)

      save

      data       one / 3*1.d0 , 3*0.0d0 /

c     Set elastic parameters for G (mu) and lambda

      G     = d(1)/(2.d0*(1.d0 + d(2)))
      K     = d(1)/(3.d0*(1.d0 - 2.d0*d(2)))


c     Compute volumetric strain and deviatoric components

      theta = (eps(1) + eps(2) + eps(3))*one3
      do i = 1,3
        ee(i  ) = eps(i) - theta
        ee(i+3) = eps(i+3)*0.5d0
      end do ! i

c     Set properties

      mu_0 = 1.0d0
      mu_r = 0.0d0
      mu_i = 0.0d0
      do n = 1,nint(d(57))
        om_n = ttim*d(2*n+50)
        mu_n = d(2*n+49)
        mu_0 = mu_0 - mu_n
        de_n = 1.d0/(1.d0 + om_n*om_n)
        mu_r = mu_r + mu_n*(om_n*om_n)*de_n
        mu_i = mu_i + mu_n*om_n*de_n
      end do ! n

c     Finish updates and save the strains

      Gr = G*(mu_r + mu_0)
      Gi = G* mu_i
      Kr = K*theta*3.0d0
      do i = 1,6
        sig(i) = 2.d0*Gr*ee(i) + Kr*one(i)
      end do ! i

c     Set tangent parameters

      Kr = K - two3*Gr
      Ki =   - two3*Gi
      do j =1,3
        do i = 1,3
          dr(i,j) = Kr
          di(i,j) = Ki
        end do ! i
        dr(j  ,j  ) = dr(j,j) + 2.d0*Gr
        di(j  ,j  ) = di(j,j) + 2.d0*Gi
        dr(j+3,j+3) = Gr
        di(j+3,j+3) = Gi
      end do ! i

      end
