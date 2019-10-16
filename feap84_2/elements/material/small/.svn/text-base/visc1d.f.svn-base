c$Id:$
      subroutine visc1d(d,eps,epsn,q, sig,ee)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: One-dimensional viscoelastic model

c     Inputs:
c       d(*)    - Material parameters
c       epsn    - Strain value at t_n
c       eps     - Strain value at t_n+1
c       q(*)    - Viscoelastic stress components at t_n

c     Outputs
c       sig     - Stress at t_n+1
c       q(*)    - Viscoelastic stress components at t_n+1
c       ee(*)   - Tangent moduli
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'tdata.h'

      integer   nv, n
      real*8    efac,mu,mu_n,exp_n,dtau,dq, sig,eps,epsn, hvisc
      real*8    d(*),q(*), ee(*)

c     Set initial values

      ee(1) = d(1)
      ee(2) = d(1)
      sig   = 0.0d0
      efac  = 0.0d0
      mu    = 0.0d0

c     Accumulate viscoelastic terms

      nv    = nint(d(57))
      do n = 1,nv

        mu_n  = d(2*n+49)

c       Exact integral

        dtau  = dt/d(2*n+50)
        exp_n = exp( -dtau)
        dq    = mu_n*hvisc(dtau,exp_n)
        q(n)  = exp_n*q(n) + dq*(eps - epsn)

c       Mid-point integral

c       dtau  = dt/d(2*n+50)*0.5d0
c       exp_n = exp( -dtau)
c       dq    = mu_n*exp_n
c       q(n)  = exp_n*exp_n*q(n) + dq*(eps - epsn)

c       Accumulate terms

        sig   = sig  + q(n)
        mu    = mu   + mu_n
        efac  = efac + dq

      end do ! n

c     Set final values and add elastic component

      mu    = 1.d0 - mu
      sig   = ee(1)*(mu*eps + sig)
      ee(1) = ee(1)*(mu + efac)
      epsn  = eps

      end
