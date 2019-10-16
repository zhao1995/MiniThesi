c$Id:$
      subroutine fpstrs(sig,aa,g33,f33,c33,ict,ictmax, conv,finite)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Check that ict > 1 for stress being zero         26/06/2010
c          Forces recompute of stress with new moduli
c       2. Modify to include displacement gradient for      10/01/2012
c          finite deformation
c       3. Change convergence to check on f33               24/02/2012
c       4. Do 3 iterations unless moduli already o.k.       14/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Plane stress driver for material models

c     Inputs:
c       sig(*)       - Full Cauchy stress values at t_n+1
c       aa(*,*)      - Full Material moduli
c       g33          - G(c33,c33) for finite deformation
c                    - eps(c33) for small  deformation
c       f33          - F(c33,c33) for finite deformation
c                    - eps(c33) for small  deformation
c       c33          - Number of component to eliminate (can be 2 or 3)
c       ict          - Iteration counter
c       ictmax       - Maximum iteration
c       finite       - True if finite deformation problem

c     Outputs
c       sig(*)       - Plane Stress Cauchy stress values at t_n+1
c       aa(*,*)      - Plane Stress Matrial moduli
c       conv         - Converged if true
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      logical    conv, finite
      integer    c33, c22, ict, ictmax
      real*8     sig(*), aa(6,6), g33, dg33,f33, dsig, a3inv, tol

      save

      data       tol    / 1.d-10 /

c     Compute stress

      dsig  = sig(c33)

c     Set other direction

      if(c33.eq.3) then
        c22 = 2
      elseif(c33.eq.2) then
        c22 = 3
      endif

c     Return if either stress or modulus is zero

      if(aa(c33,c33).eq.0.0d0 .or. (dsig.eq.0.0d0 .and. ict.gt.2)) then
        conv = .true.
        return
      endif

c     Perform iteration

      a3inv =  1.d0/aa(c33,c33)
      if(finite .and. c33.eq.3) then
        dg33  = -dsig*a3inv*f33
      else
        dg33  = -dsig*a3inv
      endif
      g33   =  g33 + dg33
      f33   =  f33 + dg33

      if(abs(dg33).le.tol*abs(f33)) then
        conv = .true.
      else
        conv = .false.
      endif

c     Reduce tangent array when converged

      if(conv .or. ict.ge.ictmax) then
        aa(1  ,1  ) = aa(1  ,1  ) - aa(1  ,c33)*a3inv*aa(c33,1  )
        aa(c22,1  ) = aa(c22,1  ) - aa(c22,c33)*a3inv*aa(c33,1  )
        aa(4  ,1  ) = aa(4  ,1  ) - aa(4  ,c33)*a3inv*aa(c33,1  )
        aa(1  ,c22) = aa(1  ,c22) - aa(1  ,c33)*a3inv*aa(c33,c22)
        aa(c22,c22) = aa(c22,c22) - aa(c22,c33)*a3inv*aa(c33,c22)
        aa(4  ,c22) = aa(4  ,c22) - aa(4  ,c33)*a3inv*aa(c33,c22)
        aa(1  ,4  ) = aa(1  ,4  ) - aa(1  ,c33)*a3inv*aa(c33,4  )
        aa(c22,4  ) = aa(c22,4  ) - aa(c22,c33)*a3inv*aa(c33,4  )
        aa(4  ,4  ) = aa(4  ,4  ) - aa(4  ,c33)*a3inv*aa(c33,4  )

        aa(1  ,c33) = 0.0d0
        aa(c22,c33) = 0.0d0
        aa(4  ,c33) = 0.0d0

        aa(c33,1  ) = 0.0d0
        aa(c33,c22) = 0.0d0
        aa(c33,4  ) = 0.0d0

        aa(c33,c33) = 0.0d0
      endif

      end
