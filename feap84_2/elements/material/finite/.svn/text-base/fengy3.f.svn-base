c$Id:$
      subroutine fengy3(d,detg,u,up,upp, ha,hp,hpp, jsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add Model 4.                                     20/12/2006
c       2. Correct Model 4 derivatives for wrong factors    20/11/2011
c          of 2.0d0 (in both up and upp). Add small strain
c       3. Remove factor of 2 from definition of U(J) in 4  24/11/2011
c       4. Return up = upp = 0 for jsw = 5                  16/02/2012
c       5. Modify series expansion form for Model 2 'up'    07/06/2013
c       6. Multiply series for 'up' in Model 1 by 0.5       22/06/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute pressure-entropy function

c     Inputs:
c       d(*)        - Material parameter array
c       detg        - Deformation gradient - 1
c       detf        - Deformation gradient
c       jsw         - Function type: 1 - K*[ 0.25*( J^2 - 1 )-0.5*ln(J)]
c                                    2 - K*[ 0.5*( J - 1 )^2 ]
c                                    3 - K*[ 0.5*( ln(J) )^2 ]
c                                    4 - K*[ (J - 1 - ln(J)) ]
c                                    5 - Incompressible

c     Outputs:
c       u           - Internal energy
c       up          - First derivative of internal energy
c       upp         - Second derivative of internal energy
c       ha          - Augmentation function
c       hp          - First derivative of augmentation function
c       hpp         - Second derivative of augmentation function
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'

      integer   jsw
      real*8    d1, dh, detg, detf, u, up, upp, ha, hp, hpp, d(*)

c     Perform (augmented) iteration on penalty

      d1    = augf*d(21)
      detf  = detg + 1.0d0

c     Free energy function
c     Psi   = U(J) -3*alpha*J*U'(J)*[T - T_ref]

c     up    = partial_J ( Psi )
c     upp   = [ partial_J ( J partial_J Psi ) - up ]/J
c           =   partial^2_J ( Psi )/J

c     Current volumetric functions are:

c     Model 1.) U(J) = K*0.25*(J^2 - 1 - 2*(log J))

      if    (jsw.eq.1) then

        dh   = d1*0.5d0
        u    = dh*( 0.5d0*detf**2 - 0.5d0  - log(abs(detf)))
        if(abs(detg).lt.0.001d0) then
          up = dh*( detg+(1.d0-detg)*(detg+detg**3+detg**5+detg**7))
        else
          up = dh*( detf - 1.d0/detf    )
        endif
        upp  = dh*( 1.d0 + 1.d0/detf**2 )

c     Model 2.) U(J) = K*0.5*(J-1)^2

      elseif(jsw.eq.2) then

        u    = d1*detg**2*0.5d0
        up   = d1*detg
        upp  = d1

c     Model 3.) U(J) = K*0.5*(log J)^2

      elseif(jsw.eq.3) then

c       up   = ( d1*log( detf ) - 3*K*alpha*(T - Tref) )/detf

        u    = d1*log(abs(detf))**2*0.5d0
        if(abs(detg).lt.0.001d0) then
          up = d1*(detg - 1.5d0*detg**2 + (110.d0*detg**3
     &       - 125.d0*detg**4 + 137.d0*detg**5 - 247.d0*detg**6)/60.d0)
        else
          up = d1*log(abs(detf)) / detf
        endif
        upp  = ( d1/detf  - up )/detf

c     Model 4.) U(J) = K*(J - 1 - log J)

      elseif(jsw.eq.4) then

        u    = d1*(detg - log(abs(detf)))
        if(abs(detg).lt.0.001d0) then
          up = d1*(detg - detg**2 + detg**3 - detg**4
     &                  + detg**5 - detg**6)
        else
          up = d1 - d1/detf
        endif
        upp  = d1/(detf*detf)

c     Model 5.) Incompressible

      else

        up  = 0.0d0
        upp = 0.0d0

      endif

c     Augmented Lagrangian function and derivatives
c     Current augmented Lagrangian function is

c     Model 1.) h(J) = (J - 1)

      ha  = detg
      hp  = 1.0d0
      hpp = 0.0d0

      end
