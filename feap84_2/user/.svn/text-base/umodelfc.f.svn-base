c$Id:$
      subroutine umodelfc(umat,f,finv,df,detf,be,gradt,ta,d,ud,hn,h1,
     &                   nh,ntm,istrt, sig,flux,dd,kt, xlamd,ha, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: User Constitutive Model for coupled problems

c     Input:
c          umat      - User material type
c          f(3,3)    - Deformation gradient
c          finv(3,3) - Inverse deformation gradient
c          df(3,3)   - Incremental deformation gradient
c          detf      - Determinant of deformation gradient
c          be(6)     - Left Cauchy-Green tensor
c          ta        - Temperature change
c          d(*)      - Program material parameters (ndd)
c          ud(*)     - User material parameters (nud)
c          hn(nh)    - History terms at point: t_n
c          h1(nh)    - History terms at point: t_n+1
c          nh        - Number of history terms
c          ntm       - Number of stress components
c          istrt     - Start state: 0 = elastic; 1 = last solution
c          isw       - Solution option from element

c     Output:
c          sig(6)    - Stresses at point.
c          dd(6,6)   - Current material tangent moduli
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'elengy.h'

      integer    umat,nh,ntm,istrt,isw
      real*8     ta, xlamd,ha
      real*8     f(*),finv(*),df(*),detf(*),be(*),d(*),ud(*)
      real*8     hn(nh),h1(nh), sig(6),dd(7,7), flux(4),kt(3,3)
      real*8     gradt(3)

      save

c     Material Model 1: ?

      if(umat.eq.1) then

c       Insert model computations here: returns 'sig' and 'dd'
c                                       (Cauchy stress values)
c                                          and  'flux' and 'kt'
      endif

      end
