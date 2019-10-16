c$Id:$
      subroutine nalp3f(d,f,finv,detf,b, hn,ntm, tau,aa,estore, isw,jsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Cast d(40) as integer for compare                07/06/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Control routine for constitutive models formulated in
c               principal stretches.

c      Inputs:
c         d(*)       Material parameters
c         f(3,3)     Deformation gradient
c         finv(3,3)  Inverse deformation gradient
c         detf       Determinant of deformation gradient
c         b(6)       Elastic Left Cauchy-Green tensor
c         hn(*)      History variables at t_n
c         ntm        Number of active stress components
c         isw        Shear Model type: 1 = Ogden
c                                      2 = Logarithmic stretch
c         jsw        Volumetric Model: 1,2,3 (see fengy3)

c      Outputs:
c         hn(*)      History variables at t_n+1
c         tau(6)     Kirchhoff stress tensor  (deviatoric)
c         aa(6,6)    Kirchhoff tangent moduli (deviatoric)
c         estore     Stored energy density
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eltran.h'
      include  'cdamag.h'

      integer   i,isw, j,jsw, k, l, nrot, ntm
      real*8    d(*), f(3,3), finv(3,3), b(6), hn(*), tau(6), aa(6,6)
      real*8    qen(3,3),q(6,6),epsd(3),taup(3),aap(6,6),detf,u,estore

      save

c     Compute principal stretches and directions

c     Initialization: Upon entry to EIG3 array qen(3,3) contains left
c     Cauchy-Green tensor b; upon exit, qen(3,3) contains eigenvectors.

      do i = 1,3
       j = 1 + mod(i,3)
       qen(i,i) = b(i)
       qen(i,j) = b(i+3)
       qen(j,i) = b(i+3)
      end do ! i

      call eig3(qen,bpr,nrot)

c     Construct matrix to change basis

      do i=1,3
        k=1+mod(i,3)
        do j=1,3
          l=1+mod(j,3)

          q(i,j)     = qen(i,j)*qen(i,j)

          q(i+3,j)   = qen(i,j)*qen(k,j)

          q(j,i+3)   = qen(j,i)*qen(j,k)*2.d0

          q(j+3,i+3) = qen(j,i)*qen(l,k) + qen(j,k)*qen(l,i)

        end do ! j
      end do ! i

c     Compute deviatoric stresses and matl tangent in principal basis

      if    (isw.eq.1) then
        call wder3f(d,bpr, taup,aap)
      elseif(isw.eq.2) then
        call wlog3f(d,bpr, epsd, taup,aap,wengy)
      endif

c     Compute deviatoric stresses and matl tangent in standard basis

      call amat3f(taup,aap,q, tau,aa)

c     Compute inelastic parts

c     1. Plasticity (logrithmic principal stretch only)

      if(nint(d(40)).eq.1 .and.isw.eq.2) then

c     2. Viscoelasticity

      elseif(nint(d(40)).eq.2) then

        estore = wengy    ! Used for damage calculations
        call viscfd(d,f,finv,detf,hn(1),hn(2),hn(2+ntm),ntm, tau,aa)

      endif

c     Add volumetric part to stress and tangent matrix

      call apre3f(d,detf, tau,aa,u, jsw)

c     Compute stored energy

      estore = u + wengy

      end
