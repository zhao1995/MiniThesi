c$Id:$
      subroutine me3fbmtx ( shx1 , shx2 , cpx1 , cpx2 , b )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FBMTX is the subroutine which computes the
c                      discrete strain-displacement operator (matrix)
c                      for the general membrane element in a cartesian
c                      reference frame.
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      shx1,shx2 ..... Nodal shape function cartesian derivatives.
c      cpx1,cpx2 ..... Current coordinate global derivatives
c                      at the Gauss points.

c      Routine Output:
c      ---------------
c      b ............. Discrete strain-displacement operator.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i , j   , nm
      real*8    shx1(4) , shx2(4) , cpx1(3) , cpx2(3) , b(3,*)

c     [Bmm] Part (Membrane)

      nm = 0
      do j = 1,4
        do i = 1,3
          b(1,nm+i) = shx1(j)*cpx1(i)
          b(2,nm+i) = shx2(j)*cpx2(i)
          b(3,nm+i) = shx2(j)*cpx1(i) + shx1(j)*cpx2(i)
        end do ! i
        nm = nm + 3
      end do ! j

      end

      subroutine me3fener ( d, dthk, ul, shp, fphi, xjw, ndf, lint)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute momenta and energy

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'erotas.h'
      include  'eldata.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'ptdat6.h'
      include  'tdata.h'

      integer   i      , j  , lint , ndf
      real*8    rhoa   , dthk , dot
      real*8    d(*)   , ul  (ndf,nen,*), shp(4,4)
      real*8    ce(3,4), fphi(3,4)      , sn(3,4), aa(3)
      real*8    pg(3)  , vd  (3)    , xjw (4)

      save

c     Compute inertial values

      rhoa = d( 4)*dthk

      do i = 1,lint

c       Interpolate Velocities and Accelerations

        do j = 1 , 3
          pg(j)   = shp(1,i)*fphi(j, 1)  + shp(2,i)*fphi(j, 2)
     &            + shp(3,i)*fphi(j, 3)  + shp(4,i)*fphi(j, 4)

          vd(j)   = shp(1,i)*ul  (j,1,4) + shp(2,i)*ul  (j,2,4)
     &            + shp(3,i)*ul  (j,3,4) + shp(4,i)*ul  (j,4,4)
        end do ! j

        call vecp ( pg , vd , aa )

c       Integrate Momenta

        do j = 1 , 3
          epl(j  ) = epl(j  ) + aa(j) * rhoa * xjw(i)
        end do ! j

c       Integrate Energy

        epl(7) = epl(7) + 0.5d0*dot(vd,vd,3) * rhoa  * xjw(i)
        epl(8) = epl(8) + 0.5d0*dot(sn(1,i),ce(1,i),3)

      end do ! i

      end

      subroutine me3flmda ( v1 , v2, xlm )

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FLMDA is a subroutine which computes a unique
c                      orthogonal array.
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      v1 ............. Tangent vector: zph1
c      v2 ............. Tangent vector: zph2

c      Routine Output:
c      ---------------
c      xlm ........... Orthogonal transformation matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i
      real*8    fac, v1(3), v2(3), xlm(3,2), w1(3), w2(3), e1(3), e2(3)

      fac = 1.d0/sqrt(v1(1)**2 +v1(2)**2 +v1(3)**2)
      do i = 1,3
        e1(i) = v1(i)*fac
      end do ! i

      fac = 1.d0/sqrt(v2(1)**2 +v2(2)**2 +v2(3)**2)
      do i = 1,3
        e2(i) = v2(i)*fac
      end do ! i

      do i = 1,3
        w1(i) = e1(i) - e2(i)
        w2(i) = e1(i) + e2(i)
      end do ! i

      fac = 1.d0/sqrt(w1(1)**2 +w1(2)**2 +w1(3)**2)
      do i = 1,3
        w1(i) = w1(i)*fac
      end do ! i

      fac = 1.d0/sqrt(w2(1)**2 +w2(2)**2 +w2(3)**2)
      do i = 1,3
        w2(i) = w2(i)*fac
      end do ! i

      do i = 1,3
        e1(i) = w1(i) + w2(i)
        e2(i) = w2(i) - w1(i)
      end do ! i

      fac = 1.d0/sqrt(e1(1)**2 +e1(2)**2 +e1(3)**2)
      do i = 1,3
        xlm(i,1) = e1(i)*fac
      end do ! i

      fac = 1.d0/sqrt(e2(1)**2 +e2(2)**2 +e2(3)**2)
      do i = 1,3
        xlm(i,2) = e2(i)*fac
      end do ! i

      end

      subroutine me3fgeom ( shx1 , shx2 , sn   , c , s )

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FGEOM is the subroutine which constructs
c                      the static geometric tangent stiffness.
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      shx1,shx2 ..... Nodal shape function cartesian derivatives.
c      sn ............ Membrane stress.

c      Routine Output:
c      ---------------
c      s ............. Element stiffness.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i  , j , k  , ni , nj
      real*8    termm  , sn1i    , sn2i
      real*8    shx1(4), shx2(4) , sn(3) , snl(3) , s(12,12) , c(3,3)

      do i = 1,3
        snl(i) = sn(i)
      end do ! i

      if(max(abs(snl(1)),abs(snl(1)),abs(snl(1))).le.0.0d0) then
        snl(1) = 1.0d-03*c(1,1)
        snl(2) = 1.0d-03*c(2,2)
      endif

      ni = 0
      do i = 1,4
        sn1i= snl(1)*shx1(i) + snl(3)*shx2(i)
        sn2i= snl(2)*shx2(i) + snl(3)*shx1(i)

        nj = 0
        do j = 1,4
          termm = sn1i*shx1(j) + sn2i*shx2(j)
          do k = 1,3
            s(ni+k,nj+k) = s(ni+k,nj+k) + termm
          end do ! k
          nj = nj + 3
        end do ! j
        ni = ni + 3
      end do ! i

      end

      subroutine me3fmasg( d, dthk, shp, xjw, s, lint, ndf, nst)

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FMASG is the subroutine which computes the
c                      mass matrix for rigid body and eigen analyses.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, j, k, l, j1, k1, lint, ndf, nst
      real*8    rhoa, dthk, fac, fac1, d(*), s(nst,*), shp(4,4), xjw(4)

c     Compute Mass Array:

      do i = 1,lint

c       Compute Weighted Density/Inertia

        rhoa = d( 4)*dthk * xjw(i)

c       Compute Mass
        j1 = 0
        do j = 1 , 4

          fac = rhoa*shp(j,i)

          k1 = 0
          do k = 1 , 4
            fac1 = fac*shp(k,i)
            do l = 1 , 3
              s(j1+l,k1+l) = s(j1+l,k1+l) + fac1
            end do ! l
            k1   = k1 + ndf
          end do ! k
          j1 = j1 + ndf
        end do ! j
      end do ! i

      end

      subroutine me3fshap ( xi , shp , shp1 , shp2 )

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FSHAP is the subroutine which computes the
c                      element shape functions and local derivatives
c                      for the 4-noded quadrilateral element.
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      xi(1),xi(2) ... Local coordinates of the current point.

c      Routine Output:
c      ---------------
c      shp ........... Nodal shape function values.
c      shp1,shp2 ..... Nodal shape function natural derivatives.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    xi(2) , shp(4) , shp1(4) , shp2(4)

c     Evaluate (xi) Derivative:

      shp1(2) =  0.25d0 * (1.d0-xi(2))
      shp1(3) =  0.25d0 * (1.d0+xi(2))
      shp1(4) = -shp1(3)
      shp1(1) = -shp1(2)

c     Evaluate (xi(2)) Derivative:

      shp2(3) =  0.25d0 * (1.d0+xi(1))
      shp2(4) =  0.25d0 * (1.d0-xi(1))
      shp2(1) = -shp2(4)
      shp2(2) = -shp2(3)

c     Evaluate Shape Function:

      shp (1) =  shp2(4) * (1.d0-xi(2))
      shp (2) =  shp2(3) * (1.d0-xi(2))
      shp (3) =  shp2(3) * (1.d0+xi(2))
      shp (4) =  shp2(4) * (1.d0+xi(2))

      end

      subroutine me3fmdsh ( shp1  , shp2 , xjinv , shx1  , shx2 )

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FMDSH is a subroutine which transforms local
c                      shape function gradients into global
c                      (Cartesian) gradients.
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      shp1,shp2 ..... Nodal shape function natural derivatives.
c      xjinv ......... Local-global coordinate Jacobian matrix.

c      Routine Output:
c      ---------------
c      shx1,shx2 ..... Nodal shape function cartesian derivatives.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   j
      real*8    shp1(4) , shp2(4) , shx1(4) , shx2(4), xjinv(2,2)

c     Transform Shape Function Derivatives:

      do j = 1 , 4
        shx1(j) = xjinv(1,1)*shp1(j) + xjinv(2,1)*shp2(j)
        shx2(j) = xjinv(1,2)*shp1(j) + xjinv(2,2)*shp2(j)
      end do ! j

      end

      subroutine me3fresd ( b      , sn     , p )

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FRESD is the routine which computes the
c                      static residual.
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      b ............. Discrete strain-displacement operator.
c      sn ............ Membrane stress.

c      Routine Output:
c      ---------------
c      p ............. Element residual.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    i
      real*8     b(3,12) , sn(3) , p(12)

c     Calculate Membrane

      do i = 1,12
        p(i) = p(i) - b(1,i)*sn(1) - b(2,i)*sn(2) - b(3,i)*sn(3)
      end do ! i

      end

      subroutine me3fstif ( b      , b1     , c      , s      )

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FSTIF is the routine which computes the
c                      material tangent, i.e., [B-t][C][B].
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      b ............. Discrete strain-displacement operator.
c      c ............. Resultant stress constitutive matrix.

c      Routine Output:
c      ---------------
c      s ............. Element stiffness.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i       , j
      real*8    b(3,12) , b1(3,12) ,c(3,3) , s(12,12) , temp(12,3)

c     Calculate Membrane [B-transpose] . [C]

      do j = 1,3
        do i = 1,12
          temp(i,j) = b(1,i)*c(1,j) + b(2,i)*c(2,j) + b(3,i)*c(3,j)
        end do ! i
      end do ! j

c     Calculate [B-transpose] . [C] . [B]

      do j = 1,12
        do i = 1,12
          s(i,j) = s(i,j) + temp(i,1)*b1(1,j) + temp(i,2)*b1(2,j)
     &                    + temp(i,3)*b1(3,j)
        end do ! i
      end do ! j

      end

      subroutine me3fplst ( sn,dt,st,shp,xjw,lint )

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:  ME3FPLST is the stress projection subroutine
c                    that computes the element contributions to the
c                    global equation system of the least-squares
c                    problem. The `mass' matrix is lumped with row-sum.
c                    The stresses are arranged as:
c                    1-3 := Membrane resultants.
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      sn ............ Membrane stress.
c      shp ........... Nodal shape function values.
c      xjw ........... Reference jacobian weighted for integration.
c      lint .......... Number quadrature points

c      Routine Output:
c      ---------------
c      dt ............ R.H.S. of stress projection equation system.
c      st ............ L.H.S. of stress projection equation system,
c                      lumped by row-sum.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'strnum.h'

      integer   i, l, lint
      real*8    sn(3,*), dt(*), st(nen,*), shp(4,*), xjw(*), xjshp

c     Loop over quadrature and nodes:

      do l = 1,lint
        do i = 1,4
          xjshp   = shp(i,l)*xjw(l)
          dt(i)   = dt(i)   + xjshp
          st(i,1) = st(i,1) + sn(1,l)*xjshp
          st(i,2) = st(i,2) + sn(2,l)*xjshp
          st(i,4) = st(i,4) + sn(3,l)*xjshp
        end do ! i
      end do ! l

      iste = 4

      end

      subroutine me3fstre ( xl, ce, sn, shp, xjw, lint, ndm, isw )

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Stress output routine

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'iofile.h'
      include  'prstrs.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   i, j, lint, ndm, isw
      real*8    xl(ndm,*),shp(4,4),zphg(3,4),ce(3,4),sn(3,4),xjw(4)

      save

c     Output Element Stresses:

      if(isw.eq.4) then
        mct = mct - 1
        if (mct.le.0) then
          if (ior.lt.0) write (*,4000)
          write (iow,4000)
          mct = 100000
        endif

c       Gauss Loop

        do i = 1 , lint

c         Interpolate Position

          do j = 1 , 3
            zphg(j,i) = shp(1,i)*xl(j,1) + shp(2,i)*xl(j,2)
     &                + shp(3,i)*xl(j,3) + shp(4,i)*xl(j,4)
          end do ! j

c         Output to Screen

          if (ior.lt.0) then
            write (*,4100)  n, i,
     &                      (zphg(j,i),        j=1,3) ,
     &                      (ce  (j,i),        j=1,3) ,
     &                      (sn  (j,i)/xjw(i), j=1,3)
          endif

c         Output to File

          write (iow,4100)  n, i,
     &                      (zphg(j,i),        j=1,3) ,
     &                      (ce  (j,i),        j=1,3) ,
     &                      (sn  (j,i)/xjw(i), j=1,3)
        end do ! i

c     Compute Projected Nodal Stresses:

      elseif(isw.eq.8) then

       call me3fplst ( sn, hr(np(35)), hr(np(36)),shp, xjw, lint )

      endif

c     Element Stress Format:

4000  format(/
     & ' ----------------------------------------',
     & '------------------------------------- '/,
     & '  F i n i t e   D e f o r m a t i o n ',
     & '  M e m b r a n e   E l e m e n t     '/,
     & ' ----------------------------------------',
     & '------------------------------------- '//,
     & '  Element      Gauss Pt # '/,
     & '  Midsrfce-X   Midsrfce-Y   Midsrfce-Z '/,
     & '  MemStrn-xx   MemStrn-yy   MemStrn-xy ',
     & '  MemStrs-xx   MemStrs-yy   MemStrs-xy  '/,
     & ' ------------ ------------ ------------',
     & ' ------------ ------------ ------------ '/)

4100  format(i7,12x,i1/,3e13.6/,6e13.6/)

      end

      subroutine me3fstrn ( cpx1 , cpx2 , ce )

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FSTRN is the subroutine which computes the
c                      strain measures for the membrane element in a
c                      cartesian reference frame.
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      cpx1,cpx2 ..... Current coordinate global derivatives
c                      at the Gauss points.

c      Routine Output:
c      ---------------
c      ce............. Current membrane strain measure.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    cpx1 (3) , cpx2 (3) , ce   (3)

c     Compute Membrane Strains:

      ce(1) = (cpx1(1)*cpx1(1)+cpx1(2)*cpx1(2)+cpx1(3)*cpx1(3))*0.5d0
      ce(2) = (cpx2(1)*cpx2(1)+cpx2(2)*cpx2(2)+cpx2(3)*cpx2(3))*0.5d0
      ce(3) =  cpx1(1)*cpx2(1)+cpx1(2)*cpx2(2)+cpx1(3)*cpx2(3)

      end

      subroutine me3fcmtx ( d , ce, dthk , xjw  , c   )

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FCMTX is the subroutine which sets-up the
c                      linear resultant constitutive relations.
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      d(*) .......... Material parameters from inpt2d
c      ce(3).......... Strains at point
c      dthk .......... Element thickness
c      xjw ........... Jacobian times quadrature weight

c      Routine Output:
c      ---------------
c      c ............. Resultant stress constitutive matrix.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i   , j
      real*8    dthk, xjw, fac, facj, d(*), ce(3), c(3,3), sl(3)

c     Fung Pseudo-Exponential Model

      if(nint(d(20)).eq.7) then
        sl(1)  = d(21)*ce(1) + d(23)*ce(2)
        sl(2)  = d(23)*ce(1) + d(22)*ce(2)
        sl(3)  = d(24)*ce(3)
        fac    = 2.d0*d(30)*exp(ce(1)*sl(1)+ce(2)*sl(2)+ce(3)*sl(3))
        fac    = fac*dthk*xjw
        c(1,1) = fac*d(21)
        c(1,2) = fac*d(23)
        c(1,3) = 0.d0
        c(2,2) = fac*d(22)
        c(2,3) = 0.d0
        c(3,3) = fac*d(24)
        fac    = 2.d0*fac
        do j = 1,3
          facj = fac*sl(j)
          do i = 1,3
            c(i,j) = c(i,j) + sl(i)*facj
          end do ! i
        end do ! j

c     St Venant-Kirchhoff Model

      else
        c(1,1) = d(21)*dthk*xjw
        c(1,2) = d(24)*dthk*xjw
        c(1,3) = 0.d0
        c(2,2) = d(22)*dthk*xjw
        c(2,3) = 0.d0
        c(3,3) = d(27)*dthk*xjw
      endif
      c(2,1) = c(1,2)
      c(3,1) = c(1,3)
      c(3,2) = c(2,3)
      end

      subroutine me3fstrs ( d , c , ce , dthk, xjw, sn )

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FSTRS is a subroutine which computes the
c                      stress resultants, given the strain measures
c                      and a linear constitutive matrix.

c      Routine Input:
c      --------------
c      d(*)........... Material parameters from inpt2d.
c      c ............. Resultant stress constitutive matrix.
c      ce............. Current membrane strain measure.
c      dthk .......... Element thickness
c      xjw ........... Jacobian times quadrature weight

c      Routine Output:
c      ---------------
c      sn ............ Membrane stress.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    fac, dthk, xjw, d(*), c(3,3), ce(3), sn(3)

c     Membrane Stresses:

      if(nint(d(20)).eq.7) then ! Fung model
        sn(1) = d(21)*ce(1) + d(23)*ce(2)
        sn(2) = d(23)*ce(1) + d(22)*ce(2)
        sn(3) = d(24)*ce(3)
        fac   = 2.d0*d(30)*exp(ce(1)*sn(1)+ce(2)*sn(2)+ce(3)*sn(3))
        fac   = fac*dthk*xjw
        sn(1) = fac*sn(1)
        sn(2) = fac*sn(2)
        sn(3) = fac*sn(3)
      else                 ! Elastic model (ce(*) multiplied by jac)
        sn(1) = c(1,1)*ce(1) + c(1,2)*ce(2) + c(1,3)*ce(3)
        sn(2) = c(2,1)*ce(1) + c(2,2)*ce(2) + c(2,3)*ce(3)
        sn(3) = c(3,1)*ce(1) + c(3,2)*ce(2) + c(3,3)*ce(3)
      endif

      end

      subroutine me3fsurf(d,xl,ul, xjw, ndf,ndm,nst, p,s)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Nodal force and tangent array for pressure loading

c      INPUT variables
c        d(10)      Value of constant pressure on face
c        xl(ndm,4)  Nodal coordinates
c        ul(ndf,4)  Nodal displacements
c        xjw(*)     Surface jacobian time quadrature weight
c        ndf        Number of DOF / node
c        ndm        Space dimension
c        nst        Dimension of residual vector

c      OUTPUT variables
c        p(ndf,4)   Contribution to residual
c        s(nst,nst) Contribution to striffness matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'eltran.h'

      logical   bflg
      integer   l, lint, ndf, ndm, nst, ii, jj, i1, j1, i, j
      real*8    pn, pp, xii(4), eti(4)
      real*8    sg(3,4), shp(3,4), xu(3,4), dx(3,2), bg(3),bl(3)
      real*8    d(*), xl(ndm,*), ul(ndf,*), p(ndf,*), s(nst,*), xjw(*)
      real*8    bf(3),bt(3,3), xx(3)

      save

      data      xii / -0.5d0, 0.5d0, 0.5d0,-0.5d0 /
      data      eti / -0.5d0,-0.5d0, 0.5d0, 0.5d0 /
      data      bl  /  3*0.0d0 /

c     Set body loading factors

      call sbodyf(d, bg)

      if(int(d(76)).gt.0) then
        bl(3) = d(10)
      else
        bl(3) = d(10)*dm
      endif

c     Compute nodal coordinates in correct reference frame

      if(d(68).eq.0.0d0) then
        do ii = 1,4
          xu(1,ii) = xl(1,ii)
          xu(2,ii) = xl(2,ii)
          xu(3,ii) = xl(3,ii)
        end do ! ii
      else
        do ii = 1,4
          xu(1,ii) = xl(1,ii) + ul(1,ii)
          xu(2,ii) = xl(2,ii) + ul(2,ii)
          xu(3,ii) = xl(3,ii) + ul(3,ii)
        end do ! ii
      endif

c     Get quadrature information

      l  = 2
      call int2d (l, lint, sg)

c     First loop over quadrature points

      bflg = d(4).gt.0.0d0 .and. d(65).gt.0.0d0  ! angular velocity

      do l = 1,lint

c       Compute shape functions and geometric factors

        do ii = 1,4
          pp        = 0.5d0+xii(ii)*sg(1,l)
          pn        = 0.5d0+eti(ii)*sg(2,l)
          shp(1,ii) = xii(ii)*pn
          shp(2,ii) = eti(ii)*pp
          shp(3,ii) = pn*pp
        end do ! ii

        do ii = 1,3
          dx(ii,1) = shp(1,1)*xu(ii,1) + shp(1,2)*xu(ii,2)
     &             + shp(1,3)*xu(ii,3) + shp(1,4)*xu(ii,4)
          dx(ii,2) = shp(2,1)*xu(ii,1) + shp(2,2)*xu(ii,2)
     &             + shp(2,3)*xu(ii,3) + shp(2,4)*xu(ii,4)
        end do ! ii

c       Angular velocity: d(4) = rho; d(65) = omega

        do i = 1,3
          bf(i) = 0.0d0
        end do ! i
        if(bflg) then
          do ii = 1,3
            xx(ii) = 0.0d0
            do jj = 1,4
              xx(ii) = xx(ii) + shp(3,jj)*xu(ii,jj)
            end do ! jj
          end do ! ii
          call sbodyw(d(4),d(65),xx, bf,bt, .true.)
          do ii = 1,3
            do jj = 1,3
              bt(jj,ii) = bt(jj,ii)*xjw(l)
            end do ! jj
          end do ! ii
        endif

c       Compute nodal loads for pressures

        pn = bl(3)*sg(3,l)
        do ii = 1,4
          pp     = shp(3,ii)*pn
          p(1,ii) = p(1,ii) + pp*(dx(2,1)*dx(3,2) - dx(3,1)*dx(2,2))
     &                      + (bf(1)+bg(1))*shp(3,ii)*xjw(l)
          p(2,ii) = p(2,ii) + pp*(dx(3,1)*dx(1,2) - dx(1,1)*dx(3,2))
     &                      + (bf(2)+bg(2))*shp(3,ii)*xjw(l)
          p(3,ii) = p(3,ii) + pp*(dx(1,1)*dx(2,2) - dx(2,1)*dx(1,2))
     &                      + (bf(3)+bg(3))*shp(3,ii)*xjw(l)
        end do ! ii

c       Compute follower surface load  tangent if necessary

        if(d(68).gt.0.0d0) then
          i1 = 0
          do ii = 1,4
            pp = shp(3,ii)*pn*ctan(1)
            j1 = 0
            do jj = 1,4
              s(i1+1,j1+2) = s(i1+1,j1+2)
     &                     - pp*(shp(1,jj)*dx(3,2) - dx(3,1)*shp(2,jj))
              s(i1+2,j1+3) = s(i1+2,j1+3)
     &                     - pp*(shp(1,jj)*dx(1,2) - dx(1,1)*shp(2,jj))
              s(i1+3,j1+1) = s(i1+3,j1+1)
     &                     - pp*(shp(1,jj)*dx(2,2) - dx(2,1)*shp(2,jj))
              s(i1+1,j1+3) = s(i1+1,j1+3)
     &                     - pp*(shp(2,jj)*dx(2,1) - dx(2,2)*shp(1,jj))
              s(i1+2,j1+1) = s(i1+2,j1+1)
     &                     - pp*(shp(2,jj)*dx(3,1) - dx(3,2)*shp(1,jj))
              s(i1+3,j1+2) = s(i1+3,j1+2)
     &                     - pp*(shp(2,jj)*dx(1,1) - dx(1,2)*shp(1,jj))
              j1 = j1 + ndf
            end do ! jj
            i1 = i1 + ndf
          end do ! ii
        endif

c       Compute rotational body force tangent if necessary

        if(bflg) then
          i1 = 0
          do ii = 1,4
            pp = shp(3,ii)*ctan(1)
            j1 = 0
            do jj = 1,4
              do i = 1,3
                do j = 1,3
                  s(i1+1,j1+2) = s(i1+1,j1+2) + pp*bt(i,j)*shp(3,jj)
                end do ! j
              end do ! i
              j1 = j1 + ndf
            end do ! jj
            i1 = i1 + ndf
          end do ! ii
        endif

      end do ! l

      end

      subroutine me3ftran ( dt   , shp  , aa   , rhoa , p    , s  )

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    ME3FTRAN is the subroutine which computes the
c                      transient (inertial) terms of the residual
c                      vector and tangent stiffness matrix
c                      for the Energy-Momentum method

c      Caution:        Must use the conservative time-stepping
c                      algorithms in conjunction with this routine,
c                      i.e., beta,cons, etc... (nop=5).

c      Routine Input:
c      --------------
c      dt ............ Time step.
c      shp ........... Nodal shape function values.
c      aa ............ Gauss point acceleration at T-n+a.

c      rhoa .......... Trans. thickness weighted density multiplyied
c                      by the jacobian and weight at the Gauss point.

c      Routine Output:
c      ---------------
c      s ............. Element stiffness.
c      p ............. Element residual.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   j, k, l, j1, k1
      real*8    dt, dtr, dtr4, rhoa, fac, fac1
      real*8    shp(4), aa(3), s(12,12), p(3,4)

c     Compute transient terms

      dtr  = 1.d0/dt
      dtr4 = 4.d0*dtr*dtr

      do j = 1 , 4

c       Set Up Residual

        do k = 1 , 3
          p(k,j) = p(k,j) - shp(j) * rhoa * aa(k)
        end do ! k

      end do ! j

      j1 = 0
      do j = 1 , 4

        fac = rhoa*shp(j) * dtr4

        k1 = 0
        do k = 1 , 4
          fac1 = fac*shp(k)
          do l = 1 , 3
            s(j1+l,k1+l) = s(j1+l,k1+l) + fac1
          end do ! l
          k1 = k1 + 3
        end do ! k

        j1 = j1 + 3
      end do ! j

      end
