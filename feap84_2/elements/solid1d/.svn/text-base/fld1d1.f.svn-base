c$Id:$
      subroutine fld1d1(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set quadrature for d(5).eq.0                     26/03/2009
c       2. Use 'quard1d' and 'interp1d' for solution        01/04/2009
c       3. Compute xref and xcur for possible  modlfd use   25/07/2009
c       4. Add flag to call list                            03/03/2010
c       5. Add 'l' to modlfd call                           05/01/2012
c       6. Add epsl on call to slcn1d, almansi strain       01/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  1-D Finite Deformation Elasticity Routine
c                Remark: This is a standard displacement model

c      Inputs:
c         d(*)      - Material set parameters
c         ul(ndf,*) - Nodal solution parameters for element
c         xl(ndm,*) - Nodal coordinates for element
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Dimension of element arrays
c         isw       - Switch to control action

c      Outputs:
c         s(nst,*)  - Element matrix
c         p(nst)    - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'elcoor.h'
      include  'eldata.h'
      include  'elengy.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'ptdat6.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'comblk.h'

      integer   ndf,ndm,nst,isw,i,i1,j,jj,j1,l,nn,nhv,istrt
      real*8    bd3,dl,dc,di,dmas0
      real*8    cfac,lfac,xx1,yy1, tempi, xlamd, ha, ta
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),body(3)
      real*8    r(ndf,*),r1(8),al,ac,vl,weng(8)
      real*8    df(3,3,8),f(9,2,8),finv(3,3,8),detf(2,8)
      real*8    ds(6,6,5),dd(6,6),sigv(9),sigl(16,8),epsv(3)
      real*8    bbd(3),bb(6),shr(8),dvol(8),epsl(3,8)
      real*8    egreen(6),ealmansi(6)

      save

      data      xlamd,ha /  2*0.0d0 /
      data      f    /144*0.0d0 /
      data      finv / 72*0.0d0/
      data      ta   /    0.0d0/

c     Set quadrature

      call quadr1d(d)

      xref(2) = 0.0d0
      xref(3) = 0.0d0
      xcur(2) = 0.0d0
      xcur(3) = 0.0d0

c     COMPUTE TANGENT STIFFNESS AND RESIDUAL FORCE VECTOR

c     Compute shape functions and derivatives in reference configuration

      do l = 1,lint

c       Set shape functions for quadrature point

        call interp1d(l, xl, ndm,nel,.false.)

        if(stype.eq.3) then
          xx1 = 0.0d0
          do i = 1,nel
            xx1 = xx1 + shp1(2,i,l)*xl(1,i)
          end do ! i
          jac(l) = jac(l)*xx1
        endif
      end do ! l

c     Compute deformation gradient and determinant; transform shape
c     functions to current configuration.

      do l = 1,lint
        call kine1d(shp1(1,1,l),xl,ul,f(1,1,l),finv(1,1,l),df(1,1,l),
     &              detf(1,l),ndm,ndf,nel,nen)
      end do ! l

c     Integration order set to static

      if((d(7).ge.0.0 .or. d(183).ne.0.0d0) .and.
     &           (ndfo(1).gt.0 .or. shflg)) then
        cfac = d(7)
        lfac = 1.0d0 - cfac
      else
        cfac = 0.0d0
        lfac = 0.0d0
      endif

      nhv   = nint(d(15))
      istrt = nint(d(84))

      if(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16 .or. isw.eq.25 ) go to 4

      call sbodyf(d, body)

c     LOOP OVER GAUSS POINTS

      nn = 0
      do l = 1,lint

c       Set reference and current coordinates

        xref(1) = 0.0d0
        xcur(1) = 0.0d0
        do i=1,nel
          xref(1) = xref(1) + shp1(2,i,l)*xl(1,i)
          xcur(1) = xcur(1) + shp1(2,i,l)*ul(1,i,1)
        end do ! i
        xcur(1) = xcur(1) + xref(1)

c       Check for axisymmetry

        if(stype.eq.3) then
          do i = 1,nel
            shr(i) = shp1(2,i,l)/xref(1)
          end do ! i
        else
          do i = 1,nel
            shr(i) = 0.0d0
          end do ! i
        end if
        dvol(l) = jac(l)*detf(1,l)
        dmas0   = jac(l)*d(4)

c       Compute Cauchy stresses and spatial tangent tensor

        call modlfd(l,d,f(1,1,l),finv(1,1,l),df(1,1,l),detf(1,l),ta,
     &             hr(nn+nh1),hr(nn+nh2),nhv,istrt, ds,sigv,bb,
     &             xlamd,ha,.false.,isw)

        if(isw.eq.13) then

          epl(8) = epl(8) + estore*jac(l)

c         Compute velocity at point

          vl = 0.0d0
          do i = 1,nel
            vl = vl + shp1(2,i,l)*ul(1,i,4)
          end do ! i

          tempi = 0.0d0
          do i = 1,nel
            tempi = tempi + ul(1,i,4)**2*shp1(2,i,l)
          end do ! i

c         Accumulate kinetic energy

          epl(7) = epl(7) + 0.5d0*(lfac*tempi + cfac*vl**2)*dmas0

        elseif(isw.ne.14) then

c         Store stress values for tplot

          j1 = 6*(l-1)
          do j = 1,4
            tt(j+j1) = sigv(j)
          end do ! j

c         Multiply tangent moduli and stresses by volume element.

          do i = 1,3
            sigv(i) = sigv(i)*dvol(l)
            do j = 1,3
              dd(i,j) = ds(i,j,1)*dvol(l)*ctan(1)
            end do ! j
          end do ! i

c         Compute accelerations

          al = 0.0d0
          do i = 1,nel
            al = al + shp1(2,i,l)*ul(1,i,5)
          end do ! i
          al = al*cfac

c         COMPUTE STRESS DIVERGENCE AND INERTIA TERMS

          xx1 = xref(1)*d(65)**2

          do i = 1,nel

c           Compute inertial and body load effects

            ac = (al + lfac*ul(1,i,5))*dmas0
     &         -  body(1)*jac(l) - xx1*dmas0

c           Stress divergence term (used in geometric stiffness)

            r1(i) = shp1(1,i,l)*sigv(1)

c           Element residual

            r(1,i) = r(1,i) - r1(i) - ac*shp1(2,i,l)
     &                      - shr(i)*sigv(3)

          end do ! i

c         COMPUTE K (s(nst,nst) = K)

          if(isw.eq.3) then

c           PART 1. - Geometric and inertial part.

            dc  = cfac*ctan(3)*dmas0
            dl  = lfac*ctan(3)*dmas0
            i1  = 1
            do i = 1,nel

              s(i1,i1) = s(i1,i1) + shp1(2,i,l)*dl

              di  = dc*shp1(2,i,l)

c             Include geometric stiffness

              if(gflag) then
                bd3 = shr(i)*sigv(3)*ctan(1)
                j1  = 1
                do j = 1,nel
                  s(i1,j1) = s(i1,j1) + r1(i)*shp1(1,j,l)*ctan(1)
     &                     + di*shp1(2,j,l) + bd3*shr(j)
                  j1 = j1 + ndf
                end do ! j

c             Include inertia only

              else
                j1  = 1
                do j = 1,nel
                  s(i1,j1) = s(i1,j1) + di*shp1(2,j,l)
                  j1       = j1 + ndf
                end do ! j
              endif

              i1 = i1 + ndf
            end do ! i

c           PART 2. - Tangent modulus part (based upon dd-array)

            i1 = 1
            do i  = 1,nel

c             Compute bmat-t * dd * dvol

              do jj = 1,3
                bbd(jj) = shp1(1,i,l)*dd(1,jj)
     &                  + shr(   i  )*dd(3,jj)
              end do ! jj

c             Compute tangent stiffness

              j1 = 1
              do j  = 1,nel
                s(i1,j1) = s(i1,j1) + bbd(1)*shp1(1,j,l)
     &                              + bbd(3)*shr(   j  )
                j1 = j1 + ndf
              end do ! j

              i1 = i1 + ndf
            end  do ! i

          endif ! end of tangent

        endif ! end of isw options

        nn = nn + nhv

      end do ! l

      return

c     OUTPUT STRESSES

   4  do i = 1,3
        sigv(i) = 0.0d0
        epsv(i) = 0.0d0
      end do ! i

c     LOOP OVER GAUSS POINTS

      nn  = 0
      xx1 = 0.d0
      yy1 = 1.d0/dble(nel)

      do l = 1,lint

        xref(1) = 0.0d0
        xcur(1) = 0.0d0
        do i=1,nel
          xref(1) = xref(1) + shp1(2,i,l)*xl(1,i)
          xcur(1) = xcur(1) + shp1(2,i,l)*ul(1,i,1)
        end do ! i
        xcur(1) = xcur(1) + xref(1)
        xx1 = xx1 + xref(1)/yy1

c       Compute Cauchy stresses and spatial tangent tensor at t-n+1

        call modlfd(l,d,f(1,1,l),finv(1,1,l),df(1,1,l),detf(1,l),ta,
     &             hr(nn+nh1),hr(nn+nh2),nhv,istrt,ds,sigl(1,l),bb,
     &             xlamd,ha,.false.,isw)
        weng(l) = estore
        dvol(l) = jac(l)*detf(1,l)

c       Compute Green-Lagrange strains and Almansi strains

        call fstrain(f(1,1,l),finv(1,1,l), egreen, ealmansi)

        epsl(1:3,l) = ealmansi(1:3)

c       Compute average stresses and jacobian for printing

        do i = 1,3
          sigv(i) = sigv(i) + yy1 * sigl(i,l)
          epsv(i) = epsv(i) + yy1 * epsl(i,l)
        end do ! i

        nn = nn + nhv

      end do ! l

c     Output stresses

      if(isw.eq.4) then
        mct = mct - 2
        if(mct.le.0) then
          write(iow,2001) o,head
          if(ior.lt.0) write(*,2001) o,head
          mct = 50
        endif

        write(iow,2002) n,ma,xx1,(sigv(jj),jj=1,3),epsv
        if(ior.lt.0) then
          write(*,2002) n,ma,xx1,(sigv(jj),jj=1,3),epsv
        end if

c     Project stress values to nodes

      elseif(isw.eq.8) then

        call slcn1d(sigl,epsl,shp1,jac,r,s,r(nen+1,1),lint,nel,16)

c     Compute fracture indices

      elseif(isw.eq.16) then

        call pfrac1f(f,detf,sigl,weng, shp1,dvol, r, lint,ndf,ndm,3)

c     Compute Z-Z projections

      elseif(isw.eq.25) then

        call stcn1z(xl,sigl,shp1,dvol,lint,ndm,nel,16)

      end if

c     Format statements

2001  format(a1,20a4//5x,'Element Stresses'//'    Elmt Mat',
     1   '     1-coord   11-stress   22-stress   33-stress'/12x,
     2   '     Almansi   11-strain   22-strain   33-strain')

2002  format(i8,i4,1p,4e12.3/24x,1p,3e12.3)

      end
