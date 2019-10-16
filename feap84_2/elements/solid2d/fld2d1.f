c$Id:$
      subroutine fld2d1(d,ul,xl,ix,tl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Return on augmented call                         26/01/2007
c       2. Shift storage for history variables +2           02/02/2007
c       3. Increase quadrature order for 3-node to l=3      21/03/2007
c       4. Add array 'bdy' to control use of body loading   20/07/2007
c          Add radial body loading
c       5. Add output of principal stresses                 04/03/2008
c       6. Add direct call to quadr2d and interp2d          11/11/2008
c       7. Remove 'nel' from call to 'quadr2d'              23/01/2009
c       8. Add accumualtion of average thickness stress     24/02/2009
c       9. Add 'defgrd.h' for deformation gradient values   27/04/2009
c      10. Dimension arrays for 64 nodes & quadrature pts.  04/05/2009
c      11. Add prints Almansi strains                       19/10/2009
c      12. Add temperature to argument list                 08/01/2010
c      13. Shift index on saving f(3,3) for stype.eq.1      26/03/2010
c      14. Reorder print of stress and strain on 2001/2002  18/05/2010
c      15. Move strain computation after call to modlfd for 17/09/2011
c          isw.eq.4 and isw.eq.8
c      16. Add 'l' to modlfd call                           05/01/2012
c      17. Modify to store displacement gradieint for       10/01/2012
c          plane stress problems.
c      18. Add average of density for multiscale            10/05/2012
c      19. Use 'jvol' for arrays 'jac' for projection       03/10/2012
c      20. Add eps on call to slcn2d                        01/01/2013
c      21. Remove use of kineps -- done by model now        29/06/2013
c          Also save of deformation gradient after modlfd
c      22. Pass strains to stcn2z for z-zhu projections     01/01/2014
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  2-D Finite Deformation Elasticity Routine
c                Remark: This is a standard displacement model

c      Inputs:
c         d(*)      - Material set parameters
c         ul(ndf,*) - Nodal solution parameters for element
c         xl(ndm,*) - Nodal coordinates for element
c         ix(*)     - Element nodal connection list
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Dimension of element arrays
c         isw       - Switch to control action

c      Outputs:
c         s(nst,*)  - Element matrix
c         p(nst)    - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'debugs.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'defgrd.h'
      include  'elcoor.h'
      include  'eldata.h'
      include  'elengy.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'oelmt.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'ptdat6.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   dynflg
      integer   ndf,ndm,nst,isw,i,i1,is,j,jj,j1,js,l,ni,nn,nhv
      integer   istrt
      integer   ix(*)
      real*8    bdb,bd3,dmas0, rr
      real*8    cfac,lfac,xx1,xx2,xx3,qfact, tempi, xlamd, ha, ta
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),tl(*)
      real*8    r(ndf,*),r1(3,64),vl(3),weng(64),body(3),bdy(3)
      real*8    ds(6,6,5),dd(6,6),sigv(9),sigl(16,64),epsl(6,64)
      real*8    bbd(6,3),bb(6),shpr(64),dvol(64),jvol(64)
      real*8    xr(2,64),ur(2,64),psig(3)

      save

      data      xlamd,ha /  2*0.0d0 /

c     Check process isw = 1

      if(isw.eq.1) then
        return

c     No Augmented Lagrangian update for displacement model

      elseif(isw.eq.10) then

        return

      endif

c     Compute tangent stiffness and residual force vector

c     Set quadrature points and weights

      call quadr2d(d,.true.)

c     COMPUTE TANGENT STIFFNESS AND RESIDUAL FORCE VECTOR

c     Compute shape functions and derivatives in reference configuration

      do l = 1,lint

        call interp2d(l, xl,ix, ndm,nel, .false.)

c       Compute coordinates at gauss points

        xr(1,l) = 0.0d0
        xr(2,l) = 0.0d0
        ur(1,l) = 0.0d0
        ur(2,l) = 0.0d0
        do i = 1,nel
          xr(1,l) = xr(1,l) + xl(1,i)  *shp2(3,i,l)
          xr(2,l) = xr(2,l) + xl(2,i)  *shp2(3,i,l)
          ur(1,l) = ur(1,l) + ul(1,i,1)*shp2(3,i,l)
          ur(2,l) = ur(2,l) + ul(2,i,1)*shp2(3,i,l)
        end do ! i

c       Axisymmetric volume

        jvol(l) = jac(l)
        if(stype.eq.3 .or. stype.eq.8) then
          jvol(l) = jvol(l)*xr(1,l)
        endif
      end do ! l
      xref(3) = 0.0d0
      xcur(3) = 0.0d0

c     Compute deformation gradient and determinant; transform shape
c     functions to current configuration.

      if(isw.eq.14) then
        call pfinit(f,df,finv,detf, lint)
      else
        call kine2d(shp2,xl,ul,f,finv,df,detf,ndm,ndf,nel,nen,lint)
      endif

      nhv   = nint(d(15))
      istrt = nint(d(84))
      ni    = 0

c     Set loop limits and consistent/lumped mass factor
      if(stype.eq.8) then
        is = 3
        js = 6
        cfac = 1.0d0 ! Axisymmetric + torsion must be consistent mass
      else
        is = 2
        js = 4
        cfac = d(7)
        vl(3) = 0.0d0
      endif
      lfac = 1.0d0 - cfac

c     Compute proper rank mass effects for 3-node triangle

      if(cfac.gt.0.0d0 .and. nel.eq.3) then
        call masst3(stype,cfac,d(4),xl,ul(1,1,5),r,s)
        cfac = 0.0d0
      endif

c     Transfer for output related values

      if(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16 .or. isw.eq.25 ) go to 4

c     Compute body forces

      call sbodyf(d, body)
      do i = 1,3
        bdy(i) = body(i)
      end do ! i

c     LOOP OVER GAUSS POINTS

      dynflg = ctan(3).ne.0.0d0
      nn     = ni
      do l = 1,lint

c       Set reference and current coordinates

        xref(1) = xr(1,l)
        xref(2) = xr(2,l)
        xcur(1) = xr(1,l) + ur(1,l)
        xcur(2) = xr(2,l) + ur(2,l)

c       Check for axisymmetry

        if(stype.eq.3 .or. stype.eq.8) then
          do i = 1,nel
            shpr(i) = shp2(3,i,l)/xcur(1)
          end do ! i
        else
          do i = 1,nel
            shpr(i) = 0.0d0
          end do ! i

c         Radial body loading

          if(nint(d(69)).eq.5) then
            rr = sqrt(xref(1)**2 + xref(2)**2)
            if(rr.ne.0.0d0) then
              bdy(1) = (xref(1)*body(1) - xref(2)*body(2))/rr
              bdy(2) = (xref(2)*body(1) + xref(1)*body(2))/rr
            endif
          endif
        end if

c       Compute temperature at point

        ta = 0.0d0
        do i = 1,nel
          ta = ta + tl(i)*shp2(3,i,l)
        end do ! i

c       Compute Cauchy stresses and spatial tangent tensor

        call modlfd(l,d,f(1,1,l),finv(1,l),df(1,l),detf(1,l),ta,
     &             hr(nn+nh1),hr(nn+nh2),nhv,istrt, ds,sigv,bb,
     &             xlamd,ha,.false.,isw)

c       Compute volume and mass factor

        dvol(l) = jvol(l)*detf(1,l)
        if((d(7).ge.0.0d0 .or. d(183).ne.0.0d0) .and.
     &             (ndfo(1).gt.0 .or. shflg)) then
          dmas0 = jvol(l)*d(4)
        else
          dmas0 = 0.0d0   ! No inertia effects
        endif

        if(isw.eq.13) then

          epl(8) = epl(8) + estore*jvol(l)

c         Compute velocity at point

          do i = 1,is
            vl(i) = 0.0d0
            do j = 1,nel
              vl(i) = vl(i) + ul(i,j,4)*shp2(3,j,l)
            end do ! j
          end do ! i

          tempi = 0.0d0
          if(stype.eq.8) then
            do i = 1,nel
              tempi = tempi
     &          + (ul(1,i,4)**2+ul(2,i,4)**2+ul(3,i,4)**2)*shp2(3,i,l)
            end do ! i
          else
            do i = 1,nel
              tempi = tempi
     &              + (ul(1,i,4)**2 + ul(2,i,4)**2)*shp2(3,i,l)
            end do ! i
          endif

c         Accumulate kinetic energy

          epl(7) = epl(7) + 0.5d0*(lfac*tempi
     &                    + cfac*(vl(1)**2 + vl(2)**2 + vl(3)**2))*dmas0

        elseif(isw.ne.14) then

c         Store stress values for tplot

          j1 = 6*(l-1)
          do j = 1,js
            tt(j+j1) = sigv(j)
          end do ! j

c         Multiply tangent moduli and stresses by volume element.

          do i = 1,js
            sigv(i) = sigv(i)*dvol(l)
            do j = 1,js
              dd(i,j) = ds(i,j,1)*dvol(l)*ctan(1)
            end do ! j
          end do ! i

c         Accumulate averaged thickness stress

          v_avg  = v_avg  + dvol(l)
          v_rho  = v_rho  + dvol(l)*d(4)
          sig_33 = sig_33 + sigv(3)

c         COMPUTE STRESS DIVERGENCE AND INERTIA TERMS

          do i = 1,nel

c           Stress divergence term (used in geometric stiffness)

            r1(1,i) = shp2(1,i,l)*sigv(1) + shp2(2,i,l)*sigv(4)
            r1(2,i) = shp2(1,i,l)*sigv(4) + shp2(2,i,l)*sigv(2)

c           Element residual

            r(1,i) = r(1,i) - r1(1,i) - shpr(i)*sigv(3)
     &                      + bdy(1)*jvol(l)*shp2(3,i,l)
            r(2,i) = r(2,i) - r1(2,i)
     &                      + bdy(2)*jvol(l)*shp2(3,i,l)
          end do ! i

c         Torsion residual

          if(stype.eq.8) then
            do i = 1,nel
              r1(3,i) = xcur(1)*(shp2(1,i,l)*sigv(6)
     &                         + shp2(2,i,l)*sigv(5))
              r(3,i)  = r(3,i) - r1(3,i) + bdy(3)*shp2(3,i,l)*jvol(l)
              r1(3,i) = r1(3,i)*2.d0 ! Geometric term only
            end do ! i
          endif

c         COMPUTE K (s(nst,nst) = K)

          if(isw.eq.3) then

c           PART 1. - Geometric part.

            if(gflag) then
              i1  = 0
              do i = 1,nel
                bd3 = shpr(i)*sigv(3)*ctan(1)
                j1  = 0
                do j = 1,nel
                  bdb          = (r1(1,i)*shp2(1,j,l)
     &                         +  r1(2,i)*shp2(2,j,l))*ctan(1)
                  s(i1+1,j1+1) = s(i1+1,j1+1) + bdb + bd3*shpr(j)
                  s(i1+2,j1+2) = s(i1+2,j1+2) + bdb
                  if(stype.eq.8) then
                    s(i1+1,j1+3) = s(i1+1,j1+3)+shpr(i)*r1(3,j)*ctan(1)
                    s(i1+3,j1+1) = s(i1+3,j1+1)+shpr(j)*r1(3,i)*ctan(1)
                    s(i1+3,j1+3) = s(i1+3,j1+3)+xcur(1)*xcur(1)*bdb
                  endif
                  j1 = j1 + ndf
                end do ! j
                i1 = i1 + ndf
              end do ! i
            endif ! gflag

c           PART 2. - Tangent modulus part (based upon dd-array)

            i1 = 0
            do i  = 1,nel

c             Compute bmat-t * dd * dvol

              do jj = 1,js
                bbd(jj,1) = shp2(1,i,l)*dd(1,jj)
     &                    + shpr(  i  )*dd(3,jj)
     &                    + shp2(2,i,l)*dd(4,jj)

                bbd(jj,2) = shp2(1,i,l)*dd(4,jj)
     &                    + shp2(2,i,l)*dd(2,jj)
              end do ! jj

              if(stype.eq.8) then
                do jj = 1,js
                  bbd(jj,3) = xcur(1)*(shp2(2,i,l)*dd(5,jj)
     &                               + shp2(1,i,l)*dd(6,jj))
                end do ! jj
              endif

c             Compute tangent stiffness

              j1 = 0
              do j  = 1,nel

                do jj = 1,is
                  s(i1+jj,j1+1) = s(i1+jj,j1+1) + bbd(1,jj)*shp2(1,j,l)
     &                                          + bbd(3,jj)*shpr(  j  )
     &                                          + bbd(4,jj)*shp2(2,j,l)

                  s(i1+jj,j1+2) = s(i1+jj,j1+2) + bbd(4,jj)*shp2(1,j,l)
     &                                          + bbd(2,jj)*shp2(2,j,l)
                end do ! jj

c               Torsion part

                if(stype.eq.8) then
                  do jj = 1,3
                    s(i1+jj,j1+3) = s(i1+jj,j1+3)
     &                            + (bbd(5,jj)*shp2(2,j,l)
     &                            +  bbd(6,jj)*shp2(1,j,l))*xcur(1)
                  end do ! jj

                endif

                j1 = j1 + ndf
              end do ! j

              i1 = i1 + ndf
            end  do ! i

          endif ! end of tangent

c         Add inertia parts

          if(dynflg) then
            call fdyn2d(ul,shp2(1,1,l),s,r,is,xcur(1),
     &                  cfac,lfac,dmas0,isw)
          endif

        endif ! end of isw options

        nn = nn + nhv

      end do ! l

c     Multiply by thickness if not unity

      if((isw.eq.3 .or. isw.eq.6) .and. d(14).ne.1.d0) then

        do j = 1,nst
          do i = 1,nst
            s(i,j) = s(i,j)*d(14)
          end do ! i
        end do ! j
        do j = 1,nel
          do i = 1,ndf
            r(i,j) = r(i,j)*d(14)
          end do ! i
        end do ! j

      endif

      return

c     OUTPUT STRESSES

   4  xx1  = 0.d0
      xx2  = 0.d0
      xx3  = 0.d0
      do i = 1,6
        sigv(i) = 0.0d0
        ebig(i) = 0.0d0
        esml(i) = 0.0d0
      end do ! i
      qfact = 1.d0/dble(lint)

c     LOOP OVER GAUSS POINTS

      nn = ni

      do l = 1,lint

c       Set reference and current coordinates

        xref(1) = xr(1,l)
        xref(2) = xr(2,l)
        xcur(1) = xr(1,l) + ur(1,l)
        xcur(2) = xr(2,l) + ur(2,l)

        xx1     = xx1 + xref(1)*qfact
        xx2     = xx2 + xref(2)*qfact

c       Compute temperature at point

        ta = 0.0d0
        do i = 1,nel
          ta = ta + tl(i)*shp2(3,i,l)
        end do ! i

c       Compute Cauchy stresses and spatial tangent tensor at t-n+1

        call modlfd(l,d,f(1,1,l),finv(1,l),df(1,l),detf(1,l),ta,
     &             hr(nn+nh1),hr(nn+nh2),nhv,istrt,ds,sigl(1,l),bb,
     &             xlamd,ha,.false.,isw)
        weng(l) = estore
        dvol(l) = jvol(l)*detf(1,l)

c       Compute Green-Lagrange strains and Almansi strains

        call fstrain(f(1,1,l),finv(1,l), egreen, ealmansi)

        do i = 1,6
          epsl(i,l) = ealmansi(i)
        end do ! i

c       Compute average stresses and jacobian for printing

        do i = 1,js
          sigv(i) = sigv(i) + qfact * sigl(i,l)
          ebig(i) = ebig(i) + qfact * egreen(i)
          esml(i) = esml(i) + qfact * ealmansi(i)
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

c       Output quadrature point values

        if(qoutfl) then

          do l = 1,lint
            xx1 = xr(1,l)
            xx2 = xr(2,l)
            call pstr3d(sigl(1,l),psig)
            write(iow,2002) n,ma,xx1,xx2,xx3,psig,(sigl(jj,l),jj=1,6),
     &                     (epsl(jj,l),jj=1,6)
            if(ior.lt.0) then
              write(*,2002) n,ma,xx1,xx2,xx3,psig,(sigl(jj,l),jj=1,6),
     &                     (epsl(jj,l),jj=1,6)
            endif
          end do ! l

c       Output averaged values

        else

          call pstr3d(sigv,psig)

          write(iow,2002) n,ma,xx1,xx2,xx3,psig,(sigv(jj),jj=1,6),esml

          if(ior.lt.0) then
            write(*,2002) n,ma,xx1,xx2,xx3,psig,(sigv(jj),jj=1,6),esml
          end if

        endif

c     Project stress values to nodes

      elseif(isw.eq.8) then

        call slcn2d(ix,sigl,epsl,r,s,r(nen+1,1),nel,16)

c     Compute fracture indices

      elseif(isw.eq.16) then

        call pfrac2f(f,detf,sigl,weng, shp2,dvol, r, lint,ndf,ndm,3)

c     Compute Z-Z projections

      elseif(isw.eq.25) then

        call stcn2z(xl,sigl,epsl,shp2,dvol,lint,ndm,nel,16)

      end if

c     Format statements

2001  format(a1,20a4//5x,'Element Stresses & Strains'//'   Elmt Matl',
     &   '   1-coord    2-coord    3-coord ',
     &   '   1-stress   2-stress   3-stress'/4x,'Cauchy  '
     &   '  11-stress  22-stress  33-stress  12-stress',
     &   '  23-stress  13-stress'/4x,'Almansi ',
     &   '  11-Strain  22-Strain  33-Strain  12-Strain',
     &   '  23-Strain  31-Strain')

2002  format(/i8,i4,1p,6e11.3/(12x,1p,6e11.3:))

      end
