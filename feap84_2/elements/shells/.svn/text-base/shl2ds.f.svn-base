c$Id:$
      subroutine shl2ds(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: SHL2DS is a 2-node element for shells

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         r(ndf,*)  - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'elplot.h'
      include   'eltran.h'
      include   'part0.h'
      include   'rdata.h'

      integer    ndf,ndm,nst,isw, a,b, i,j, i1,j1
      real*8     b11,b12,gh, d11,d12, r1,r2, rh1,rh3, dvol, ctan3
      real*8     ur,du,dw,dk, ra,za, th, hh, cosphi,sinphi
      real*8     gam,e_11,e_22,k_11,k_22, n_11,n_22,s_13,m_11,m_22
      real*8     d(*),ul(ndf,nen,*),xl(ndm,*)
      real*8     s(nst,*),r(ndf,*), ss(6,6),rr(3,2)
      real*8     d_t(5,5),aa(3,2),bb(5,3,2),bd(3,5),uu(3,2),shp(2,2)

      if(isw.eq.1) then
      else

c       Plate properties

        b11 = d(1)*d(14)/(1.d0 - d(2)*d(2))    ! Et/(1.d0 - nu*nu)
        b12 = d(2)*b11
        gh  = d(1)*d(14)/(1.d0+d(2))*0.5*d(37) ! kappa * Gt
        d11 = d(14)**2/12.d0*b11               ! D
        d12 = d(2)*d11

c       Shape functions for one point quadrature

        hh       =   sqrt((xl(1,2) - xl(1,1))**2
     &                  + (xl(2,2) - xl(2,1))**2)

        shp(1,1) = - 1.0d0/hh
        shp(2,1) =   0.5d0

        shp(1,2) =   1.0d0/hh
        shp(2,2) =   0.5d0

c       Compute local displacements at nodes

        cosphi   =  (xl(1,2) - xl(1,1))/hh
        sinphi   =  (xl(2,2) - xl(2,1))/hh

        do i = 1,2
          uu(1,i) =   cosphi*ul(1,i,1) + sinphi*ul(2,i,1)
          uu(2,i) = - sinphi*ul(1,i,1) + cosphi*ul(2,i,1)
          uu(3,i) =   ul(3,i,1)
          aa(1,i) =   cosphi*ul(1,i,5) + sinphi*ul(2,i,5)
          aa(2,i) = - sinphi*ul(1,i,5) + cosphi*ul(2,i,5)
          aa(3,i) =   ul(3,i,5)
        end do ! i

c       Compute local derivatives

        ra = (xl(1,1) + xl(1,2))*0.5d0
        za = (xl(2,1) + xl(2,2))*0.5d0
        ur = (ul(1,1,1) + ul(1,2,1))*0.5d0
        du = (uu(1,2) - uu(1,1))/hh
        dw = (uu(2,2) - uu(2,1))/hh
        dk = (uu(3,2) - uu(3,1))/hh
        th = (uu(3,1) + uu(3,2))*0.5d0

c       Compute strains

        e_11 = du
        e_22 = ur/ra
        k_11 = dk
        k_22 = th/ra
        gam  = dw + th

c       Set the volume of the element

        dvol = hh*ra

c       Compute the force resultants

        n_11 =  (b11*e_11 + b12*e_22)*dvol
        n_22 =  (b12*e_11 + b11*e_22)*dvol
        s_13 =   gh*gam *dvol
        m_11 =  (d11*k_11 + d12*k_22)*dvol
        m_22 =  (d12*k_11 + d11*k_22)*dvol

        if(isw.eq.3 .or. isw.eq.6) then

c         Store time history plot data for element

          tt(1) = n_11
          tt(2) = n_22
          tt(3) = s_13
          tt(4) = m_11
          tt(5) = m_22

c         Elastic 'tangent moduli'

          do i = 1,5
            do j = 1,5
              d_t(j,i) = 0.0d0
            end do ! j
          end do ! i
          d_t(1,1) = b11*dvol
          d_t(1,2) = b12*dvol
          d_t(2,1) = b12*dvol
          d_t(2,2) = b11*dvol
          d_t(3,3) = gh*dvol
          d_t(4,4) = d11*dvol
          d_t(4,5) = d12*dvol
          d_t(5,4) = d12*dvol
          d_t(5,5) = d11*dvol

c         Strain-displacement matrix

          do i = 1,2
            bb(1,1,i) =  shp(1,i)
            bb(1,2,i) =  0.0d0
            bb(1,3,i) =  0.0d0

            bb(2,1,i) =  cosphi/ra*shp(2,i)
            bb(2,2,i) = -sinphi/ra*shp(2,i)
            bb(2,3,i) =  0.0d0

            bb(3,1,i) =  0.0d0
            bb(3,2,i) =  shp(1,i)
            bb(3,3,i) =  shp(2,i)

            bb(4,1,i) =  0.0d0
            bb(4,2,i) =  0.0d0
            bb(4,3,i) =  shp(1,i)

            bb(5,1,i) =  0.0d0
            bb(5,2,i) =  0.0d0
            bb(5,3,i) =  shp(2,i)/ra
          end do ! i

c         Compute residual and tangent

          do a = 1,3
            rr(a,1) = 0.0d0
            rr(a,2) = 0.0d0
          end do ! a

          i1 = 0
          do i = 1,2

c           Compute the residual term

            do a = 1,3
              rr(a,i) = rr(a,i) - bb(1,a,i)*n_11
     &                          - bb(2,a,i)*n_22
     &                          - bb(3,a,i)*s_13
     &                          - bb(4,a,i)*m_11
     &                          - bb(5,a,i)*m_22
              do b = 1,5
                bd(a,b) = bb(1,a,i)*d_t(1,b)
     &                  + bb(2,a,i)*d_t(2,b)
     &                  + bb(3,a,i)*d_t(3,b)
     &                  + bb(4,a,i)*d_t(4,b)
     &                  + bb(5,a,i)*d_t(5,b)
              end do ! b
            end do ! a

            j1 = 0
            do j = 1,2
              do a = 1,3
                do b = 1,3
                  ss(i1+a,j1+b) = bd(a,1)*bb(1,b,j)
     &                          + bd(a,2)*bb(2,b,j)
     &                          + bd(a,3)*bb(3,b,j)
     &                          + bd(a,4)*bb(4,b,j)
     &                          + bd(a,5)*bb(5,b,j)
                end do ! b
              end do ! a
              j1 = j1 + 3
            end do ! j
            i1 = i1 + 3
          end do ! i

          do j = 1,6
            do i = 1,6
              ss(i,j) = ss(i,j)*ctan(1)
            end do ! i
          end do ! j

c         Inertial terms

          if(ndfo(1).gt.0 .or. shflg) then
            ctan3 = ctan(3)
          else
            ctan3 = 0.0d0
          endif
          call shms2d(d,xl,r1,r2,hh,rh1,rh3,shp,ctan3,rr,ss,ndm,
     &               .false.)

c         Inertia residual

          do j = 1,2
            rr(1,j) = rr(1,j) - (shp(i,1)*r1*aa(1,1)
     &                        +  shp(i,2)*r2*aa(1,2))*rh1
            rr(2,j) = rr(2,j) - (shp(i,1)*r1*aa(2,1)
     &                        +  shp(i,2)*r2*aa(2,2))*rh1
     &                        - (shp(i,1)*r1 + shp(i,2)*r2)*d(10)*hh
            rr(3,j) = rr(3,j) - (shp(i,1)*r1*aa(3,1)
     &                        +  shp(i,2)*r2*aa(3,2))*rh3
          end do ! j

c         Transform to global coordinates

          call trsh2d(cosphi,sinphi,rr,ss, r,s, ndf,nst)

c         Add body forces

          j = 2
          call shl2bd(d,xl,r,ndm,ndf,j)

c       Stress outputs

        elseif(isw.eq.4) then

          call shst2d(ra,za,e_11,e_22, gam,k_11,k_22,
     &                      n_11,n_22,s_13,m_11,m_22)

c       Mass

        elseif(isw.eq.5) then

c         Form a mass matrix

          call shms2d(d,xl,r1,r2,hh,rh1,rh3,shp,1.0d0,rr,ss,ndm,
     &               .true.)

c         Transform to global coordinates

          call trsh2d(cosphi,sinphi,rr,ss, r,s, ndf,nst)

c       Project stress resultants to nodes

        elseif(isw.eq.8) then
          call shcn2d(e_11,e_22, gam,k_11,k_22,
     &                n_11,n_22,s_13,m_11,m_22, r,s,nen)

        endif

      endif

      end

      subroutine trsh2d(cosphi,sinphi,rr,ss, r,s, ndf,nst)

      implicit   none

      integer    ndf,nst, i,i1,i2, j,j1,j2, a
      real*8     rr(3,*), ss(6,*), r(ndf,*), s(nst,*)
      real*8     bd(3,3), cosphi,sinphi

c     Transform to global coordinates

      i1 = 0
      i2 = 0
      do i = 1,2
        r(1,i) = cosphi*rr(1,i) - sinphi*rr(2,i)
        r(2,i) = sinphi*rr(1,i) + cosphi*rr(2,i)
        r(3,i) = rr(3,i)
        j1 = 0
        j2 = 0
        do j = 1,2
          do a = 1,3
            bd(1,a) = cosphi*ss(i1+1,j1+a) - sinphi*ss(i1+2,j1+a)
            bd(2,a) = sinphi*ss(i1+1,j1+a) + cosphi*ss(i1+2,j1+a)
            bd(3,a) = ss(i1+3,j1+a)
          end do ! a
          do a = 1,3
            s(i2+a,j2+1) = bd(a,1)*cosphi - bd(a,2)*sinphi
            s(i2+a,j2+2) = bd(a,1)*sinphi + bd(a,2)*cosphi
            s(i2+a,j2+3) = bd(a,3)
          end do ! a
          j1 = j1 + 3
          j2 = j2 + ndf
        end do ! j
        i1 = i1 + 3
        i2 = i2 + ndf
      end do ! i

      end
