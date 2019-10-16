      subroutine elmt25(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Bogner-Fox-Schmit Plate Bending Element for FEAP
c                Degrees of Freedom:  w, theta_x, theta_y, w_xy
c                                     theta_x = w_y; theta_y = - w_x

c     Control data:  ndm = 2(3); ndf = 4; nen = 4.

c     Inputs:    MATErial #
c                 E, nu, rho, quad pts arrays(3,4), quad pts output(3,4)
c                 h, q, alpha, T_0 

c     Reference: F.K. Bogner, R.L. Fox and L.A. Schmit
c                The generation of interelement-compatible stiffness
c                and mass matrices by the use of interpolation formulae.
c                Proc. 1st Conf. Matrix Methods in Structural Mechanics
c                Volume: AFFDITR-66-80, pp 397 ff.
c
c     modified for FEAP UKA WW 1/07 
c
c     open:
c     mass+geom matrix not tested
c     ac\
c     dr, ctan\
c     shear forces not coded, from derivatives of N!
c     nx,ny,nxy?
c     set of boundary cond. w_xy?
c     w_x? w_y?
c     alpha T_0?
c     choice of Gauss Points?
c
c-----[--+---------+---------+---------+---------+---------+---------+-]
      USE bdata      
      USE cdata
      USE eldata
      USE evdata
      USE iofile
      USE pdata7
      USE strnam 

      implicit  none
      integer iplma

      character wd*12 

      integer   ndf,ndm,nst,isw
      integer   i,j,k,l,ii,i1,jj,j1,lint

      real*8    e,xnu,alp,t0,dv
      real*8    dr
      real*8    b1
      real*8    xx,yy

      integer   ix(*),ia(16)

      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),tl(*),s(nst,*),p(*)
      real*8    shp(4,10,4,16),sg(3,16),sig(7),eps(3),td(6)
cww   real*8    ac(4)
      real*8    dd(3,3),b(3,16),bl(16),bdi(3),xsj(16),shp4(4)
      real*8    nm(16)
      real*8    h1(*),h2(*),h3(*)

      save

      data      wd/'Plate Bend.'/

c     Go to correct array processor

      if(isw.eq.0 .and. ior.lt.0) then
        write(*,*) '   Elmt  25:  Bogner-Fox Rectangular Plate Bending'

c     Input record 1 of material properties

      elseif(isw.eq.1) then

        if(ior.lt.0) write(*,3000)
        call dinput(td,6)

c       Move properties

        e    = td(1)
        xnu  = td(2)
        d(4) = td(3)
        l    = td(4)
        k    = td(5)

c       Input record 2 of material properties

        if(ior.lt.0) write(*,3001)
        call dinput(td,5)

c       Move properties

        d(14) = td(1)
        d(11) = td(2)
        alp   = td(3)
        t0    = td(4)
        if(d(14).eq.0.0d0) d(14) = 1.0d0
        l = min(4,max(3,l))
        k = min(4,max(3,k))

        if(ior.lt.0) write(*,3002)
        call dinput(td,5)

        d(15) = td(1)
        d(16) = td(2)
        d(17) = td(3)

c       Output current parameters

        write(iow,2000) wd,e,xnu,d(4),l,k,d(14),d(11),alp,t0,
     &                  d(15),d(16),d(17)
        if(ior.lt.0) then
          write(*,2000) wd,e,xnu,d(4),l,k,d(14),d(11),alp,t0,
     &                  d(15),d(16),d(17)
        endif

c       Set properties into d-array for subsequent use

        d(1)  = e/12.d0/(1.d0-xnu*xnu)*d(14)**3
        d(2)  = xnu*d(1)
        d(3)  = e/24.d0/(1.d0+xnu)*d(14)**3
        d(4)  = d(4)*d(14)
        d(5)  = l
        d(6)  = k
        d(9)  = t0
        d(10) = alp
        lint  = 0

        call pzero(dd,9)
        call pzero(b ,48)
        call pzero(bl,16)

c       Set assembly vector

        do i = 1,4
          ia(i   ) = i
          ia(i+ 4) = i + ndf
          ia(i+ 8) = i + ndf*2
          ia(i+12) = i + ndf*3
        end do ! i

        if(ndm.eq.3) ipla = 1

c....   description of stresses  
        strsus( 1) = '  MOMENT m_xx  '
        strsus( 2) = '  MOMENT m_xy  '
        strsus( 3) = '  MOMENT m_yy  '
        strsus( 4) = '               '
        strsus( 5) = '  MOMENT m_1   '
        strsus( 6) = '  MOMENT m_2   '
        strsus( 7) = '  ANGLE Phi_1  '
        do i = 8,25
          strsus(i) = '               '
        end do

c     Check element for errors in input data

      elseif(isw.eq.2) then

c       call ckisop(ix,xl,shp,ndm)
        return

c     Compute stress-divergence vector (p) and stiffness matrix (s)

      elseif(isw.eq.3 .or. isw.eq.6) then

c       Compute gauss quadrature points and weights

        l    = d(5)
        call int2d(l,lint,sg)

c       Compute shape functions for all gauss points

        call shapbf(sg,xl,shp,xsj,ndm,lint)
        call pzero(shp4,4)

c       Stiffness and residual computation

        do l = 1,lint

c         Compute stresses and strains

          call strs25(d,xl,ul,tl,shp(1,1,1,l),shp4,ndf,ndm,
     &                xx,yy,eps,sig,dd)

          dv = xsj(l)*sg(3,l)
cww       dr = d(4)*ctan(3)*xsj(l)  !??????????????
          b1 = d(11)*dv  

c         Material moduli and stresses multiplied by area

          do i = 1,3
            sig(i) = sig(i)*dv
            do j = 1,3
              dd(i,j) = dd(i,j)*dv
            end do ! j
          end do ! i

c         Form B-matrix  and N-matrix for bending element

          i1 = 1
          do i = 1,4
            call bmat25( b(1,i1),nm(i1),shp(1,1,i,l))
            do j = 1,4
              bl(i1+j-1) = shp(j,1,i,l)
            end do ! j
            i1 = i1 + 4
          end do ! i

cwwc         Compute acceleration * density * area
cww
cww          do i = 1,4
cww            ac(i) = 0.0d0
cww            do j = 1,4
cww              ac(i) = ac(i) + shp(i,1,j,l)*ul(i,j,5)
cww            end do ! j
cww            ac(i) = ac(i)*d(4)*dv
cww          end do ! i

c         Loop over rows

          do i = 1,16

            p(ia(i)) = p(ia(i)) - b(1,i)*sig(1)
     &                          - b(2,i)*sig(2)
     &                          - b(3,i)*sig(3)
cww  &                          - shp(1,1,i,l)*ac(1)
cww  &                          - shp(2,1,i,l)*ac(2)
cww  &                          - shp(3,1,i,l)*ac(3)
cww  &                          - shp(4,1,i,l)*ac(4)
     &                          + bl(i)*b1

            if(isw.eq.3) then

c             Compute B_trans * D * j * w

              do k = 1,3
                bdi(k) = b(1,i)*dd(1,k)+b(2,i)*dd(2,k)+b(3,i)*dd(3,k)
              end do ! k

c             Loop over columns (symmetry noted)

              do j = 1,16
                s(ia(i),ia(j)) = s(ia(i),ia(j)) + bdi(1)*b(1,j)
     &                                          + bdi(2)*b(2,j)
     &                                          + bdi(3)*b(3,j)
cww  &                                          + nm(i)*nm(j)*dr
              end do ! j

            end if

          end do ! i
        end do ! l

      elseif(isw.eq.4) then

c       Compute gauss quadrature points and weights

        l = d(6)
        call int2d(l,lint,sg)

c       Compute shape functions for all gauss points

        call shapbf(sg,xl,shp,xsj,ndm,lint)

        do l = 1,lint

c         Compute stresses and strains

          shp4(1) = 0.25d0*(1.d0 - sg(1,l))*(1.d0 - sg(2,l))
          shp4(2) = 0.25d0*(1.d0 + sg(1,l))*(1.d0 - sg(2,l))
          shp4(3) = 0.25d0*(1.d0 + sg(1,l))*(1.d0 + sg(2,l))
          shp4(4) = 0.25d0*(1.d0 - sg(1,l))*(1.d0 + sg(2,l))

          call strs25(d,xl,ul,tl,shp(1,1,1,l),shp4,ndf,ndm,
     &                xx,yy,eps,sig,dd)
          sig(4) = sig(3)
          sig(3) = 0.0d0
cww       call pstr2d(sig,sig(5))
          call pstres(sig,sig(5),sig(6),sig(7))

c         Output stresses and strains

          mct = mct - 2
          if(mct.le.0) then
          write(iow,2001) o,head
            if(ior.lt.0) then
              write(*,2001) o,head
            endif
            mct = 50
          endif
          write(iow,2002)  n,ma,xx,yy, sig(1), sig(2), sig(4),
     &                     (sig(i),i=5,7)
          if(ior.lt.0) then
            write(*,2002)  n,ma,xx,yy, sig(1), sig(2), sig(4),
     &                     (sig(i),i=5,7)
          endif

        end do ! l

c     Mass and geometric stiffness matrix

      elseif(isw.eq.5) then

c       Compute gauss quadrature points and weights

        l    = d(5)
        call int2d(l,lint,sg)

c       Compute shape functions for all gauss points

        call shapbf(sg,xl,shp,xsj,ndm,lint)

c       Mass

        if    (imtyp.eq.1) then

          do l = 1,lint
  
            dr = d(4)*xsj(l)

c           Compute lumped mass

            shp4(1) = 0.25d0*(1.d0 - sg(1,l))*(1.d0 - sg(2,l))
            shp4(2) = 0.25d0*(1.d0 + sg(1,l))*(1.d0 - sg(2,l))
            shp4(3) = 0.25d0*(1.d0 + sg(1,l))*(1.d0 + sg(2,l))
            shp4(4) = 0.25d0*(1.d0 - sg(1,l))*(1.d0 + sg(2,l))

            i1 = 1
            do i = 1,4
              p(i1) = p(i1) + shp4(i)*dr
              i1 = i1 + ndf
            end do ! i

c           Form N-matrix for bending element

            i1 = 1
            do i = 1,4
              call bmat25(b(1,i1),nm(i1),shp(1,1,i,l))
              i1 = i1 + 4
            end do ! i

c           Compute consistent mass

            do i = 1,16
              do j = 1,16
                s(ia(i),ia(j)) = s(ia(i),ia(j)) + nm(i)*nm(j)*dr
              end do ! j
            end do ! i

          end do ! l

c       Geometric stiffness

        elseif(imtyp.eq.2) then

c         Stiffness computation

          do l = 1,lint
  
            dv = xsj(l)*sg(3,l)

c           Initial stress matrix

            dd(1,1) = d(15)*dv
            dd(2,2) = d(16)*dv
            dd(1,2) = d(17)*dv
            dd(2,1) = d(17)*dv

            i1 = 0
            do ii = 1,4

              do i = 1,4
                b(1,i) = dd(1,1)*shp(i,2,ii,l) + dd(1,2)*shp(i,3,ii,l)
                b(2,i) = dd(2,1)*shp(i,2,ii,l) + dd(2,2)*shp(i,3,ii,l)
              end do !i

              j1 = 0
              do jj = 1,4
                do j = 1,4
                  do i = 1,4
                    s(i+i1,j+j1) = s(i+i1,j+j1) + b(1,i)*shp(j,2,jj,l)
     &                                          + b(2,i)*shp(j,3,jj,l)
                  end do ! i
                end do ! j
                j1 = j1 + ndf
              end do ! jj
              i1 = i1 + ndf
            end do ! ii

          end do ! l

        endif

c     Lumped projection of moments

      elseif(isw.eq.8) then

        istv = 4
        if(iplma(ma).eq.0)  return ! only if MATN
        call stcn25(ix,d,xl,ul,tl,shp,xsj,
     &              strea,strea(1+numnp),ndf,ndm,numnp)


      endif

c     Formats for input-output

2000  format(/5x,a12,' Linear Elastic Element'//
     1 10x,'Modulus  ',e16.5/10x,'Poisson ratio',f8.5/
     2 10x,'Density  ',e16.5/10x,'Gauss pts/dir',i3/
     3 10x,'Stress pts  ',i4/10x,'Thickness',e16.5/
     3 10x,'Pressure ',e16.5/
     4 10x,'Alpha    ',e16.5/10x,'Base temp',e16.5/
     5 10x,'N_xx     ',e16.5/10x,'N_yy     ',e16.5/
     6 10x,'N_xy     ',e16.5/)

2001  format(a1,20a4//5x,'Bogner Fox Plate Element Stresses'//
     &   '  Elem Mate',
     &   '  1-coord  2-coord',   
     &   '   11-moment   12-moment   22-moment',
     &   '    1-moment    2-moment       angle')

2002  format(i6,i5,2f9.3,5e12.3,f12.2)

3000  format(' Input: e, nu, rho, #-pts K, #-pts sig'/'   >',$)

3001  format(' Input: thick, pressure, alpha, T0'    /'   >',$)

3002  format(' Input: N_xx , N_yy , N_xy (buckling)' /'   >',$)

      end

      subroutine strs25(d,xl,ul,tl,shp,shp4,ndf,ndm, xx,yy,eps,sig,dd)
c-----[--+---------+---------+---------+---------+---------+---------+-]
c
c     calculate strains and moments
c
c-----[--+---------+---------+---------+---------+---------+---------+-]

      implicit none

      integer  ndf,ndm, i,j
      real*8   xx,yy,ta
      real*8   d(*), xl(ndm,*), ul(ndf,*), tl(*), shp(4,10,*)
      real*8   eps(4), sig(4), dd(3,3), shp4(4)

      save

c     Compute strains and coordinates

      do j = 1,3
        eps(j) = 0.0
      end do ! j

      xx = 0.0
      yy = 0.0
      ta = -d(9)
      do j = 1,4
        do i = 1,4
          eps(1) = eps(1) + shp(i,4,j)*ul(i,j)
          eps(2) = eps(2) + shp(i,5,j)*ul(i,j)
          eps(3) = eps(3) + shp(i,6,j)*ul(i,j)
        end do ! i
        xx = xx + shp4(j)*xl(1,j)
        yy = yy + shp4(j)*xl(2,j)
        ta = ta + shp4(j)*tl(1)
      end do ! j

      eps(3) = eps(3)*2.d0

c     Moduli

      dd(1,1) = d(1)
      dd(1,2) = d(2)
      dd(2,1) = d(2)
      dd(2,2) = d(1)
      dd(3,3) = d(3)

c     Compute stresses

      sig(1) = dd(1,1)*(eps(1) - d(10)*ta) + dd(1,2)*(eps(2) - d(10)*ta)
      sig(2) = dd(2,1)*(eps(1) - d(10)*ta) + dd(2,2)*(eps(2) - d(10)*ta)
      sig(3) = dd(3,3)*eps(3)

      end

      subroutine bmat25(b,nm,shp)
c-----[--+---------+---------+---------+---------+---------+---------+-]
c
c     Strain-displacement matrix
c
c-----[--+---------+---------+---------+---------+---------+---------+-]

      implicit none

      integer  j
      real*8   b(3,4), nm(4), shp(4,*)

      do j = 1,4
        b(1,j) = shp(j,4)
        b(2,j) = shp(j,5)
        b(3,j) = shp(j,6)*2.d0
        nm(j)  = shp(j,1)
      end do ! j

      end

      subroutine shapbf(sg,xl,shp,xsj,ndm,lint)
c-----[--+---------+---------+---------+---------+---------+---------+-]
c
c     shape functions and derivatives for Bogner/Fox/Schmit plate
c
c-----[--+---------+---------+---------+---------+---------+---------+-]
      USE eldata
      implicit none
  
      integer  i, j, l, lint, ndm
      real*8   xsj0, detxs, stxy, xyd, dshp, shps(4,4), shpt(4,4)
      real*8   sg(3,*),xl(ndm,*), shp(4,10,4,*),xsj(*)
      real*8   xs(2,2), sx(2,2), ssx(2,2)

c     Coordinate derivatives

      xs(1,1) = 0.25d0*(-xl(1,1)+xl(1,2)+xl(1,3)-xl(1,4))
      xs(1,2) = 0.25d0*(-xl(1,1)-xl(1,2)+xl(1,3)+xl(1,4))
      xs(2,1) = 0.25d0*(-xl(2,1)+xl(2,2)+xl(2,3)-xl(2,4))
      xs(2,2) = 0.25d0*(-xl(2,1)-xl(2,2)+xl(2,3)+xl(2,4))

      xsj0    = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)
      detxs   = 1.d0/xsj0

      sx(1,1) =  xs(2,2)*detxs
      sx(2,2) =  xs(1,1)*detxs
      sx(1,2) = -xs(1,2)*detxs
      sx(2,1) = -xs(2,1)*detxs

      xyd     =  xs(1,1)*xs(2,2) + xs(1,2)*xs(2,1)
      stxy    =  1.d0/xyd
      detxs   =  detxs*stxy

      ssx(1,1)=  xs(2,2)*xs(2,2)*detxs
      ssx(2,2)=  xs(1,1)*xs(1,1)*detxs
      ssx(1,2)= -xs(1,2)*xs(1,2)*detxs
      ssx(2,1)= -xs(2,1)*xs(2,1)*detxs

c     Form beam shape functions

      do l = 1,lint
        call beams(sg(1,l),shps)
        call beams(sg(2,l),shpt)
      
c       Set jacobian

        xsj(l) = xsj0

c       Shape functions

        shp(1,1,1,l) =  shps(1,1)*shpt(1,1)
        shp(2,1,1,l) =  shps(1,3)*shpt(1,1)
        shp(3,1,1,l) =  shps(1,1)*shpt(1,3)
        shp(4,1,1,l) =  shps(1,3)*shpt(1,3)
      
        shp(1,1,2,l) =  shps(1,2)*shpt(1,1)
        shp(2,1,2,l) =  shps(1,4)*shpt(1,1)
        shp(3,1,2,l) =  shps(1,2)*shpt(1,3)
        shp(4,1,2,l) =  shps(1,4)*shpt(1,3)
      
        shp(1,1,3,l) =  shps(1,2)*shpt(1,2)
        shp(2,1,3,l) =  shps(1,4)*shpt(1,2)
        shp(3,1,3,l) =  shps(1,2)*shpt(1,4)
        shp(4,1,3,l) =  shps(1,4)*shpt(1,4)
      
        shp(1,1,4,l) =  shps(1,1)*shpt(1,2)
        shp(2,1,4,l) =  shps(1,3)*shpt(1,2)
        shp(3,1,4,l) =  shps(1,1)*shpt(1,4)
        shp(4,1,4,l) =  shps(1,3)*shpt(1,4)
      
c       xi-derivative of shape functions

        shp(1,2,1,l) =  shps(2,1)*shpt(1,1)
        shp(2,2,1,l) =  shps(2,3)*shpt(1,1)
        shp(3,2,1,l) =  shps(2,1)*shpt(1,3)
        shp(4,2,1,l) =  shps(2,3)*shpt(1,3)
      
        shp(1,2,2,l) =  shps(2,2)*shpt(1,1)
        shp(2,2,2,l) =  shps(2,4)*shpt(1,1)
        shp(3,2,2,l) =  shps(2,2)*shpt(1,3)
        shp(4,2,2,l) =  shps(2,4)*shpt(1,3)
      
        shp(1,2,3,l) =  shps(2,2)*shpt(1,2)
        shp(2,2,3,l) =  shps(2,4)*shpt(1,2)
        shp(3,2,3,l) =  shps(2,2)*shpt(1,4)
        shp(4,2,3,l) =  shps(2,4)*shpt(1,4)
      
        shp(1,2,4,l) =  shps(2,1)*shpt(1,2)
        shp(2,2,4,l) =  shps(2,3)*shpt(1,2)
        shp(3,2,4,l) =  shps(2,1)*shpt(1,4)
        shp(4,2,4,l) =  shps(2,3)*shpt(1,4)
      
c       eta-derivative of shape functions

        shp(1,3,1,l) =  shps(1,1)*shpt(2,1)
        shp(2,3,1,l) =  shps(1,3)*shpt(2,1)
        shp(3,3,1,l) =  shps(1,1)*shpt(2,3)
        shp(4,3,1,l) =  shps(1,3)*shpt(2,3)
      
        shp(1,3,2,l) =  shps(1,2)*shpt(2,1)
        shp(2,3,2,l) =  shps(1,4)*shpt(2,1)
        shp(3,3,2,l) =  shps(1,2)*shpt(2,3)
        shp(4,3,2,l) =  shps(1,4)*shpt(2,3)
      
        shp(1,3,3,l) =  shps(1,2)*shpt(2,2)
        shp(2,3,3,l) =  shps(1,4)*shpt(2,2)
        shp(3,3,3,l) =  shps(1,2)*shpt(2,4)
        shp(4,3,3,l) =  shps(1,4)*shpt(2,4)
      
        shp(1,3,4,l) =  shps(1,1)*shpt(2,2)
        shp(2,3,4,l) =  shps(1,3)*shpt(2,2)
        shp(3,3,4,l) =  shps(1,1)*shpt(2,4)
        shp(4,3,4,l) =  shps(1,3)*shpt(2,4)
      
c       xi-xi-derivative of shape functions

        shp(1,4,1,l) =  shps(3,1)*shpt(1,1)
        shp(2,4,1,l) =  shps(3,3)*shpt(1,1)
        shp(3,4,1,l) =  shps(3,1)*shpt(1,3)
        shp(4,4,1,l) =  shps(3,3)*shpt(1,3)
      
        shp(1,4,2,l) =  shps(3,2)*shpt(1,1)
        shp(2,4,2,l) =  shps(3,4)*shpt(1,1)
        shp(3,4,2,l) =  shps(3,2)*shpt(1,3)
        shp(4,4,2,l) =  shps(3,4)*shpt(1,3)
      
        shp(1,4,3,l) =  shps(3,2)*shpt(1,2)
        shp(2,4,3,l) =  shps(3,4)*shpt(1,2)
        shp(3,4,3,l) =  shps(3,2)*shpt(1,4)
        shp(4,4,3,l) =  shps(3,4)*shpt(1,4)
      
        shp(1,4,4,l) =  shps(3,1)*shpt(1,2)
        shp(2,4,4,l) =  shps(3,3)*shpt(1,2)
        shp(3,4,4,l) =  shps(3,1)*shpt(1,4)
        shp(4,4,4,l) =  shps(3,3)*shpt(1,4)
      
c       eta-eta-derivative of shape functions

        shp(1,5,1,l) =  shps(1,1)*shpt(3,1)
        shp(2,5,1,l) =  shps(1,3)*shpt(3,1)
        shp(3,5,1,l) =  shps(1,1)*shpt(3,3)
        shp(4,5,1,l) =  shps(1,3)*shpt(3,3)
      
        shp(1,5,2,l) =  shps(1,2)*shpt(3,1)
        shp(2,5,2,l) =  shps(1,4)*shpt(3,1)
        shp(3,5,2,l) =  shps(1,2)*shpt(3,3)
        shp(4,5,2,l) =  shps(1,4)*shpt(3,3)
      
        shp(1,5,3,l) =  shps(1,2)*shpt(3,2)
        shp(2,5,3,l) =  shps(1,4)*shpt(3,2)
        shp(3,5,3,l) =  shps(1,2)*shpt(3,4)
        shp(4,5,3,l) =  shps(1,4)*shpt(3,4)
      
        shp(1,5,4,l) =  shps(1,1)*shpt(3,2)
        shp(2,5,4,l) =  shps(1,3)*shpt(3,2)
        shp(3,5,4,l) =  shps(1,1)*shpt(3,4)
        shp(4,5,4,l) =  shps(1,3)*shpt(3,4)
      
c       xi-eta-derivative of shape functions

        shp(1,6,1,l) =  shps(2,1)*shpt(2,1)
        shp(2,6,1,l) =  shps(2,3)*shpt(2,1)
        shp(3,6,1,l) =  shps(2,1)*shpt(2,3)
        shp(4,6,1,l) =  shps(2,3)*shpt(2,3)
      
        shp(1,6,2,l) =  shps(2,2)*shpt(2,1)
        shp(2,6,2,l) =  shps(2,4)*shpt(2,1)
        shp(3,6,2,l) =  shps(2,2)*shpt(2,3)
        shp(4,6,2,l) =  shps(2,4)*shpt(2,3)
      
        shp(1,6,3,l) =  shps(2,2)*shpt(2,2)
        shp(2,6,3,l) =  shps(2,4)*shpt(2,2)
        shp(3,6,3,l) =  shps(2,2)*shpt(2,4)
        shp(4,6,3,l) =  shps(2,4)*shpt(2,4)
      
        shp(1,6,4,l) =  shps(2,1)*shpt(2,2)
        shp(2,6,4,l) =  shps(2,3)*shpt(2,2)
        shp(3,6,4,l) =  shps(2,1)*shpt(2,4)
        shp(4,6,4,l) =  shps(2,3)*shpt(2,4)

c       xi-xi-xi derivative of shape function

        shp(1,7,1,l) =  shps(4,1)*shpt(1,1)
        shp(2,7,1,l) =  shps(4,3)*shpt(1,1)
        shp(3,7,1,l) =  shps(4,1)*shpt(1,3)
        shp(4,7,1,l) =  shps(4,3)*shpt(1,3)
      
        shp(1,7,2,l) =  shps(4,2)*shpt(1,1)
        shp(2,7,2,l) =  shps(4,4)*shpt(1,1)
        shp(3,7,2,l) =  shps(4,2)*shpt(1,3)
        shp(4,7,2,l) =  shps(4,4)*shpt(1,3)
      
        shp(1,7,3,l) =  shps(4,2)*shpt(1,2)
        shp(2,7,3,l) =  shps(4,4)*shpt(1,2)
        shp(3,7,3,l) =  shps(3,2)*shpt(1,4)
        shp(4,7,3,l) =  shps(4,4)*shpt(1,4)
      
        shp(1,7,4,l) =  shps(4,1)*shpt(1,2)
        shp(2,7,4,l) =  shps(4,3)*shpt(1,2)
        shp(3,7,4,l) =  shps(4,1)*shpt(1,4)
        shp(4,7,4,l) =  shps(4,3)*shpt(1,4)

c       xi-xi-eta derivative of shape function

        shp(1,8,1,l) =  shps(3,1)*shpt(2,1)
        shp(2,8,1,l) =  shps(3,3)*shpt(2,1)
        shp(3,8,1,l) =  shps(3,1)*shpt(2,3)
        shp(4,8,1,l) =  shps(3,3)*shpt(2,3)
      
        shp(1,8,2,l) =  shps(3,2)*shpt(2,1)
        shp(2,8,2,l) =  shps(3,4)*shpt(2,1)
        shp(3,8,2,l) =  shps(3,2)*shpt(2,3)
        shp(4,8,2,l) =  shps(3,4)*shpt(2,3)
      
        shp(1,8,3,l) =  shps(3,2)*shpt(2,2)
        shp(2,8,3,l) =  shps(3,4)*shpt(2,2)
        shp(3,8,3,l) =  shps(3,2)*shpt(2,4)
        shp(4,8,3,l) =  shps(3,4)*shpt(2,4)
      
        shp(1,8,4,l) =  shps(3,1)*shpt(2,2)
        shp(2,8,4,l) =  shps(3,3)*shpt(2,2)
        shp(3,8,4,l) =  shps(3,1)*shpt(2,4)
        shp(4,8,4,l) =  shps(3,3)*shpt(2,4)

c       xi-eta-eta derivative of shape function

        shp(1,9,1,l) =  shps(2,1)*shpt(3,1)
        shp(2,9,1,l) =  shps(2,3)*shpt(3,1)
        shp(3,9,1,l) =  shps(2,1)*shpt(3,3)
        shp(4,9,1,l) =  shps(2,3)*shpt(3,3)
      
        shp(1,9,2,l) =  shps(2,2)*shpt(3,1)
        shp(2,9,2,l) =  shps(2,4)*shpt(3,1)
        shp(3,9,2,l) =  shps(2,2)*shpt(3,3)
        shp(4,9,2,l) =  shps(2,4)*shpt(3,3)
      
        shp(1,9,3,l) =  shps(2,2)*shpt(3,2)
        shp(2,9,3,l) =  shps(2,4)*shpt(3,2)
        shp(3,9,3,l) =  shps(2,2)*shpt(3,4)
        shp(4,9,3,l) =  shps(2,4)*shpt(3,4)
      
        shp(1,9,4,l) =  shps(2,1)*shpt(3,2)
        shp(2,9,4,l) =  shps(2,3)*shpt(3,2)
        shp(3,9,4,l) =  shps(2,1)*shpt(3,4)
        shp(4,9,4,l) =  shps(2,3)*shpt(3,4)

c       eta-eta-eta derivative of shape function

        shp(1,10,1,l) =  shps(1,1)*shpt(4,1)
        shp(2,10,1,l) =  shps(1,3)*shpt(4,1)
        shp(3,10,1,l) =  shps(1,1)*shpt(4,3)
        shp(4,10,1,l) =  shps(1,3)*shpt(4,3)
      
        shp(1,10,2,l) =  shps(1,2)*shpt(4,1)
        shp(2,10,2,l) =  shps(1,4)*shpt(4,1)
        shp(3,10,2,l) =  shps(1,2)*shpt(4,3)
        shp(4,10,2,l) =  shps(1,4)*shpt(4,3)
      
        shp(1,10,3,l) =  shps(1,2)*shpt(4,2)
        shp(2,10,3,l) =  shps(1,4)*shpt(4,2)
        shp(3,10,3,l) =  shps(1,2)*shpt(4,4)
        shp(4,10,3,l) =  shps(1,4)*shpt(4,4)
      
        shp(1,10,4,l) =  shps(1,1)*shpt(4,2)
        shp(2,10,4,l) =  shps(1,3)*shpt(4,2)
        shp(3,10,4,l) =  shps(1,1)*shpt(4,4)
        shp(4,10,4,l) =  shps(1,3)*shpt(4,4)

c       Transform to global rotation/twist quantities

        do i = 1,4
          do j = 1,10
            dshp         =  shp(2,j,i,l)*xs(1,1) + shp(3,j,i,l)*xs(1,2)
            shp(2,j,i,l) =  shp(2,j,i,l)*xs(2,1) + shp(3,j,i,l)*xs(2,2)
            shp(3,j,i,l) = -dshp
            shp(4,j,i,l) =  shp(4,j,i,l)*xyd
          end do ! j
        end do ! i

c       Transform to global derivatives

        do i = 1,4
          do j = 1,4

c           First derivatives

            dshp         = shp(j,2,i,l)*sx(1,1) + shp(j,3,i,l)*sx(1,2)
            shp(j,3,i,l) = shp(j,2,i,l)*sx(2,1) + shp(j,3,i,l)*sx(2,2)
            shp(j,2,i,l) = dshp

c           Second derivatives

            dshp         = shp(j,4,i,l)*ssx(1,1) + shp(j,5,i,l)*ssx(1,2)
            shp(j,5,i,l) = shp(j,4,i,l)*ssx(2,1) + shp(j,5,i,l)*ssx(2,2)
            shp(j,4,i,l) = dshp

            shp(j,6,i,l) = shp(j,6,i,l)*stxy

          end do ! j
        end do ! i

      end do ! l

      end

      subroutine beams(ss, shp)
c-----[--+---------+---------+---------+---------+---------+---------+-]
c
c     Shape functions for Hermite cubic interpolations
c
c-----[--+---------+---------+---------+---------+---------+---------+-]

      implicit none

      real*8   ss, s2, s3, shp(4,4)

      save

      s2       =  1.0d0 - ss*ss
      s3       =  ss*3.0d0

      shp(1,1) =  0.25d0*(2.d0 - s3 + ss**3)
      shp(1,2) =  0.25d0*(2.d0 + s3 - ss**3)
      shp(1,3) = -0.25d0*(ss - 1.d0)*s2
      shp(1,4) = -0.25d0*(ss + 1.d0)*s2

      shp(2,1) = -0.75d0*s2
      shp(2,2) =  0.75d0*s2
      shp(2,3) =  0.25d0*(ss - 1.d0)*(s3 + 1.d0)
      shp(2,4) =  0.25d0*(ss + 1.d0)*(s3 - 1.d0)

      shp(3,1) =  1.5d0*ss
      shp(3,2) = -1.5d0*ss
      shp(3,3) =  0.50d0*(s3 - 1.d0)
      shp(3,4) =  0.50d0*(s3 + 1.d0)

      shp(4,1) =  1.5d0
      shp(4,2) = -1.5d0
      shp(4,3) =  1.5d0
      shp(4,4) =  1.5d0

      end

      subroutine stcn25(ix,d,xl,ul,tl,shp,xsj,dt,st,ndf,ndm,numnp)
c-----[--+---------+---------+---------+---------+---------+---------+-]
c
c     Lumped projection of moments
c
c-----[--+---------+---------+---------+---------+---------+---------+-]

      implicit   none

      integer   ndf,ndm,numnp
      integer   j,l,ll,lint
      real*8    xg, xx,yy

      integer   ix(*)

      real*8    d(*),ul(ndf,*),xl(ndm,*),tl(*)
      real*8    shp(4,10,4,16),sg(3,16),sig(6),eps(3)
      real*8    dd(3,3),xsj(16)
      real*8    dt(numnp),st(numnp,*), sg0(4),tg0(4), shp4(4)

      save

      data      sg0/ -0.5d0 , 0.5d0 , 0.5d0 , -0.5d0 /
      data      tg0/ -0.5d0 ,-0.5d0 , 0.5d0 ,  0.5d0 /
      data      shp4 / 4*0.0d0 /

c     Lumped projection routine

      l = d(6)

c     Compute gauss quadrature points and weights

      call int2d(l,lint,sg)

c     Compute shape functions for all gauss points

      call shapbf(sg,xl,shp,xsj,ndm,lint)

      do l = 1,lint

c       Compute stresses and strains

        shp4(1) = 0.25d0*(1.d0 - sg(1,l))*(1.d0 - sg(2,l))
        shp4(2) = 0.25d0*(1.d0 + sg(1,l))*(1.d0 - sg(2,l))
        shp4(3) = 0.25d0*(1.d0 + sg(1,l))*(1.d0 + sg(2,l))
        shp4(4) = 0.25d0*(1.d0 - sg(1,l))*(1.d0 + sg(2,l))

        call strs25(d,xl,ul,tl,shp(1,1,1,l),shp4,ndf,ndm,
     &              xx,yy,eps,sig,dd)

c       Compute lumped projection and assemble stress integrals

        do j = 1,4
          ll = ix(j)
          if(ll.gt.0) then
            xg     = xsj(l)*(0.5d0+sg0(j)*sg(1,l))
     &                     *(0.5d0+tg0(j)*sg(2,l))
            dt(ll) = dt(ll) + xg

c           Stress projections

            st(ll,1) = st(ll,1) + sig(1)*xg
            st(ll,2) = st(ll,2) + sig(3)*xg
            st(ll,3) = st(ll,3) + sig(2)*xg

          endif
        end do ! j
      end do ! l

      end

