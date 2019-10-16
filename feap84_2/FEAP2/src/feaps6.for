      subroutine shap2(s,t,shp,ix,nel)
c----------------------------------------------------------------------
c
c      Purpose: Adds quadratic functions to quadrilaterals for any
c               non-zero mid-side or central node
c
c      Inputs:
c         s,t      - Natural coordinates
c         ix(*)    - List of nodes attached to element (0 = no node)
c         nel      - Maximum number of local node on element <= 9
c
c      Outputs:
c         shp(3,*) - Shape functions and derivatives w/r natural coords
c                    shp(1,i) = dN_i/dxi_1
c                    shp(2,i) = dN_i/dxi_2
c                    shp(3,i) = N_i
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ix(nel),shp(3,9),shp9(3) !in case of generation nel is the number of nodes of patch
      s2 = (1.d0-s*s)/2.d0
      t2 = (1.d0-t*t)/2.d0
c...  for midside/center node
      do  i=5,9
        do  j = 1,3
          shp(j,i) = 0.0d0
        end do
      end do
c.... midside nodes (serendipity)
      if(ix(5).ne.0) then
        shp(1,5) = -s* (1.d0-t)
        shp(2,5) = -s2
        shp(3,5) =  s2*(1.d0-t)
      end if
      if(nel.lt.6) go to 100
      if(ix(6).ne.0) then
        shp(1,6) =  t2
        shp(2,6) = -t* (1.d0+s)
        shp(3,6) =  t2*(1.d0+s)
      end if
      if(nel.lt.7) go to 100
      if(ix(7).ne.0) then
        shp(1,7) = -s* (1.d0+t)
        shp(2,7) =  s2
        shp(3,7) =  s2*(1.d0+t)
      end if
      if(nel.lt.8) go to 100
      if(ix(8).ne.0) then
        shp(1,8) = -t2
        shp(2,8) = -t* (1.d0-s)
        shp(3,8) =  t2*(1.d0-s)
      end if
      if(nel.lt.9)   go to 100   ! ww separated,check bounds
      if(ix(9).eq.0) go to 100    
c.... interior node (lagrangian) ! ww pos. changed, check bounds
      shp(1,9) = -4.d0*s*t2
      shp(2,9) = -4.d0*t*s2
      shp(3,9) =  4.d0*s2*t2
c
c.... correct edge nodes for interior node (lagrangian)
      do 106 j= 1,3
        do 105 i = 1,4
105     shp(j,i) = shp(j,i) - 0.25d0*shp(j,9)
        do 106 i = 5,8
106   if(ix(i).ne.0) shp(j,i) = shp(j,i) - 0.5d0*shp(j,9)
c
c.... correct corner nodes for presense of midside nodes
100   k = 8
      do i = 1,4
        l = i + 4
        do j = 1,3
          shp(j,i) = shp(j,i) - 0.5d0*(shp(j,k)+shp(j,l))
        end do
        k = l
      end do
C
      return
      end
c
      subroutine shape(ss,tt,x,shp,xsj,ndm,nel,ix,flg)
c----------------------------------------------------------------------
c
c      Purpose: Computes shape function and derivatives for
c               quadrilateral elements
c
c      Inputs:
c         ss        - Natural coordinates for point
c         tt        - Natural coordinates for point
c         x(ndm,*)  - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c         nel       - Number of nodes on element
c         ix(*)     - Nodes attached to element
c         flg       - Flag, compute global x/y derivatives if false,
c                           else derivatives are w/r natural coords.
c
c      Outputs:
c         shp(3,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dx or dN_i/dxi_1
c                     shp(2,i) = dN_i/dy or dN_i/dxi_2
c                     shp(3,i) = N_i
c         xsj       - Jacobian determinant at point
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical flg
      dimension shp(3,nel),x(ndm,nel),s(4),t(4),xs(3,2),sx(2,2),ix(*)
      data s/-0.5d0,0.5d0,0.5d0,-0.5d0/,t/-0.5d0,-0.5d0,0.5d0,0.5d0/
c.... form 4-node quadrilateral shape functions
      if(nel.eq.4.and. .not.flg) then
          call shapef(ss,tt,x,shp,xsj,ndm,flg)
      else
          do 100 i = 1,4
            shp(3,i) = (0.5+s(i)*ss)*(0.5+t(i)*tt)
            shp(1,i) = s(i)*(0.5+t(i)*tt)
100         shp(2,i) = t(i)*(0.5+s(i)*ss)
          if(nel.ge.4) go to 120
c....   form linear triangle by adding third and fourth together
          do 110 i = 1,3
110         shp(i,3) = shp(i,3)+shp(i,4)
c....   add quadratic terms if necessary, quadratic triangle not possible!
120       if(nel.gt.4) call shap2(ss,tt,shp,ix,nel)
c....   construct jacobian and its inverse
          do 130 i = 1,ndm
           do 130 j = 1,2
            xs(i,j) = 0.0
            do 130 k = 1,nel
130          xs(i,j) = xs(i,j) + x(i,k)*shp(j,k)
        xsj = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
        if(flg) return
        if(xsj.le.0.0d0) then
          call drawmess('negative Jacobian in shape',1,-2)
          return
        end if
        sx(1,1) = xs(2,2)/xsj
        sx(2,2) = xs(1,1)/xsj
        sx(1,2) =-xs(1,2)/xsj
        sx(2,1) =-xs(2,1)/xsj
c....   form global derivatives
        do 140 i = 1,nel
            tp = shp(1,i)*sx(1,1)+shp(2,i)*sx(2,1)
            shp(2,i)  = shp(1,i)*sx(1,2)+shp(2,i)*sx(2,2)
140         shp(1,i) = tp
      end if
      return
      end
c
      subroutine shapef(s,t,xl,shp,xsj,ndm,flg)
c----------------------------------------------------------------------
c
c      Purpose: Shape function routine for 4-node isoparametric
c               quadrilaterals
c
c      Inputs:
c         s,t       - Natural coordinates of point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c         flg       - Flag, Compute global derivatives if true,
c                           else compute derivatives w/r natural coords.

c      Outputs:
c         shp(3,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dx  or dN_i/dxi_1
c                     shp(2,i) = dN_i/dy  or dN_i/dxi_2
c                     shp(3,i) = N_i
c         xsj       - Jacobian determinant at point
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical flg
      dimension xl(ndm,4),shp(3,4)
c.... set up interpolations
      sh = 0.5d0*s
      th = 0.5d0*t
      sp = 0.5d0 + sh
      tp = 0.5d0 + th
      sm = 0.5d0 - sh
      tm = 0.5d0 - th
      shp(3,1) =   sm*tm
      shp(3,2) =   sp*tm
      shp(3,3) =   sp*tp
      shp(3,4) =   sm*tp
c.... set up natural coordinate functions (times 4)
      xo =  xl(1,1)-xl(1,2)+xl(1,3)-xl(1,4)
      xs = -xl(1,1)+xl(1,2)+xl(1,3)-xl(1,4) + xo*t
      xt = -xl(1,1)-xl(1,2)+xl(1,3)+xl(1,4) + xo*s
      yo =  xl(2,1)-xl(2,2)+xl(2,3)-xl(2,4)
      ys = -xl(2,1)+xl(2,2)+xl(2,3)-xl(2,4) + yo*t
      yt = -xl(2,1)-xl(2,2)+xl(2,3)+xl(2,4) + yo*s
c.... compute jacobian (times 16)
      xsj1 = xs*yt - xt*ys
      if(xsj1.le.0.0d0) then
        call drawmess('negative Jacobian in shapef',1,-2)
        return
      end if
c.... divide jacobian by 16 (multiply by .0625)
      xsj = 0.0625d0*xsj1
      if(flg) return
      if(xsj1.eq.0.0d0) xsj1 = 1.0d0
c.... divide functions by jacobian
      xs  = (xs+xs)/xsj1
      xt  = (xt+xt)/xsj1
      ys  = (ys+ys)/xsj1
      yt  = (yt+yt)/xsj1
c.... multiply by interpolations
      ytm =  yt*tm
      ysm =  ys*sm
      ytp =  yt*tp
      ysp =  ys*sp
      xtm =  xt*tm
      xsm =  xs*sm
      xtp =  xt*tp
      xsp =  xs*sp
c.... compute shape functions
      shp(1,1) = - ytm+ysm
      shp(1,2) =   ytm+ysp
      shp(1,3) =   ytp-ysp
      shp(1,4) = - ytp-ysm
      shp(2,1) =   xtm-xsm
      shp(2,2) = - xtm-xsp
      shp(2,3) = - xtp+xsp
      shp(2,4) =   xtp+xsm
      return
      end
c
      subroutine shp2dc(ss,xl, shp, xsj, ord, flg)
c-----[--.----+----.----+----.-----------------------------------------]
c       Original version feap 8.3                            01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Shape function routine for cubic (16-node) elements
c     Inputs:
c       ss(*)    - Gauss point
c       xl(*)    - Element coordinates
c       ord      - Order to generate:
c       flg      - .false. for global derivatives

c     Outputs:
c       shp(*)   - Shape functions    (ord.ge.0) and
c                  first derivatives  (ord.ge.1)
c-----[--+---------+---------+---------+---------+---------+---------+-]
      USE sdata
      implicit   none
      logical    flg
      integer    ord, i,j,k, xi1(16),xi2(16)
      real*8     ss(2),xl(ndm,*),xds(2,2), shp(3,16), xsj
      real*8     n1(4),dn1(4), n2(4),dn2(4)
      real*8     mn1, xi1s2,xi2s2,xi1s9,xi2s9
      real*8     one9

      save

      data       xi1/1,2,2,1,3,4,2,2,4,3,1,1,3,4,4,3/
      data       xi2/1,1,2,2,1,1,3,4,2,2,4,3,3,3,4,4/

      one9=1.d0/9.d0

c     Do Shape functions

      xi1s9  = one9 - ss(1)*ss(1)
      xi2s9  = one9 - ss(2)*ss(2)
      xi1s2  = 1.d0 - ss(1)*ss(1)
      xi2s2  = 1.d0 - ss(2)*ss(2)

      n1(1)  = -9.d0*(1.d0 - ss(1))*xi1s9*0.0625d0
      n1(2)  = -9.d0*(1.d0 + ss(1))*xi1s9*0.0625d0
      n1(3)  = 27.d0*xi1s2*(1.d0/3.d0 - ss(1))*0.0625d0
      n1(4)  = 27.d0*xi1s2*(1.d0/3.d0 + ss(1))*0.0625d0

      n2(1)  = -9.d0*(1.d0 - ss(2))*xi2s9*0.0625d0
      n2(2)  = -9.d0*(1.d0 + ss(2))*xi2s9*0.0625d0
      n2(3)  = 27.d0*xi2s2*(1.d0/3.d0 - ss(2))*0.0625d0
      n2(4)  = 27.d0*xi2s2*(1.d0/3.d0 + ss(2))*0.0625d0

      dn1(1) = (  1.d0 + (18.d0 - 27.d0*ss(1))*ss(1))*0.0625d0
      dn1(2) = ( -1.d0 + (18.d0 + 27.d0*ss(1))*ss(1))*0.0625d0
      dn1(3) = (-27.d0 - (18.d0 - 81.d0*ss(1))*ss(1))*0.0625d0
      dn1(4) = ( 27.d0 - (18.d0 + 81.d0*ss(1))*ss(1))*0.0625d0

      dn2(1) = (  1.d0 + (18.d0 - 27.d0*ss(2))*ss(2))*0.0625d0
      dn2(2) = ( -1.d0 + (18.d0 + 27.d0*ss(2))*ss(2))*0.0625d0
      dn2(3) = (-27.d0 - (18.d0 - 81.d0*ss(2))*ss(2))*0.0625d0
      dn2(4) = ( 27.d0 - (18.d0 + 81.d0*ss(2))*ss(2))*0.0625d0

      do k = 1,16
        shp(3,k) = n1(xi1(k))*n2(xi2(k))
      end do ! k

c     Do first derivatives

      if(ord.ge.1) then

c       Local derivatives

        do k = 1,16
          shp(1,k) = dn1(xi1(k))* n2(xi2(k))
          shp(2,k) =  n1(xi1(k))*dn2(xi2(k))
        end do ! k

c       Jacobian matrix

        do j = 1,2
          do i = 1,2
            xds(i,j) = 0.0d0
            do k = 1,16
              xds(i,j) = xds(i,j) + xl(i,k)*shp(j,k)
            end do ! k
          end do ! i
        end do ! j
        xsj = xds(1,1)*xds(2,2) - xds(1,2)*xds(2,1)

c       Global derivatives

        if(.not.flg) then
          do k = 1,16
            mn1      = ( xds(2,2)*shp(1,k) - xds(2,1)*shp(2,k))/xsj
            shp(2,k) = (-xds(1,2)*shp(1,k) + xds(1,1)*shp(2,k))/xsj
            shp(1,k) = mn1
          end do ! k
        end if

      end if

      end
c
      subroutine shp3d(ss,xsj,shp,xl,ndm)
c----------------------------------------------------------------------
c
c      Purpose: Compute 3-d isoparametric 8-node element shape
c               functions and their derivatives w/r x,y,z
c
c      Inputs:
c         ss(3)     - Natural coordinates of point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c
c      Outputs:
c         xsj       - Jacobian determinant at point
c         shp(4,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dx
c                     shp(2,i) = dN_i/dy
c                     shp(3,i) = dN_i/dz
c                     shp(4,i) =  N_i
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      real   s1(8),s2(8),s3(8)
      real*8 ss(3),shp(4,*),xl(ndm,*),xs(3,3),ad(3,3)
      data s1/-0.5d0, 0.5d0, 0.5d0,-0.5d0,-0.5d0, 0.5d0,0.5d0,-0.5d0/
      data s2/-0.5d0,-0.5d0, 0.5d0, 0.5d0,-0.5d0,-0.5d0,0.5d0, 0.5d0/
      data s3/-0.5d0,-0.5d0,-0.5d0,-0.5d0, 0.5d0, 0.5d0,0.5d0, 0.5d0/
c.... compute shape functions and their natural coord. derivatives
      do 100 i = 1,8
        a1 = 0.5d0 + s1(i)*ss(1)
        a2 = 0.5d0 + s2(i)*ss(2)
        a3 = 0.5d0 + s3(i)*ss(3)
        shp(1,i) = s1(i)*a2*a3
        shp(2,i) = s2(i)*a1*a3
        shp(3,i) = s3(i)*a1*a2
        shp(4,i) = a1*a2*a3
100   continue
      if(ndm.lt.3) return
c.... compute jacobian transformation
      do 115 i = 1,3
      do 115 j = 1,3
        xs(i,j) = 0.0
        do 110 k = 1,8
          xs(i,j) = xs(i,j) + xl(j,k)*shp(i,k)
110     continue
115   continue
c.... compute adjoint to jacobian
      ad(1,1) = xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2)
      ad(1,2) = xs(1,3)*xs(3,2) - xs(1,2)*xs(3,3)
      ad(1,3) = xs(1,2)*xs(2,3) - xs(1,3)*xs(2,2)
      ad(2,1) = xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3)
      ad(2,2) = xs(1,1)*xs(3,3) - xs(1,3)*xs(3,1)
      ad(2,3) = xs(1,3)*xs(2,1) - xs(1,1)*xs(2,3)
      ad(3,1) = xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1)
      ad(3,2) = xs(1,2)*xs(3,1) - xs(1,1)*xs(3,2)
      ad(3,3) = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)
c.... compute determinant of jacobian
      xsj     = xs(1,1)*ad(1,1)+xs(1,2)*ad(2,1)+xs(1,3)*ad(3,1)
      if(xsj.le.0.0d0) then
        call drawmess('negative Jacobian in shp3d',1,-2)
        return
      end if
c.... compute jacobian inverse
      do 120 i = 1,3
      do 120 j = 1,3
        xs(i,j) = ad(i,j)/xsj
120   continue
c.... compute the derivatives with repect to global coords.
      do 125 k = 1,8
        c1 = xs(1,1)*shp(1,k) + xs(1,2)*shp(2,k) + xs(1,3)*shp(3,k)
        c2 = xs(2,1)*shp(1,k) + xs(2,2)*shp(2,k) + xs(2,3)*shp(3,k)
        c3 = xs(3,1)*shp(1,k) + xs(3,2)*shp(2,k) + xs(3,3)*shp(3,k)
        shp(1,k) = c1
        shp(2,k) = c2
        shp(3,k) = c3
125   continue
      return
      end
c
      subroutine shp3dc(ss, shp)
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute cubic shape functions for 64-node Lagrangian
c               element.

c      Inputs:
c         ss(*)    - Gauss points

c      Outputs:
c         shp(4,*) - Cubic shape functions

c      Original version                             FEAP 8.3 07/02/2009
c
c      nodal numbering on 64-node brick  by layers from bottom to top:
c      N.B. vertex nodes different from other bricks.

c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

cww      logical    newfl
      integer    i, ir(64), is(64), it(64), or(64), os(64), ot(64)
      real*8     nr(4),dr(4), ns(4),ds(4), nt(4),dt(4), ss(3), shp(4,*)

cww      data       newfl / .true. /

      data       ir /1,3,4,2, 1,3,4,2, 1,3,4,2, 1,3,4,2,
     &               1,3,4,2, 1,3,4,2, 1,3,4,2, 1,3,4,2,
     &               1,3,4,2, 1,3,4,2, 1,3,4,2, 1,3,4,2,
     &               1,3,4,2, 1,3,4,2, 1,3,4,2, 1,3,4,2/

      data       is /1,1,1,1, 3,3,3,3, 4,4,4,4, 2,2,2,2,
     &               1,1,1,1, 3,3,3,3, 4,4,4,4, 2,2,2,2,
     &               1,1,1,1, 3,3,3,3, 4,4,4,4, 2,2,2,2,
     &               1,1,1,1, 3,3,3,3, 4,4,4,4, 2,2,2,2/

      data       it /1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
     &               3,3,3,3, 3,3,3,3, 3,3,3,3, 3,3,3,3,
     &               4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
     &               2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2/

      data       or /1,2,2,1, 1,2,2,1, 1,2,2,1, 1,2,2,1,
     &               3,4,2,2,4,3,1,1,3,4,4,3,
     &               3,4,2,2,4,3,1,1,3,4,4,3,
     &               3,4,2,2,4,3,1,1,3,4,4,3,
     &               3,4,2,2,4,3,1,1,3,4,4,3/

      data       os /1,1,2,2, 1,1,2,2, 1,1,2,2, 1,1,2,2,
     &               1,1,3,4,2,2,4,3,3,3,4,4,
     &               1,1,3,4,2,2,4,3,3,3,4,4,
     &               1,1,3,4,2,2,4,3,3,3,4,4,
     &               1,1,3,4,2,2,4,3,3,3,4,4/

      data       ot /1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4,
     &               1,1,1,1,1,1,1,1,1,1,1,1,
     &               2,2,2,2,2,2,2,2,2,2,2,2,
     &               3,3,3,3,3,3,3,3,3,3,3,3,
     &               4,4,4,4,4,4,4,4,4,4,4,4/

c     Set 1-d shape functions for each local direction

      nr(3) =  9.0d0*(1.d0 - ss(1)*ss(1))*(1.d0 - 3.d0*ss(1))/16.d0
      nr(4) =  9.0d0*(1.d0 - ss(1)*ss(1))*(1.d0 + 3.d0*ss(1))/16.d0
      nr(1) =  0.5d0 * (1.d0 - ss(1)) - (2.d0*nr(3) + nr(4))/3.d0
      nr(2) =  0.5d0 * (1.d0 + ss(1)) - (2.d0*nr(4) + nr(3))/3.d0

      ns(3) =  9.0d0*(1.d0 - ss(2)*ss(2))*(1.d0 - 3.d0*ss(2))/16.d0
      ns(4) =  9.0d0*(1.d0 - ss(2)*ss(2))*(1.d0 + 3.d0*ss(2))/16.d0
      ns(1) =  0.5d0 * (1.d0 - ss(2)) - (2.d0*ns(3) + ns(4))/3.d0
      ns(2) =  0.5d0 * (1.d0 + ss(2)) - (2.d0*ns(4) + ns(3))/3.d0

      nt(3) =  9.0d0*(1.d0 - ss(3)*ss(3))*(1.d0 - 3.d0*ss(3))/16.d0
      nt(4) =  9.0d0*(1.d0 - ss(3)*ss(3))*(1.d0 + 3.d0*ss(3))/16.d0
      nt(1) =  0.5d0 * (1.d0 - ss(3)) - (2.d0*nt(3) + nt(4))/3.d0
      nt(2) =  0.5d0 * (1.d0 + ss(3)) - (2.d0*nt(4) + nt(3))/3.d0

      dr(3) =  9.0d0*( 9.d0*ss(1)*ss(1) - 2.d0*ss(1) - 3.d0)/16.d0
      dr(4) =  9.0d0*(-9.d0*ss(1)*ss(1) - 2.d0*ss(1) + 3.d0)/16.d0
      dr(1) = -0.5d0 - (2.d0*dr(3) + dr(4))/3.d0
      dr(2) =  0.5d0 - (2.d0*dr(4) + dr(3))/3.d0

      ds(3) =  9.0d0*( 9.d0*ss(2)*ss(2) - 2.d0*ss(2) - 3.d0)/16.d0
      ds(4) =  9.0d0*(-9.d0*ss(2)*ss(2) - 2.d0*ss(2) + 3.d0)/16.d0
      ds(1) = -0.5d0 - (2.d0*ds(3) + ds(4))/3.d0
      ds(2) =  0.5d0 - (2.d0*ds(4) + ds(3))/3.d0

      dt(3) =  9.0d0*( 9.d0*ss(3)*ss(3) - 2.d0*ss(3) - 3.d0)/16.d0
      dt(4) =  9.0d0*(-9.d0*ss(3)*ss(3) - 2.d0*ss(3) + 3.d0)/16.d0
      dt(1) = -0.5d0 - (2.d0*dt(3) + dt(4))/3.d0
      dt(2) =  0.5d0 - (2.d0*dt(4) + dt(3))/3.d0

c     Set local 3-d shape functions

cww      if(newfl) then           ! New numbering order
        do i = 1,64
          shp(1,i) = dr(ir(i)) * ns(is(i)) * nt(it(i))
          shp(2,i) = nr(ir(i)) * ds(is(i)) * nt(it(i))
          shp(3,i) = nr(ir(i)) * ns(is(i)) * dt(it(i))
          shp(4,i) = nr(ir(i)) * ns(is(i)) * nt(it(i))
        end do ! i
cww      else                     ! Old numbering order
cww        do i = 1,64
cww          shp(1,i) = dr(or(i)) * ns(os(i)) * nt(ot(i))
cww          shp(2,i) = nr(or(i)) * ds(os(i)) * nt(ot(i))
cww          shp(3,i) = nr(or(i)) * ns(os(i)) * dt(ot(i))
cww          shp(4,i) = nr(or(i)) * ns(os(i)) * nt(ot(i))
cww        end do ! i
cww      end if

      end
c
      subroutine shp3dt(ss,xsj,shp,xl,ndm,n)
c-----------------------------------------------------------------------
c
c      Purpose: Compute 3-d isoparametric 4-node tetrahedron element shape
c               functions and their derivatives wrt. x,y,z
c
c      Inputs:
c         ss(3)     - Natural coordinates of point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c         n         - Element number
c
c      Outputs:
c         xsj       - Jacobian determinant at point
c         shp(4,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dx
c                     shp(2,i) = dN_i/dy
c                     shp(3,i) = dN_i/dz
c                     shp(4,i) =  N_i
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      real*8 s1(4),s2(4),s3(4)
      real*8 ss(3),shp(4,*),xl(ndm,*),xs(3,3),ad(3,3)
      character*40 yyy
      data s1/-1.0d0, 1.0d0, 0.0d0, 0.0d0/
      data s2/-1.0d0, 0.0d0, 1.0d0, 0.0d0/
      data s3/-1.0d0, 0.0d0, 0.0d0, 1.0d0/
c.... compute shape functions and their natural coord. derivatives
      shp(1,1) = -1.d0
      shp(2,1) = -1.d0
      shp(3,1) = -1.d0
      shp(4,1) = 1.d0 - ss(1) - ss(2) - ss(3)
      do 100 i = 2,4
        shp(1,i) = s1(i)
        shp(2,i) = s2(i)
        shp(3,i) = s3(i)
        shp(4,i) = ss(i-1)
100   continue
      if(ndm.lt.3) return
c.... compute jacobian transformation
c     do 115 i = 1,3
c     do 115 j = 1,3
c       xs(i,j) = 0.0
c       do 110 k = 1,4
c         xs(i,j) = xs(i,j) + xl(j,k)*shp(i,k)
c110     continue
c115   continue
      xs(1,1) = xl(1,2) - xl(1,1)
      xs(1,2) = xl(2,2) - xl(2,1)
      xs(1,3) = xl(3,2) - xl(3,1)
      xs(2,1) = xl(1,3) - xl(1,1)
      xs(2,2) = xl(2,3) - xl(2,1)
      xs(2,3) = xl(3,3) - xl(3,1)
      xs(3,1) = xl(1,4) - xl(1,1)
      xs(3,2) = xl(2,4) - xl(2,1)
      xs(3,3) = xl(3,4) - xl(3,1)
c.... compute adjoint to jacobian
      ad(1,1) = xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2)
      ad(1,2) = xs(1,3)*xs(3,2) - xs(1,2)*xs(3,3)
      ad(1,3) = xs(1,2)*xs(2,3) - xs(1,3)*xs(2,2)
      ad(2,1) = xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3)
      ad(2,2) = xs(1,1)*xs(3,3) - xs(1,3)*xs(3,1)
      ad(2,3) = xs(1,3)*xs(2,1) - xs(1,1)*xs(2,3)
      ad(3,1) = xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1)
      ad(3,2) = xs(1,2)*xs(3,1) - xs(1,1)*xs(3,2)
      ad(3,3) = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)
c.... compute determinant of jacobian
      xsj     = xs(1,1)*ad(1,1)+xs(1,2)*ad(2,1)+xs(1,3)*ad(3,1)
      if(xsj.le.0.0d0) then
        write(yyy,'(a33,i7)') 'negative Jacobian in shp3dt in EL',n
        call drawmess(yyy,1,-2)
        return
      end if
c.... compute jacobian inverse
      do 120 i = 1,3
      do 120 j = 1,3
        xs(i,j) = ad(i,j)/xsj
120   continue
c.... compute derivatives with respect to global coordinates
      do 125 k = 1,4
        c1 = xs(1,1)*shp(1,k) + xs(1,2)*shp(2,k) + xs(1,3)*shp(3,k)
        c2 = xs(2,1)*shp(1,k) + xs(2,2)*shp(2,k) + xs(2,3)*shp(3,k)
        c3 = xs(3,1)*shp(1,k) + xs(3,2)*shp(2,k) + xs(3,3)*shp(3,k)
        shp(1,k) = c1
        shp(2,k) = c2
        shp(3,k) = c3
125   continue
      return
      end
c
      subroutine sphere(x,ndm,prt)
c----------------------------------------------------------------------
c
c      Purpose: Converts spherical to cartesian coordinates
c
c      Inputs:
c         x(ndm,*)  - Spherical coordinates of point
c         ndm       - Spatial dimension of mesh
c         prt       - Flag, output results if true

c      Outputs:
c         x(ndm,*)  - Cartesian coordinates of point
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension x(ndm,1),td(6)
      character yyy*80
      if(ndm.lt.3) then
        write(iow,4000)
        if(ior.lt.0) write(*,4000)
        return
      else
        mct = 0
        th = datan(1.0d0)/45.0
100     if(ior.lt.0) write(*,3001)
        call dinput(td,6)
        if(errck) go to 100
        ni  = td(1)
        ne  = td(2)
        inc = td(3)
        x0  = td(4)
        y0  = td(5)
        z0  = td(6)
        if(ni.le.0) return
        if(ni.gt.numnp.or.ne.gt.numnp) go to 300
        inc = isign(max0(iabs(inc),1),ne-ni)
        if(ne.eq.0) ne = ni
        n = ni
200     r = x(1,n)
        snp    = sin(x(3,n)*th)
        x(1,n) = x0 + r*cos(x(2,n)*th)*snp
        x(2,n) = y0 + r*sin(x(2,n)*th)*snp
        x(3,n) = z0 + r*cos(x(3,n)*th)
        if(mct.gt.0) go to 250
        if(prt) then
cww                    write(iow,2000) o,head,x0,y0,z0,(i,i=1,ndm)
cww       if(ior.lt.0) write(*  ,2000) o,head,x0,y0,z0,(i,i=1,ndm)
                       write(iow,2000)        x0,y0,z0,(i,i=1,ndm)
          if(ior.lt.0) write(*  ,2000)        x0,y0,z0,(i,i=1,ndm)
        end if
        mct = 50
250     if(prt) then
          write(iow,2001) n,(x(i,n),i=1,ndm)
          if(ior.lt.0) write(*,2001) n,(x(i,n),i=1,ndm)
        end if
        mct = mct - 1
        n = n + inc
        if((ne-n)*inc.ge.0) go to 200
        if(mod(ne-ni,inc).eq.0) go to 100
        ni = ne
        n = ne
        go to 200
      end if
c.... error
cww300   write(iow,3000) ni,ne
cww      if(ior.lt.0) write(*,3000) ni,ne
cww      stop
300   write(yyy,3000) ni,ne
      call drawmess(yyy,1,0)
      return
c.... formats
cww2000  format(a1,19a4,a3//
2000  format(/
     1 '   cartesian coordinates computed from spherical input.'
     2 /5x,'x0 = ',e12.4,'    y0 = ',e12.4,'    z0 = ',e12.4/
     3   4x,'node',6(i6,'-coord')/(8x,6(i6,'-coord')))
2001  format(i8,6f12.4/(8x,6f12.4))
3000  format('  **ERROR** attempt to convert nodes ni = ',i6,
     1 ' to ne = ',i6)
3001  format(' Input: node-1,node-2,inc, x0, y0, z0'/'   >',$)
4000  format(' **ERROR** attempt to convert spherical coordinates'/
     1       '             for a problem with less than 3-dimensions.')
      end
c
      subroutine storeh(h,s,nst,ni)
c----------------------------------------------------------------------
c
c      Purpose: Saves static condensation part of element tangent
c               in history (enhanced strain or mixed type elements)
c
c      Inputs:
c         s(nst,*)  - Array to save in history
c         nst       - Dimension of s array
c         ni        - Dimension of h array
c
c      Outputs:
c         h(ni,*)   - Saved history array
c
c----------------------------------------------------------------------
      double precision h(ni,nst),s(nst,nst)
      do 110 i = 1,ni
         do 100 j = 1,nst
            h(i,j) = s(i,j)
100      continue
110   continue
      return
      end
c
      subroutine summariz(ct,b,badd,str,stradd,nneq,npp,fac)
c----------------------------------------------------------------------
c
c      Purpose: Summarize nodal displacement and stress values from disk
c          b(i)   = b(i)   + fac*badd(i)
c          str(i) = str(i) + fac*stradd(i)
c
c      Inputs:
c         ct        - Name of array, file rewind, or file name
c         b(*)      - Nodal displacement values
c         badd(*)   - Nodal displacement array to add
c         str(*)    - Nodal stress values
c         stradd(*) - Nodal stress array to add
c         nneq      - Length of displacement array    nneq = numnp*ndf
c         npp       - Length of stress       array    npp  = numnp*npstr
c         fac       - Factor for adding arrays
c
c      Outputs:
c         b(*)      - Nodal displacement values
c         str(*)    - Nodal stress values
c
c----------------------------------------------------------------------
      USE iodata
      USE iofile
      USE strnam
      logical lflg
      character*4 ct,cc,fname1,yyy*80
      real *8 b(*),badd(*),str(*),stradd(*),fac
      integer*4 nneq,npp
      save lflg
      data lflg /.false./
c.... go to processor
      if(ct.eq.'zero') go to 500
c.... open file set with name 'ct'
      fname1 = ct
      inquire(file=fname1,exist=lflg)
      if(lflg) then
                       write(iow,2001) fname1
          if(ior.lt.0) write(*  ,2001) fname1
          open(ird,file=fname1,status='old',form='unformatted',err=990)
c....   read displacement state from disk
          read(ird,end=980,err=990) cc
          if(cc.ne.'disp') go to 970
          read(ird,end=980,err=990) (badd(i),i=1,nneq)
                       write(iow,2100) fname1,fac
          if(ior.lt.0) write(*  ,2100) fname1,fac
        do i=1,nneq
          b(i) = b(i) + fac*badd(i)
        end do
c....   read nodal stress state from disk
          read(ird,end=980,err=990) cc
          if(cc.ne.'stre') go to 970
          read(ird,end=980,err=990) istv,(stradd(i),i=1,npp)
                       write(iow,2200) fname1,fac
          if(ior.lt.0) write(*  ,2200) fname1,fac
        do i=1,npp
          str(i) = str(i) + fac*stradd(i)
        end do
        close(ird)
        lflg = .false.
      else
          write(yyy,2060) fname1
          call drawmess(yyy,1,0)
      end if
      return
c.... clear all fields
500   call pzero(b,nneq)
        call pzero(str,npp)
                       write(iow,2500)
          if(ior.lt.0) write(*  ,2500)
      return
970   write(yyy,2070) ct,cc
        call drawmess(yyy,1,0)
      return
980   write(yyy,2080) ct
        call drawmess(yyy,1,0)
      return
990   write(yyy,2090) ct
        call drawmess(yyy,1,0)
      return
2001  format('File name for a read has been set to ',a4)
2060  format('File ',a4,' does not exist')
2070  format('SUMM requested ',a4,' but found ',a4)
2080  format('End of file on a summ command for ',a4)
2090  format('on a summ command for ',a4)
2100  format('Add DISP from file ',a4,' with factor ',f10.5)
2200  format('Add STRE from file ',a4,' with factor ',f10.5)
2500  format('Clear fields for DISP and STRE')
      end
c
c----------------------------------------------------------------------
c
      subroutine tgauss(l,lint,el,wg)
c----------------------------------------------------------------------
c
c      Purpose: triangle gauss point integration
c
c      Inputs:
c         l       - Number of points: 1,3,4,7
c
c      Outputs:
c         lint     - Total number of points
c         el(3,7)  - 1-3 coordinates for each Gauss point
c         wg(7)    - Gauss weight
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension el(3,7),wg(7)
c.... 1-point gauss integration
      if(l.eq.1) then
        el(1,1) = 1.d0/3.d0
        el(2,1) = el(1,1)
        el(3,1) = el(1,1)
        wg(1)   = 1.d0
        lint    = 1
c.... 3-point gauss integration
      else if(l.eq.3) then
        el(1,1) = 0.d0
        el(2,1) = 0.5d0
        el(3,1) = 0.5d0
        el(1,2) = 0.5d0
        el(2,2) = 0.d0
        el(3,2) = 0.5d0
        el(1,3) = 0.5d0
        el(2,3) = 0.5d0
        el(3,3) = 0.d0
        wg(1)   = 1.d0/3.d0
        wg(2)   = wg(1)
        wg(3)   = wg(1)
        lint    = 3
c.... 4-point gauss integration
      else if(l.eq.4) then
        el(1,1) =  1.d0/3.d0
        el(2,1) =  el(1,1)
        el(3,1) =  el(1,1)
        el(1,2) =  0.6
        el(2,2) =  0.2
        el(3,2) =  el(2,2)
        el(1,3) =  el(2,2)
        el(2,3) =  el(1,2)
        el(3,3) =  el(2,2)
        el(1,4) =  el(2,2)
        el(2,4) =  el(2,2)
        el(3,4) =  el(1,2)
        wg(1)   = -27.d0/48.d0
        wg(2)   =  25.d0/48.d0
        wg(3)   =  wg(2)
        wg(4)   =  wg(2)
        lint    =  4
c.... 7-point gauss integration
      else if(l.eq.7) then
        r0      =  sqrt(15.0d0)
        r1      =  3.d0/7.d0
        r2      =  (r0 + r0)/21.d0
        el(1,1) =  1.d0/3.d0
        el(2,1) =  el(1,1)
        el(3,1) =  el(1,1)
        el(1,2) =  r1 + r2
        el(2,2) =  0.5d0 - 0.5d0*el(1,2)
        el(3,2) =  el(2,2)
        el(1,3) =  el(2,2)
        el(2,3) =  el(1,2)
        el(3,3) =  el(2,2)
        el(1,4) =  el(2,2)
        el(2,4) =  el(2,2)
        el(3,4) =  el(1,2)
        el(1,5) =  r1 - r2
        el(2,5) =  0.5d0 - 0.5d0*el(1,5)
        el(3,5) =  el(2,5)
        el(1,6) =  el(2,5)
        el(2,6) =  el(1,5)
        el(3,6) =  el(2,5)
        el(1,7) =  el(2,5)
        el(2,7) =  el(2,5)
        el(3,7) =  el(1,5)
        wg(1)   =  0.225d0
        wg(2)   =  (155.d0 - r0)/1200.d0
        wg(3)   =  wg(2)
        wg(4)   =  wg(2)
        wg(5)   =  (155.d0 + r0)/1200.d0
        wg(6)   =  wg(5)
        wg(7)   =  wg(5)
        lint    =  7
      else
        write(*,2000) l
        lint    = -1
      end if
      return
2000  format(' * * ERROR * * call with wrong quadrature, l =',i3)
      end
c
      subroutine tienod(id,ix,x,f,ip,ndm,ndf,nen,nen1,
     1                  numnp,numel,ntied,n1,n2,n3,n4,nt,td,prt)
c----------------------------------------------------------------------
c
c      Purpose: Procedure to connect nodes which have same coordinates.

c      Inputs:
c         id(ndf,*)   - Boundary conditions +1/-1=fixed 0=free
c     [later: id(ndf,*)   - Equation numbers for each active dof]
c         ix(nen1,*)  - Element nodal connection list
c         x(ndm,*)    - Nodal coordinates for mesh
c         f(ndf,*)    - Nodal force values
c         ip(*)       - Node number list for tied nodes
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         nen         - Number of nodes/element
c         nen1        - Dimension of ix array
c         numnp       - Number of nodes in mesh
c         numel       - Number of elements in mesh
c         n1          - First  node number for search (Node option)
c         n2          - Second node number for search (Node option)
c         n3          - First  region number for tie (Region option)
c         n4          - Second region number for tie (Region option)
c         nt          - type of calculation
c         tol=td(1)   - Tolerance for tie
c
c         Node - Option nt=1 -----------------------
c         n1             - node 1
c         n2             - node 2
c
c         Line - Option nt=2 -----------------------
c         x1 = td(2)     - x-Coordinates of point 1
c         y1 = td(3)     - y-Coordinates of point 1
c        <z1 = td(4)     - z-Coordinates of point 1>
c         x2 = td(ndm+2) - x-Coordinates of point 2
c         y2 = td(ndm+3) - y-Coordinates of point 2
c        <z2 = td(ndm+4) - z-Coordinates of point 2>
c
c         Mlin - Option nt=3 -----------------------
c         like Line
c         +
c         n3             - mate 1
c         n4             - mate 2
c
c
c         Dir - Option nt=4 ------------------------
c         mt = td(2)     - Direction
c         xt = td(3)     - Coordinate value
c
c         Mate - Option nt=5 -----------------------
c         n3             - mate 1
c         n4             - mate 2


c         prt         - Print Flag
c
c      Outputs:
c         ix(nen1,*)  - Element nodal connection list with tied node
c                       eliminated
c         x(ndm,*)    - Nodal coordinates, tied nodes x1=-999
c         ntied       - Number of tied nodes: numnp_tied = numnp-ntied
c
c
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      logical fflg,prt
      integer id(ndf,*),ix(nen1,numel),ip(numnp),il(10)
      real*8  x(ndm,numnp),f(ndf,numnp),tol(3),td(7)
      real*8  p(3),p1(3),p2(3),dp(3),a(3)
      data bl/-999.d0/

c     Set tolerance for each direction to tie nodes
      if(td(1).gt.0.0d0) then
        sn = td(1)
      else
        sni = numnp
        sni = sqrt(sni)*1000.d0
        sn= 1.d0/sni
      endif
      if(prt) write(iow,2004) sn
c
      do 110 i = 1,ndm
          xmx = x(i,1)
          xmn = x(i,1)
          do 100 n = 1,numnp
            if(x(1,n).eq.bl) go to 100
            xmx = max(xmx,x(i,n))
            xmn = min(xmn,x(i,n))
100       continue
110   tol(i) = (xmx - xmn)*sn

c     Set circular tolerance for nodes on line to tie
      if(ndm.eq.2) tol(3) = 0.d0
      tolc = dsqrt(dot(tol,tol,3))

cww   only once! for repetition of macro tie! this is done in pcontr with pconsi
cwwc.... set up original numbers
cww      do n = 1,numnp
cww        ip(n) = n
cww      end do


c.... Set line values for tie,line and tie,mlin
      if(nt.eq.2.or.nt.eq.3) then
c....   Point 1
        p1(1) = td(2)
        p1(2) = td(3)
        p1(3) = 0.d0
        if(ndm.eq.3) p1(3) = td(4)
c....   Point 2
        p2(1) = td(ndm+2)
        p2(2) = td(ndm+3)
        p2(3) = 0.d0
        if(ndm.eq.3) p2(3) = td(ndm+4)
c....   vector a = p2-p1, aa = (a*a)
        a(1) = p2(1) - p1(1)
        a(2) = p2(2) - p1(2)
        a(3) = p2(3) - p1(3)
        aa   = dot(a,a,3)
      end if

c.... Set search coordinate for tie,dir
      nr = 0
      xt = 0.d0
      if(nt.eq.4) then
        nr = int(td(2))
        nr = min(ndm,max(1,nr))
        xt = td(3)
      end if

c.... do an n-squared search to tie
      if(prt) write(iow,2000)
      if(prt) write(  *,2000)

      call perform(0,numnp,1,31)

c.... loop node k
cww   do 250 k = 1,numnp - 1
      do 250 k = n1,n2-1   ! only between node n1 and n2

c....   node has coordinate x1=-999.d0
        if(x(1,k).eq.bl) go to 250

c....   tie,dir: check if node is in direction
      if(nr.gt.0) then
        if(abs(x(nr,k)-xt).gt.tol(nr)) go to 250
      end if

c....   tie,line, tie,mlin: check if node is on line
        if(nt.eq.2.or.nt.eq.3) then
          call pzero(dp,3)
          call pzero(p,3)
          do i = 1,ndm
            p(i)  = x(i,k)
            dp(i) = p(i) - p1(i)
          end do
c....     test if node P_k is near line dp=p-p1 =(!) dlamb*a, a = p2-p1
          dlamb = dot(dp,a,3)/aa
          if(dlamb.lt.-tolc.or.dlamb.gt.(1.d0+tolc)) goto 250 ! 0<dlamb<1
c.....    vector dp = dp - dlamb`*a
          do i= 1,3
            dp(i) = dp(i) - dlamb*a(i)
          end do
c.....    length of d
          dd = dsqrt(dot(dp,dp,3))

          if(dd.gt.tolc) goto 250

c         write(*,*) 'node on line',k,dd

        end if

c....   loop node j
cww     do 240 j = k+1,numnp
        do 240 j = k+1,n2   ! only between node k+1 and n2

c....     node has coordinate x1=-999.d0
          if(x(1,j).eq.bl) go to 240

c....     distance between k and j in each dir
          do 200 i = 1,ndm
            if(abs(x(i,k)-x(i,j)).gt.tol(i)) go to 240
200       continue

c....     look for mate option mate1=n3, mate2=n4
          if(nt.eq.3.or.nt.eq.5) then
            jn3 = -1   ! set to a never used value
            jn4 = -1
            kn3 = -1
            kn4 = -1
            do jj = 1,numel  ! find mate of node j
              do i = 1,nen
                if(ix(i,jj).eq.j) then
                  if(ix(nen1,jj).eq.n3) jn3 = n3
                  if(ix(nen1,jj).eq.n4) jn4 = n4
                  if(jn3.eq.n3 .or.jn4.eq.n4) goto 201 ! node found
                  jn3 = -1
                  jn4 = -1
                end if
              end do
            end do

201         do kk = 1,numel  ! find mate of node k
              do i = 1,nen
                if(ix(i,kk).eq.k) then
                  if(ix(nen1,kk).eq.n3) kn3 = n3
                  if(ix(nen1,kk).eq.n4) kn4 = n4
                  if(kn3.eq.n3 .or.kn4.eq.n4) goto 202 ! node found
                  kn3 = -1
                  kn4 = -1
                end if
              end do
            end do

202         if(jn3.eq.n3 .and.kn4.eq.n4) goto 203  ! nodes to tie
            if(jn4.eq.n4 .and.kn3.eq.n3) goto 203  ! nodes to tie
            goto 240
          end if

c....     connect k-j
c....     connect node-j to node-k in the id-list (eliminate node-j)
203       do i = 1,numnp
            if(ip(i).eq.ip(j)) then
              ip(i) = ip(k)
            end if
          end do
c....     look at each dof - fix if either is restrained
          do i = 1,ndf
            if(id(i,j).ne.0 .or. id(i,k).ne.0) then
              id(i,j) = 1
              id(i,k) = 1
            end if
          end do
c....     set sum of loads on both nodes
          do i = 1,ndf
            if(id(i,j).eq.0 .and. id(i,k).eq.0) then ! no b.c.s
              fjk = f(i,j)+f(i,k)
              f(i,j) = fjk     ! set sum of load on node j and k
              f(i,k) = fjk
            end if
          end do

c     old version WW+FG+SLAU 10-2013
cc....     check if a force is applied to both of the tied nodes
c          do i = 1,ndf
c            if(f(i,j).ne.f(i,k)) then
c              if(f(i,j).ne.0.0d0 .and. f(i,k).eq.0.0d0) then
c                f(i,k) = f(i,j)  ! set load on node k
c              else if(f(i,j).eq.0.0d0 .and. f(i,k).ne.0.0d0) then
c                f(i,j) = f(i,k)  ! set load on node i
c              else                ! different loads
c                        write(iow,3000) j,k
c                if(prt) write(*  ,3000) j,k
ccww             if(ior.lt.0.and.prt) write(*  ,3000) j,k
ccww...added
c                fik = f(i,j)+f(i,k)
c                f(i,j) = fik     ! set sum of load on node i and k
c                f(i,k) = fik
c              end if
c            end if
c          end do
c3000  format(5x,' **WARNING** Inconsistent Loads applied to tied',
c     1     ' nodes',i5,' and',i5,/,8x,'Loads are added!' ,/1x)
c

240     continue ! loop j

        call perform(k,numnp,2,31)

250   continue   ! loop k

      call perform(numnp,numnp,3,31)

c.... eliminate node j from element connections for solution
      do 270 n = 1,numel
        do 260 i = 1,nen
          j = abs(ix(i,n))
          if(j.gt.0) ix(i,n) = ip(j)
260     continue
270   continue

c.... output the list of nodes removed
      if(prt) write(  *,2001)
      if(prt) write(iow,2001)
      i1 = 0
      fflg = .true.

      do 290 n = 1,numnp
        if(ip(n).ne.n) then
c....     delete equations for all unused nodes from a tie
          x (1,n) = bl
          do 280 i = 1,ndf
            id(i,n) = 1
            f (i,n) = 0.0
280       continue
          i1 = i1 + 1
          il(i1) = n
        end if
        if(i1.ge.10 .or. (i1.gt.0 .and. n.eq.numnp)) then
          if(prt) write(  *,2002) (il(i),i=1,i1)
          if(prt) write(iow,2002) (il(i),i=1,i1)
          i1 = 0
          fflg = .false.
        end if
290   continue

      if(fflg) then
        if(prt) write(*  ,2003)
        if(prt) write(iow,2003)
      end if

c.... count number of tied nodes
      ntied = 0
      do n = 1, numnp
        if(x(1,n).eq.bl) ntied=ntied+1
      end do

      return
2000  format(/6x,' t i e  - - -  n o d a l   c o o r d i n a t e s'/
     1    6x,'    Nodes with same coordinates will be connected.'/)
2001  format(6x,'    The following nodes have been deleted:'/)
2002  format(6x,10i5)
2003  format(10x,'No nodes have been removed by the tie command'/1x)
2004  format(/,6x,'Tie: Gap Tolerance =',1p,1e12.5/1x)
      end
c
      subroutine transx(x,ndm,numnp,prt)
c----------------------------------------------------------------------
c
c      Purpose: move cartesian coordinates  Xnew = Xold - dx
c
c      Input:
c         x(ndm,*)    - Nodal coordinates of mesh
c         ndm         - Spatial dimension of mesh
c         numnp       - Number of nodes in mesh
c
c      Output:
c         x(ndm,*)    - Nodal coordinates of mesh
c
c----------------------------------------------------------------------
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension x(ndm,*),td(6)
      character yyy*80
      mct = 0
100   if(ior.lt.0) write(*,3001)
      call dinput(td,6)
      if(errck) go to 100
      ni  = td(1)
      ne  = td(2)
      inc = td(3)
      dx  = td(4)
      dy  = td(5)
      dz  = td(6)
      if(ni.le.0) return
      if(ni.gt.numnp.or.ne.gt.numnp) go to 300
      inc = isign(max(iabs(inc),1),ne-ni)
      if(ne.eq.0) ne = ni
      n = ni
c...  calculate new values
200     x(1,n) = x(1,n) - dx
        x(2,n) = x(2,n) - dy
        if(ndm.eq.3) x(3,n) = x(3,n) - dz
      if(mct.gt.0) go to 250
c
      if(prt)              write(iow,2000) dx,dy,dz,(i,i=1,ndm)
      if(prt.and.ior.lt.0) write(*,  2000) dx,dy,dz,(i,i=1,ndm)
      mct = 50
250   if(prt)              write(iow,2001) n,(x(i,n),i=1,ndm)
      if(prt.and.ior.lt.0) write(*,  2001) n,(x(i,n),i=1,ndm)
      mct = mct - 1
      n = n + inc
      if((ne-n)*inc.ge.0) go to 200
      if(mod(ne-ni,inc).eq.0) go to 100
      ni = ne
      n  = ne
      go to 200
c.... error
300   write(yyy,3000) ni,ne
      call drawmess(yyy,1,0)
      return
c.... formats
2000  format(/,
     1 '  New cartesian coordinates computed from TRANS',/,
     2 '  with dx = ',f12.4,' dy = ',f12.4,' dz = ',f12.4,/,
     3 4x,'node',3(i6,'-coord'))
2001  format(i8,3f12.4)
3000  format(
     +' Try to use nodes in TRANS from',i4,' to ',i4,'(> max node)')
3001  format(' Input: node-1,node-2,inc, dx, dy, (dz)'/'   >',$)
      end
c
      subroutine trblk(nr,xl,ixl,x,ix,ndm,nod1,nuel1,
     1                nel1,nm,ma,prt)
c----------------------------------------------------------------------
c
c      Purpose: Generate a triangular block of 3-node triangular elements
c
c      Inputs:
c         nr        - Number elements in 1-local coordinate dir.
c         xl(ndm,*) - Block nodal coordinate array
c         ixl(*)    - Block nodal connection list
c         ndm       - Spatial dimension of mesh
c         nod1      - Initial node number for block
c         nuel1     - Initial element number for block
c         nel1      - Dimension of ix array
c         nm
c         ma        - Material set number for block
c         prt       - Output generated data if true
c
c      Outputs:
c         x(ndm,*)  - Nodal coordinates for block
c         ix(*)     - Element nodal connection list for block
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE iofile
c..... Declare variable types
      logical prt
      character*6 xh
      integer ni,nn,mct,k,i,j,i1,j1,nei,n1,n2,
     1        nr,ndm,nod1,nuel1,nel1,nm,ma
      real*8  dl
c..... Declare array types
      integer ixl(*),ix(nel1,*)
      real*8  xl(3,*),x(ndm,*),tshp(6), el(3)
c
      data xh/' coord'/
      ni  = nr
      mct = 0
c
      do 50 i=1,3
        j = mod(i,3) + 1
        if((ixl(i+3)).ne.0) then
          do 40 i1=1,3
            xl(i1,i+3) = xl(i1,i+3) - 0.5d0*(xl(i1,i) + xl(i1,j))
   40     continue
        end if
   50 continue

c.... generate nodes
      nn = nod1
      dl  = 1.0d0/ni
      do 130 i=0,ni
        el(3) = i*dl
c
        do 120 j=0,ni-i
          el(2) = j*dl
          el(1) = 1.0d0 - el(3) - el(2)
c.... form the shape functions
          do 105 i1=1,3
            j1 = mod(i1,3) + 1
            tshp(i1)   = el(i1)
            tshp(i1+3) = 4.0d0*el(i1)*el(j1)
  105     continue
c
cww       if(nn.gt.numnp) then
cww         write(*,*) ' trying to generate node',nn
cww         stop
cww       end if
          if(nn.gt.numnp) goto 131
          do 110 i1 = 1,ndm
            x(i1,nn) = 0.0d0
            do 100 k=1,6
              x(i1,nn) = x(i1,nn) + tshp(k)*xl(i1,k)
  100       continue
  110     continue
          if(prt) then
            mct = mct + 1
cww         if(mod(mct,50).eq.1) write(iow,2003) o,head,(k,xh,k=1,ndm)
            if(mod(mct,50).eq.1) write(iow,2003) (k,xh,k=1,ndm)
            write(iow,2004) nn,(x(k,nn),k=1,ndm)
            if(ior.lt.0) then
cww           if(mod(mct,50).eq.1) write(*,2003) o,head,(k,xh,k=1,ndm)
              if(mod(mct,50).eq.1) write(*,2003) (k,xh,k=1,ndm)
              write(*,2004) nn,(x(k,nn),k=1,ndm)
            end if
          end if
          nn = nn + 1
  120   continue
  130 continue
  131 continue
c
c.... generate elements
      nei = nuel1
      n1 = nod1
      do 220 i=1,ni
        do 210 j=1,ni-i+1
          n2 = n1 + ni - i + 2
c
          ix(1,nei)    = n1 + j - 1
          ix(2,nei)    = n1 + j
          ix(3,nei)    = n2 + j - 1
          ix(nel1,nei) = ma
          if(j.lt.(ni-i+1)) then
            nei = nei + 1
          if(nei.gt.numel) goto 221
            ix(1,nei)    = n2 + j - 1
            ix(2,nei)    = n1 + j
            ix(3,nei)    = n2 + j
            ix(nel1,nei) = ma
          end if
          nei = nei + 1
        if(nei.gt.numel) goto 221
  210   continue
        n1 = n2
  220 continue
  221 continue
c
cww2003  format(a1,19a4,a3/'  n o d a l   c o o r d i n a t e s'/
cww     1         6x,'node',5(i7,a6))
2003  format(/'  n o d a l   c o o r d i n a t e s'/
     1    6x,'node',5(i7,a6))
2004  format(i10,5f13.4)
      end
c
c-----------------------------------------------------------------------
c
      subroutine umacr1(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jp,u,dr,
     1           lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,plo)
c-----------------------------------------------------------------------
c
c      Purpose: User defined macro 1
c
c....  Inputs:
c       ul(nst,6)      - element displacements
c       xl(ndm,nen)    - element coordinates
c       tl(nen)        - element temperatures
c       ld(nst)        -
c       p(nst)         - element load vector
c       s(nst,nst)     - element stiffness matrix
c       ie(nie,numat)  - assembly information for material set nie=ndf+2
c       d(ndd,numat)   - material set parameters
c       id(ndf,numnp)  - equation numbers for each active dof
c       x(ndm,numnp)   - nodal coordinates of mesh
c       ix(nen1,numel) - element nodal connections of mesh
c       f0(ndf,numnp)  - nodal initial force values
c       f (ndf,numnp)  - load vector
c       t(numnp)       - temperature vector
c       jp(*)          - pointer array for row/columns of tangent
c       u(ndf,3*numnp) - displacement vector
c       dr(nneq)       - working array
c       lct            - second macro
c       ct(3,200)      - parameter array for 200 macros
c       ndf            - number dof/node
c       ndm            - spatial dimension of mesh
c       nen            - max. number of nodes/element
c       nen1           - dimension for ix array: nen+4
c       nst            - dimension for element array: ndf*nen
c       nneq           - number of equations numnp*ndf
c       prt            - print option
c       plo(7,nplo)    - array for plotting TPLO-data
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      call drawmess('**WARNING** call to dummy procedure UMACR1 *',1,0)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine umacr2(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jp,u,dr,
     1           lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,plo)
c-----------------------------------------------------------------------
c
c      Purpose: User defined macro 2
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      call drawmess('**WARNING** call to dummy procedure UMACR2 *',1,0)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine umacr3(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jp,u,dr,
     1           lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,plo)
c-----------------------------------------------------------------------
c
c      Purpose: User defined macro 3
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      call drawmess('**WARNING** call to dummy procedure UMACR3 *',1,0)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine umacr4(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jp,u,dr,
     1           lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,plo)
c-----------------------------------------------------------------------
c
c      Purpose: User defined macro 4
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      call drawmess('**WARNING** call to dummy procedure UMACR4 *',1,0)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine umacr5(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jp,u,dr,
     1           lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,plo)
c-----------------------------------------------------------------------
c
c      Purpose: User defined macro 5
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      call drawmess('**WARNING** call to dummy procedure UMACR5 *',1,0)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine umesh1(idl,ie,d,id,x,ix,f,t,ndd,nie,ndf,ndm,nen1,iii,
     1                  prt)
c----------------------------------------------------------------------
c
c      Purpose: user Data input routine UMESH 1 for mesh description
c
c      Inputs:
c         idl(nst)    - scratch array
c         ie(nie,*)   - Assembly information for material set
c         d(ndd,*)    - Material set parameters
c         id(ndf,*)   - Equation numbers for each active dof
c         x(ndm,*)    - Nodal coordinates of mesh
c         ix(nen1,*)  - Element nodal connections of mesh
c         f(ndf,*,2)  - Nodal force and displacement values
c         t(*)        - Nodal temperature values
c         ndd         - Dimension for d array
c         nie         - Dimension for ie array
c         ndf         - Number dof/node
c         ndm         - Spatial dimension of mesh
c         nen1        - Dimension for ix array
c         iii         - Initialization indicator
c         prt         - Flag, print input data if true
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      integer idl(ndf*ndf)
      integer ie(nie,*)
      real*8  d(ndd,*)
      integer id(ndf,*)
      real*8  x(ndm,*)
      integer ix(nen1,*)
      real*8  f(ndf,*)
      real*8  t(*)
      integer ndd,nie,ndf,ndm,nen1,iii
      logical prt

      call drawmess('**WARNING**call to dummy input macro UMESH1 *',1,0)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine umesh2(idl,ie,d,id,x,ix,f,t,ndd,nie,ndf,ndm,nen1,iii,
     1                  prt)
c----------------------------------------------------------------------
c
c      Purpose: user Data input routine UMESH 2 for mesh description
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      call drawmess('**WARNING**call to dummy input macro UMESH2 *',1,0)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine umesh3(idl,ie,d,id,x,ix,f,t,ndd,nie,ndf,ndm,nen1,iii,
     1                  prt)
c----------------------------------------------------------------------
c
c      Purpose: user Data input routine UMESH 3 for mesh description
c
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      call drawmess('**WARNING**call to dummy input macro UMESH3 *',1,0)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine umesh4(idl,ie,d,id,x,ix,f,t,ndd,nie,ndf,ndm,nen1,iii,
     1                  prt)
c----------------------------------------------------------------------
c
c      Purpose: user Data input routine UMESH 4 for mesh description
c
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      call drawmess('**WARNING**call to dummy input macro UMESH4 *',1,0)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine umesh5(idl,ie,d,id,x,ix,f,t,ndd,nie,ndf,ndm,nen1,iii,
     1                  prt)
c----------------------------------------------------------------------
c
c      Purpose: user Data input routine UMESH 5 for mesh description
c
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      call drawmess('**WARNING**call to dummy input macro UMESH5 *',1,0)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine vblk(nr,ns,nt,xl,x,ixl,ix,id,dr,ds,dt,
     1                ni,ne,ndm,ndf,nel1,ma,ntyp,prt)
c----------------------------------------------------------------------
c
c      Purpose: Generate a block of 3-d 8-node hexahedron elements
c      ntyp = 10   8 - node hexahedrons
c      ntyp = 11   4 - node tetrahedrons
c      ntyp = 12  20 - node hexahedrons (21 nodes!)
c      ntyp = 13  27 - node hexahedrons
c      ntyp = 14  20 - node hexahedrons
c      ntyp = 15  18 - node hexahedrons
c      ntyp = 19  64 - node hexahedrons
c
c      Inputs:
c         nr        - Number elements in 1-local coordinate dir.
c         ns        - Number elements in 2-local coordinate dir.
c         nt        - Number elements in 3-local coordinate dir.
c         xl(ndm,*) - Block nodal coordinate array
c         ixl(*)    - Block nodal connection list
c         id(ndf,*) - b.c. conditions
c         dr        - 1-local coordinate increment
c         ds        - 2-local coordinate increment
c         dt        - 3-local coordinate increment
c         ni        - Initial node number for block
c         ne        -
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         nen1      - Dimension for ix array
c         ma        - Material type of block
c         ntyp      - Type of block coordinates
c         prt       - Output generated data if true
c
c         Outputs:
c         x(ndm,*)  - Nodal coordinates for block
c         ix(nen1,*)  - Element nodal connections of mesh
c
c      Comments:
c         15.05.15 WW IBS KIT set b.c. for unused nodes for ntyp=12/14
c
c-------------------------------------------------------------
c..... type declaration for variables
      USE bdata
      USE cdata
      USE iofile
      logical prt,phd
      character xh*6
      integer me,ni,ne,ndm,nel1,ma
      integer nr,ns,nt,nrs, i, j,jj, k, l, n, mct,ntyp,ixb(7)
      real*8  dr,ds,dt
c..... type declaration for arrays
      integer ixl(27),ix(nel1,*),id(ndf,*)
      real*8  ss(3),xl(3,*),x(ndm,*),shp(27)
      data xh/' coord'/
c.... check that all corners of brick are defined
      do 10 k = 1,8
        if(ixl(k).ne.k) go to 900
10    continue
      n = ni
      mct = 0
      ss(3) = -1.0d0
      do 300 k = 1,nt
        ss(2) = -1.0d0
        do 200 j = 1,ns
          ss(1) = -1.0d0
          do 100 i = 1,nr
c....       compute shape functions and coordinates for each point
            call gshape3(ss, ixl, shp)
            do 40 jj = 1,3
              x(jj,n) = 0.0d0
40          continue
            do 60 l = 1,27
              if(ixl(l).ne.0) then
                do 50 jj = 1,3
                  x(jj,n) = x(jj,n) + shp(l)*xl(jj,l)
50              continue
              end if
60          continue
c....       here modification of x possible
c           x = x + dx
c....       e.g.parabola  x_3 = a * x_1**2
c           x(3,n) = x(3,n) + a * x(1,n)**2
c....       output the point
            if(prt) then
              mct = mct + 1
              phd = mod(mct,50).eq.1
cww           if(phd) write(iow,2000) o,head,(l,xh,l=1,ndm)
              if(phd) write(iow,2000) (l,xh,l=1,ndm)
              write(iow,2001) n,(x(l,n),l=1,ndm)
              if(ior.lt.0) then
cww             if(phd) write(*,2000) o,head,(l,xh,l=1,ndm)
                if(phd) write(*,2000) (l,xh,l=1,ndm)
                write(*,2001) n,(x(l,n),l=1,ndm)
              end if
            end if
            n = n + 1
            if(n.gt.numnp) goto 301
            ss(1) = ss(1) + dr
100       continue
          ss(2) = ss(2) + ds
200     continue
        ss(3) = ss(3) + dt
300   continue
301   continue
c.... compute the element connections
      if(ne.le.0) return
      if(ntyp.eq.10 ) then
c....  for bricks  (ntyp = 10)
       me = ne - 1
       nrs = nr*ns
         do 600 k = 1,nt-1
          do 500 j = 1,ns-1
            n = nr*(j-1 + ns*(k-1)) + ni
            do 400 i = 1,nr-1
              n = n + 1
              me = me + 1
            if(me.gt.numel) goto 601
              ix(nel1,me) = ma
              ix(1,me)    = n - 1
              ix(2,me)    = n
              ix(3,me)    = n + nr
              ix(4,me)    = n + nr  - 1
              ix(5,me)    = n + nrs - 1
              ix(6,me)    = n + nrs
              ix(7,me)    = n + nrs + nr
              ix(8,me)    = n + nrs + nr - 1
400         continue
500       continue
600      continue
601    continue
      else if(ntyp.eq.11 ) then
c....  for tetrahedrons (ntyp = 11)
       me = ne - 1
       if(me.gt.numel) goto 6011
         nrs = nr*ns
         do 611 k = 1,nt-1
          do 511 j = 1,ns-1
            n = nr*(j-1 + ns*(k-1)) + ni
            do 411 i = 1,nr-1
              n = n + 1
              me = me + 1
            if(me.gt.numel) goto 6011
              ix(nel1,me) = ma
              ix(1,me)    = n
              ix(2,me)    = n + nr
              ix(3,me)    = n - 1
              ix(4,me)    = n + nrs - 1
              me = me + 1
            if(me.gt.numel) goto 6011
              ix(nel1,me) = ma
              ix(1,me)    = n + nrs - 1
              ix(2,me)    = n
              ix(3,me)    = n + nr
              ix(4,me)    = n + nrs + nr
              me = me + 1
            if(me.gt.numel) goto 6011
              ix(nel1,me) = ma
              ix(1,me)    = n + nrs + nr
              ix(2,me)    = n + nrs - 1
              ix(3,me)    = n
              ix(4,me)    = n + nrs
              me = me + 1
            if(me.gt.numel) goto 6011
              ix(nel1,me) = ma
              ix(1,me)    = n - 1
              ix(2,me)    = n + nr
              ix(3,me)    = n + nr  - 1
              ix(4,me)    = n + nrs + nr - 1
              me = me + 1
            if(me.gt.numel) goto 6011
              ix(nel1,me) = ma
              ix(1,me)    = n + nr
              ix(2,me)    = n + nrs + nr - 1
              ix(3,me)    = n - 1
              ix(4,me)    = n + nrs - 1
              me = me + 1
            if(me.gt.numel) goto 6011
              ix(nel1,me) = ma
              ix(1,me)    = n + nrs + nr - 1
              ix(2,me)    = n + nrs - 1
              ix(3,me)    = n + nr
              ix(4,me)    = n + nrs + nr
411         continue
511       continue
611      continue
6011   continue
      else if(ntyp.eq.12.or.ntyp.eq.13) then
c....  for 20(21!)/27 node bricks
       me = ne - 1
       if(me.gt.numel) goto 6012
         nrs = nr*ns
         do 612 k = 1,nt-1,2
          do 512 j = 1,ns-1,2
            n = nr*(j-1 + ns*(k-1)) + ni
            do 412 i = 1,nr-1,2
              n = n + 1
              me = me + 1
            if(me.gt.numel) goto 6012
              ix(nel1,me) = ma
c....         bottom  t=-1
                             ix( 1,me)    = n - 1
                             ix(13,me)    = n
                             ix( 2,me)    = n + 1
c....         bottom  t= 0
                             ix(16,me)    = n - 1 + nr
              if(ntyp.eq.13) ix(17,me)    = n     + nr
              if(ntyp.eq.12) ixb(1)       = n     + nr
                             ix(14,me)    = n + 1 + nr
c....         bottom  t=+1
                             ix( 4,me)    = n - 1 + nr*2
                             ix(15,me)    = n     + nr*2
                             ix( 3,me)    = n + 1 + nr*2
c....         midside t=-1
                             ix( 9,me)    = n - 1        + nrs
              if(ntyp.eq.13) ix(23,me)    = n            + nrs
              if(ntyp.eq.12) ixb(2)       = n            + nrs
                             ix(10,me)    = n + 1        + nrs
c....         midside t= 0
              if(ntyp.eq.13) then
                             ix(26,me)    = n - 1 + nr   + nrs
                             ix(27,me)    = n     + nr   + nrs
                             ix(24,me)    = n + 1 + nr   + nrs
              else if(ntyp.eq.12) then
                             ixb(3)       = n - 1 + nr   + nrs
                             ixb(4)       = n     + nr   + nrs
                             ixb(5)       = n + 1 + nr   + nrs
              end if
c....         midside t=+1
                             ix(12,me)    = n - 1 + nr*2 + nrs
              if(ntyp.eq.13) ix(25,me)    = n     + nr*2 + nrs
              if(ntyp.eq.12) ixb(6)       = n     + nr*2 + nrs
                             ix(11,me)    = n + 1 + nr*2 + nrs
c....         top     t=-1
                             ix( 5,me)    = n - 1        + nrs*2
                             ix(18,me)    = n            + nrs*2
                             ix( 6,me)    = n + 1        + nrs*2
c....         top     t= 0
                             ix(21,me)    = n - 1 + nr   + nrs*2
              if(ntyp.eq.13) ix(22,me)    = n     + nr   + nrs*2
              if(ntyp.eq.12) ixb(7)       = n     + nr   + nrs*2
                             ix(19,me)    = n + 1 + nr   + nrs*2
c....         top     t=+1
                             ix( 8,me)    = n - 1 + nr*2 + nrs*2
                             ix(20,me)    = n     + nr*2 + nrs*2
                             ix( 7,me)    = n + 1 + nr*2 + nrs*2

c....         set b.c. for unused nodes 
              if(ntyp.eq.12) then
                do ibc = 1,7
                  do idf=1,ndf
                    id(idf,ixb(ibc)) = 1
                  end do 
                end do                   
              end if 
              n = n + 1
412         continue
512       continue
612     continue
6012    continue
      else if(ntyp.eq.14 ) then
c....  for 20 node bricks
       me = ne - 1
       if(me.gt.numel) goto 6014
         nrs = nr*ns
         do k = 1,nt-1,2
          do j = 1,ns-1,2
            n = nr*(j-1 + ns*(k-1)) + ni
            do i = 1,nr-1,2
              n = n + 1
              me = me + 1
            if(me.gt.numel) goto 6014
              ix(nel1,me) = ma
c....         bottom  t=-1
                             ix( 1,me)    = n - 1
                             ix(13,me)    = n
                             ix( 2,me)    = n + 1
c....         bottom  t= 0
                             ix(16,me)    = n - 1 + nr
                             ixb(1)       = n     + nr
                             ix(14,me)    = n + 1 + nr
c....         bottom  t=+1
                             ix( 4,me)    = n - 1 + nr*2
                             ix(15,me)    = n     + nr*2
                             ix( 3,me)    = n + 1 + nr*2
c....         midside t=-1
                             ix( 9,me)    = n - 1        + nrs
                             ixb(2)       = n            + nrs
                             ix(10,me)    = n + 1        + nrs
c....         midside t= 0
                             ixb(3)       = n - 1  + nr  + nrs
                             ixb(4)       = n      + nr  + nrs
                             ixb(5)       = n + 1  + nr  + nrs
c....         midside t=+1
                             ix(12,me)    = n - 1 + nr*2 + nrs
                             ixb(6)       = n     + nr*2 + nrs
                             ix(11,me)    = n + 1 + nr*2 + nrs
c....         top     t=-1
                             ix( 5,me)    = n - 1        + nrs*2
                             ix(17,me)    = n            + nrs*2
                             ix( 6,me)    = n + 1        + nrs*2
c....         top     t= 0
                             ix(20,me)    = n - 1 + nr   + nrs*2
                             ixb(7)       = n     + nr   + nrs*2
                             ix(18,me)    = n + 1 + nr   + nrs*2
c....         top     t=+1
                             ix( 8,me)    = n - 1 + nr*2 + nrs*2
                             ix(19,me)    = n     + nr*2 + nrs*2
                             ix( 7,me)    = n + 1 + nr*2 + nrs*2

c....         set b.c. for unused nodes 
              do ibc = 1,7
                do idf=1,ndf
                  id(idf,ixb(ibc)) = 1
                end do 
              end do                   
              n = n + 1
            end do
          end do
        end do
6014    continue
      else if(ntyp.eq.15) then
c....  for 18 node bricks
       me = ne - 1
       if(me.gt.numel) goto 6015
         nrs = nr*ns
          do 615 k = 1,nt-1,2
c2         do 615 k = 1,nt-1,1
c3         do 615 k = 1,nt
          do 515 j = 1,ns-1,2
            n = nr*(j-1 + ns*(k-1)) + ni
            do 415 i = 1,nr-1,2
              n = n + 1
              me = me + 1
            if(me.gt.numel) goto 6015
              ix(nel1,me) = ma
c....         bottom  t=-1
                             ix( 1,me)    = n - 1
                             ix( 9,me)    = n
                             ix( 2,me)    = n + 1
c....         bottom  t= 0
                             ix(12,me)    = n - 1 + nr
                             ix(13,me)    = n     + nr
                             ix(10,me)    = n + 1 + nr
c....         bottom  t=+1
                             ix( 4,me)    = n - 1 + nr*2
                             ix(11,me)    = n     + nr*2
                             ix( 3,me)    = n + 1 + nr*2
c....         top     t=-1
                             ix( 5,me)    = n - 1        + nrs
                             ix(14,me)    = n            + nrs
                             ix( 6,me)    = n + 1        + nrs
c....         top     t= 0
                             ix(17,me)    = n - 1 + nr   + nrs
                             ix(18,me)    = n     + nr   + nrs
                             ix(15,me)    = n + 1 + nr   + nrs
c....         top     t=+1
                             ix( 8,me)    = n - 1 + nr*2 + nrs
                             ix(16,me)    = n     + nr*2 + nrs
                             ix( 7,me)    = n + 1 + nr*2 + nrs
              n = n + 1
415         continue
515       continue
615     continue
6015    continue
      else if(ntyp.eq.19 ) then
c....   for 64 node hexahedrons
        nf = ne !??
        call pqrblk(3,3,3, nr,ns,nt, nf,ni, ma, ix,nel1) !nel1!
      end if
      return
c.... error
900   write(iow,3000) k
      if(ior.lt.0) then
        write(*,3000) k
          return
      end if
cww      stop
      return
c.... formats
cww2000  format(a1,19a4,a3/'  n o d a l   c o o r d i n a t e s'/
cww     1         6x,'node',5(i7,a6))
2000  format(/'  n o d a l   c o o r d i n a t e s'/
     1    6x,'node',5(i7,a6))
2001  format(i10,5f13.4)
3000  format(' **ERROR** Block node',i3,' is undefined')
      end
c
      subroutine gshape3(xi, ixl, shp)
c----------------------------------------------------------------------
c
c      Purpose: Shape functions for 3-d mesh generation by block.
c               8 to 27 nodes.  No derivatives computed.
c
c      Inputs:
c         xi(3)   - Natural coordinates for point
c         ixl(*)  - List of active nodes
c
c                   b = bottom
c                   t = top
c                   m = mid
c
c      Outputs:
c         shp(*)  - Shape functions for point
c
c
c      Bug:       - combination of nodes 9-12/13-16/18-21 with 23-26 leads to errors WW 27.04.06
c
c      Comment:   - modification for node 27, seems to be ok. WW 27.04.06
c----------------------------------------------------------------------
c..... type declaration for variables
      integer i, j
      real*8  f1,f2,f3
c..... type declaration for arrays
      real*8  xi(3),shp(27)
      integer nxi(3,27),ixl(27)
      data nxi/-1,-1,-1,   1,-1,-1,   1, 1,-1,  -1, 1,-1,  !  1- 4
     1         -1,-1, 1,   1,-1, 1,   1, 1, 1,  -1, 1, 1,  !  5- 8
     2         -1,-1, 0,   1,-1, 0,   1, 1, 0,  -1, 1, 0,  !  9-12
     3          0,-1,-1,   1, 0,-1,   0, 1,-1,  -1, 0,-1,  ! 13-16
     5          0, 0,-1,                                   ! 17
     6          0,-1, 1,   1, 0, 1,   0, 1, 1,  -1, 0, 1,  ! 18-21
     7          0, 0, 1,                                   ! 22
     8          0,-1, 0,   1, 0, 0,   0, 1, 0,  -1, 0, 0,  ! 23-26
     7          0, 0, 0/                                   ! 27

      call pzero(shp,27)  ! added ww

c.... generate first functions (all corner nodes b,t,m)
      do i = 1,4
        f1 = 0.25*(1.+nxi(1,i)*xi(1))*(1.+nxi(2,i)*xi(2))
        f2 = 1.-xi(3)
        f3 = 1.+xi(3)
        shp(i  ) = 0.50*f1*f2
        shp(i+4) = 0.50*f1*f3
        shp(i+8) =      f1*f2*f3
      end do

c.... form relative shape functions - level 1    midside nodes b,t
      f1      = 0.25*(1.-xi(1)**2)
      if(ixl(13).ne.0)
     1   shp(13) = f1*(1.+nxi(2,13)*xi(2))*(1.+nxi(3,13)*xi(3))
      if(ixl(15).ne.0)
     1   shp(15) = f1*(1.+nxi(2,15)*xi(2))*(1.+nxi(3,15)*xi(3))

      if(ixl(18).ne.0)
     1   shp(18) = f1*(1.+nxi(2,18)*xi(2))*(1.+nxi(3,18)*xi(3))
      if(ixl(20).ne.0)
     1   shp(20) = f1*(1.+nxi(2,20)*xi(2))*(1.+nxi(3,20)*xi(3))

      f1      = 0.25*(1.-xi(2)**2)
      if(ixl(14).ne.0)
     1   shp(14) = f1*(1.+nxi(1,14)*xi(1))*(1.+nxi(3,14)*xi(3))
      if(ixl(16).ne.0)
     1   shp(16) = f1*(1.+nxi(1,16)*xi(1))*(1.+nxi(3,16)*xi(3))

      if(ixl(19).ne.0)
     1   shp(19) = f1*(1.+nxi(1,19)*xi(1))*(1.+nxi(3,19)*xi(3))
      if(ixl(21).ne.0)
     1   shp(21) = f1*(1.+nxi(1,21)*xi(1))*(1.+nxi(3,21)*xi(3))

c.... form relative shape functions - level 2   midface nodes b,t
      f1      = 1.-xi(1)**2
      f2      = 1.-xi(2)**2
      f3      = 1.-xi(3)**2
      if(ixl(17).ne.0) shp(17) = 0.50*f1*f2*(1.+nxi(3,17)*xi(3))
      if(ixl(22).ne.0) shp(22) = 0.50*f1*f2*(1.+nxi(3,22)*xi(3))

c                                                  midside nodes m
      if(ixl(23).ne.0) shp(23) = 0.50*f3*f1*(1.+nxi(2,23)*xi(2))
      if(ixl(25).ne.0) shp(25) = 0.50*f3*f1*(1.+nxi(2,25)*xi(2))

      if(ixl(24).ne.0) shp(24) = 0.50*f2*f3*(1.+nxi(1,24)*xi(1))
      if(ixl(26).ne.0) shp(26) = 0.50*f2*f3*(1.+nxi(1,26)*xi(1))

c.... form relative shape functions - level 3
c.... convert to absolute shape functions
c
c.... (1) modify for center node
      if(ixl(27).ne.0) then
         shp(27) = f1*f2*f3
         f1 = 0.125d0*shp(27)
         f2 = 0.250d0*shp(27)
         do i = 1,4
           shp(i  ) = shp(i  ) - f1                        ! corner b 1,2,3,4
           shp(i+4) = shp(i+4) - f1                        ! corner t 5,6,7,8
           if(ixl(i+ 8).ne.0) shp(i+ 8) = shp(i+ 8) - f2   ! corner m 9,10,11,12
cww        if(ixl(i+18).ne.0) shp(i+18) = shp(i+18) - f2   ! mid    m 19,20,21,22 ???? original
           if(ixl(i+17).ne.0) shp(i+17) = shp(i+17) - f2   ! mid    m 18,19,20,21  modified ww
           if(ixl(i+12).ne.0) shp(i+12) = shp(i+12) - f2   ! mid    m 13,14,15,16  added ww

           if(ixl(  17).ne.0) shp(  17) = shp(  17) - f2   ! warum in Schleife ???
cww        if(ixl(  26).ne.0) shp(  26) = shp(  26) - f2   ! node number und warum in Schleife ???
           if(ixl(  22).ne.0) shp(  22) = shp(  22) - f2   ! node number passend zu 17,   warum in Schleife ???
        end do
      end if
c
c.... (2) modify for mid-face nodes at bottom/top
      if(ixl(17).ne.0) then
        f1 = 0.5d0*shp(17)
        if(ixl(13).ne.0) shp(13) = shp(13) - f1
        if(ixl(14).ne.0) shp(14) = shp(14) - f1
        if(ixl(15).ne.0) shp(15) = shp(15) - f1
        if(ixl(16).ne.0) shp(16) = shp(16) - f1
        f1 = 0.5d0*f1
        shp(1) = shp(1) - f1
        shp(2) = shp(2) - f1
        shp(3) = shp(3) - f1
        shp(4) = shp(4) - f1
      end if

      if(ixl(22).ne.0) then
        f1 = 0.5d0*shp(22)
        if(ixl(18).ne.0) shp(18) = shp(18) - f1
        if(ixl(19).ne.0) shp(19) = shp(19) - f1
        if(ixl(20).ne.0) shp(20) = shp(20) - f1
        if(ixl(21).ne.0) shp(21) = shp(21) - f1
        f1 = 0.5d0*f1
        shp(5) = shp(5) - f1
        shp(6) = shp(6) - f1
        shp(7) = shp(7) - f1
        shp(8) = shp(8) - f1
      end if

      do i = 1,4
        j  = (mod(i,4)) + 1                               ! 2,3,4,1

        if(ixl(22+i).ne.0) then                           ! 23-26
          f1 = 0.5d0*shp(22+i)
          if(ixl(i+ 8).ne.0) shp(i+ 8) = shp(i+ 8) - f1   !  9-12 cm
          if(ixl(j+ 8).ne.0) shp(j+ 8) = shp(j+ 8) - f1   ! 10,11,12,9 cm
          if(ixl(i+12).ne.0) shp(i+12) = shp(i+12) - f1   ! 13-16 mb
          if(ixl(i+17).ne.0) shp(i+17) = shp(i+17) - f1   ! 18-21 mt
          f1 = 0.5d0*f1
          shp(i  ) = shp(i  ) - f1                        ! 1,2,3,4
cww       shp(i  ) = shp(i  ) - f1 ! original but wrong
          shp(j  ) = shp(j  ) - f1 ! ok                   ! 2,3,4,1
          shp(i+4) = shp(i+4) - f1                        ! 5,6,7,8
          shp(j+4) = shp(j+4) - f1                        ! 6,7,8,5
        end if
c
c.... (3) modify for mid-edge nodes
        if(ixl( 8+i).ne.0) then                           ! 9-12 cm
          f1 = 0.25d0*shp( 8+i)
          if(ixl(i+22).ne.0) shp(i+22) = shp(i+22) - f1
          shp(i  ) = shp(i  ) - f1
          shp(i+4) = shp(i+4) - f1
        end if

        if(ixl( 8+j).ne.0) then                            ! 10,11,12,9 cm
          f1 = 0.25d0*shp( 8+j)
          if(ixl(i+22).ne.0) shp(i+22) = shp(i+22) - f1
          shp(j  ) = shp(j  ) - f1
          shp(j+4) = shp(j+4) - f1
        end if

        if(ixl(12+i).ne.0) then                            ! 13-16 mb
          f1 = 0.5d0*shp(12+i)
          if(ixl(i+22).ne.0) shp(i+22) = shp(i+22) - f1
          shp(i  ) = shp(i  ) - f1
          shp(j  ) = shp(j  ) - f1
        end if

        if(ixl(17+i).ne.0) then                            ! 18-21 mt
          f1 = 0.5d0*shp(17+i)
          if(ixl(i+22).ne.0) shp(i+22) = shp(i+22) - f1
          shp(i+4) = shp(i+4) - f1
          shp(j+4) = shp(j+4) - f1
        end if
      end do
      end
c
      subroutine pqrblk(op,oq,or, np1,nq1,nr1, ne,ni, ma, ix,nen1)
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate a block of p x q x r Lagrangian elements
c
c      Inputs:
c         op         - Order of element in p-direction
c         oq         - Order of element in q-direction
c         or         - Order of element in r-direction
c         np1        - Number nodes in r-direction
c         nq1        - Number nodes in s-direction
c         nr1        - Number nodes in t-direction
c         ne         - Initial element number
c         ni         - Initial node number
c         ma         - Material number of block

c      Outputs:
c        ix(nen1,*)  - Element connection list
c        ne          - Final element number
c
c      Comments:     - copy from FEAP8.3 (ww-09/2012)
c
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    op,oq,or, np,nq,nr, ne,ni, ma, nen1
      integer    ix(nen1,*)

      integer    i,j,k, ii,jj,kk, nii,njj,nkk, ij
      integer    ne1,ne2,ne3, np1,nq1,nr1

c     Generate elements pxqxr type

      np  = np1 - 1
      nq  = nq1 - 1
      nr  = nr1 - 1

      ne1 = np/op
      ne2 = nq/oq
      ne3 = nr/or

      nii = ni - 1
      do k = 1,ne3
        do j = 1,ne2
          do i = 1,ne1
            njj = np1*(nq1*or*(k-1) + oq*(j-1)) + op*(i-1) + nii
            ij  = 0
            do kk = 1,or+1
              nkk = njj
              do jj = 1,oq+1
                do ii = 1,op+1
                  ix(ii+ij,ne) = nkk + ii
                end do ! ii
                ij = ij + op + 1
                nkk = nkk + np1
              end do ! jj
              njj = njj + np1*nq1
            end do ! kk
            ix(nen1  ,ne) = ma   ! Material number
            ne = ne + 1
          end do ! i
        end do ! j
      end do ! k
      ne = ne - 1

      end
c
      subroutine writer(ct,b,nneq)
c----------------------------------------------------------------------
c
c      Purpose: Save nodal displacement and stress values for later use
c
c      Inputs:
c         ct        - Name of array, file rewind, or file name
c         b(*)      - displacement array to write
c         nneq     - Length of array b
c
c      Outputs:
c         none
c
c----------------------------------------------------------------------
      USE cdata
      USE fdata
      USE iodata
      USE iofile
      USE pdata3
      USE psize
      USE strnam
      logical lflg
      character*4 ct,fname1,y*1
      double precision b(1)
      save lflg
      data lflg/.false./
      if(ct.eq.'disp') go to 100
      if(ct.eq.'stre') go to 200
      if(ct.eq.'wind') go to 300
      if(ct.eq.'clos') go to 400
c.... set filename
      fname1 = ct
      inquire(file=fname1,exist=lflg)
      if(lflg) then
c....   this is a old file name
          if(ior.lt.0) then
            write(*,2060) fname1
10          read (*,1000,err=11,end=12) y
            goto 13
11          call errclr ('WRITER')
            goto 10
12          call endclr ('WRITER',y)
13          if(y.ne.'y' .or. y.ne.'Y') return
          else
            write(iow,2061) fname1
cww         stop
          return
          end if
          open(iwd,file=fname1,status='old',form='unformatted')
          rewind iwd
      else
c....   this is a new file name
          if(ior.lt.0) write(*  ,2001) fname1
                       write(iow,2001) fname1
          open(iwd,file=fname1,status='new',form='unformatted')
      end if
      lflg = .true.
      return
c.... save current displacement state
100   if(lflg) then
          write(iwd,err=980) ct
          write(iwd,err=980) (b(i),i=1,nneq)
      else
          go to 990
      end if
      return
c.... save current nodal stress state
200   if(lflg) then
          if(fl(11)) then
            write(iwd,err=980) ct
            write(iwd,err=980) istv,(strea,i=1,size(strea))
          else
                         write(iow,2070)
            if(ior.lt.0) write(*  ,2070)
          end if
      else
          go to 990
      end if
      return
c.... rewind file
300   if(lflg) then
          rewind iwd
      end if
      return
c.... close file
400   close(iwd)
      lflg = .false.
      return
980   write(iow,2080) ct
      if(ior.gt.0) stop 'SR WRITER 980'
      write(*,2080)
      return
990   write(iow,2090)
      if(ior.gt.0) stop 'SR WRITER 990'
      write(*,2090)
      return
c.... format statements
1000  format(a1)
2001  format('   Output file for write operations is named ',a4)
2060  format(' ** WARNING ** File ',a4,' exists. Erase? (y or n) >',$)
2061  format(' ** ERROR ** File ',a4,' exists.')
2070  format(' ** ERROR ** Nodal stresses do not exist for tape write')
2080  format(' ** ERROR ** on a tape write command for ',a4)
2090  format(' ** ERROR ** No write file is open.')
      end
c
      subroutine writerb(x,ix,ndm,numnp,nen,nen1,numel)
c----------------------------------------------------------------------
c
c      Purpose: Write coordinates and elements to binary file Bname
c               on ios=12
c
c      Inputs:
c          x(ndm,numnp)   - coordinate array
c         ix(nen1,numel)  - element array
c
c      Outputs:
c         none
c
c----------------------------------------------------------------------
      USE comfil
      USE iodata
      USE iofile
      character*229   fbout
      double precision x(ndm,numnp)
      integer ix(nen1,numel)
c
c.... open file bname=iname+1=B
      call dochar2(finp,ipos)
      fbout = finp
      call dochar1(fbout,'B',ipos)
      open(ios,file=fbout,status='unknown',form='unformatted')
c
c.... write coordinates
      write(ios)  'coor'
      write(ios) ((x(i,j),j=1,numnp),i=1,ndm)
c
c.... write elements
      write(ios) 'elem'
      write(ios) ((ix(   i,j),j=1,numel),i=1,nen) ! nodes
      write(ios)  (ix(nen1,j),j=1,numel)          ! mate

c
c.... close file
      close(ios)
      return
      end
c
c
      subroutine xmodif(xg,ndm,numnp,prt)
c----------------------------------------------------------------------
c
c      Purpose:  modify  cartesian coordinates  due to type
c
c....  usage:
c      cmod, 1,na,ne,r    : Type1: sphere         z = sqrt(r**2-x**2-y**2)
c                                  from           x**2+y**2+z**2=r**2
c      cmod, 2,na,ne,r,f  : Type2: parabola       z = [-f/(r**2)](x**2+y**2)+f
c                                  from parabola on circle
c      cmod, 3,na,ne,a    : Type3: hypar          z = a*x*y
c
c      cmod, 4,na,ne,a,c  : Type4: hyperboloid    r = a/c*sqrt(c**2+z**2)
c                                  from           r**2/a**2-z**2/c**2=1
c                                                 coor must be polar
c      cmod, 5,na,ne,r,a  : Type5: hyperboloid    z = r*(sqrt(1+x**2/a**2)-1)
c
c      cmod, 6,na,ne,r,a  : Type6: sphere         z = z + sqrt(r*r-y*y)-a
c
c      cmod, 7,na,ne,a    : Type7: twisted beam p = pi/2/a, y=y*cos p, z=y*sin p
c
c      cmod, 8,na,ne,r,m  : Type1: sphere         1/8 of sphere using 3 blocs,each m*m elements
c                                  bloc1          recompute x,y,z
c                                  bloc2          from spherical coordinates
c                                  bloc3          recompute z from
c                                  from           z = sqrt(r**2-x**2-y**2)
c      use only after coordinates have been defined
c      use eload, edge etc before!!
c
c      Input:
c         xg(ndm,*)  - Original nodal coordinates of mesh
c         ndm        - Spatial dimension of mesh
c         numnp      - Number of nodes in mesh
c         prt        - Print Flag
c
c      Output:
c         xg(ndm,*)  - modified nodal coordinates of mesh
c
c----------------------------------------------------------------------
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension xg(ndm,*),td(5),tc(2)
      character yyy*80,cname(8)*15
      data cname /'sphere         ','parabola       ','hypar          ',
     +            'hyperboloid    ','hyperboloid    ','sphere         ',
     +            'twisted beam   ','sphere-3 blocs '/
c
      pi = 4.d0*datan(1.d0)
      if(ndm.ne.3) then
        write(yyy,'(a)') 'cmod: ndm .ne. 3'
        call drawmess(yyy,1,0)
      end if
100   if(ior.lt.0) write(*,3000)
      call dinput(td,5)
      tc(1) = td(4)
      tc(2) = td(5)
      if(errck) go to 100
      ityp = td(1)
      na   = td(2)
      ne   = td(3)
      if(na.le.0) na=1
      if(ne.le.0) ne=numnp
      na = min(na,numnp)
      ne = min(ne,numnp)
      if(ityp.eq.1) then
c....   input sphere
        r  = tc(1)
      else if(ityp.eq.2) then
c....   input parabola
        r  = tc(1)
        f  = tc(2)
      else if(ityp.eq.3) then
c....   input hypar
        a  = tc(1)
      else if(ityp.eq.4) then
c....   input hyperboloid
        a  = tc(1)
        c  = tc(2)
      else if(ityp.eq.5) then
c....   input hyperboloid
        r  = tc(1)
        a  = tc(2)
      else if(ityp.eq.6) then
c....   input sphere
        r  = tc(1)
        a  = tc(2)
      else if(ityp.eq.7) then
c....   input twisted beam
        a  = tc(1)
      else if(ityp.eq.8) then
c....   input sphere with 3 blocs
        r  = tc(1)
        m  = tc(2)
        j  = m+1
        j3 = j*j
        j1 = j3-m
        j2 = j1+m/2
        j4 = 2*j3-j
        theta1 = dacos(xg(3,j1)/r)
        theta2 = dacos(xg(3,j2)/r)
        theta3 = dacos(xg(3,j3)/r)
        b = 4.d0*theta2-3.d0*theta1-theta3
        c = theta3-theta1-b
        pi4 = datan(1.d0)
      end if
c.... input error: wrong type
        if(ityp.lt.1.or.ityp.gt.8) then
          write(yyy,'(a)') 'cmod wrong ityp'
          call drawmess(yyy,1,0)
        end if
c...  calculate new values
      do n = na,ne
        x = xg(1,n)
        y = xg(2,n)
        z = xg(3,n)
        if(ityp.eq.1) then                     ! sphere
          d =  r*r - x*x - y*y
          if(d.lt.0.d0) d=0.d0
          xg(3,n) = dsqrt(d)
        else if(ityp.eq.2) then                ! parabola
          xg(3,n) = -f/(r*r) * (x*x + y*y) + f
        else if(ityp.eq.3) then                ! hypar
          xg(3,n) = a*x*y
        else if(ityp.eq.4) then                ! hyperboloid (polar)
          xg(1,n) = a/c*dsqrt(c*c+z*z)
        else if(ityp.eq.5) then                ! hyperboloid
          xg(3,n) = r*(sqrt(1+x*x/a**2)-1)
        else if(ityp.eq.6) then                ! sphere
          xg(3,n) = xg(3,n)+ sqrt(r*r-xg(2,n)**2)-a
        else if(ityp.eq.7) then                ! twisted beam
          p = pi/2.d0/a*x                      ! actual angle
          xg(2,n) = y*dcos(p)
          xg(3,n) = y*dsin(p)
        else if(ityp.eq.8) then                ! sphere 3 blocs
          if(n.le.j3)then                      ! bloc1
            k = (n-1)/j + 1                    ! row
            l =  n - (k-1)*j                   ! column
            phi = (45.d0/m)*(l-1)
            xi = phi/45.d0
            phib = phi *pi4/45.d0
            theta0 = theta1 + b*xi +c*xi*xi
            dtheta = (2.d0*pi4 -theta0)/m
            thetab = (2.d0*pi4)-(k-1)*dtheta

            xg(1,n) = r*dsin(thetab)*dcos(phib)
            xg(2,n) = r*dsin(thetab)*dsin(phib)
            xg(3,n) = r*dcos(thetab)

          else if(n.gt.j3.and.n.le.j4)then     ! bloc2

            nj3 = n - j3
            k = (nj3-1)/m  + 1                 ! row
            l =  nj3 - (k-1)*m +1              ! column
            phi = (45.d0/m)*(l-1)
            xi = 1.d0-phi/45.d0
            phi = phi + 45.d0
            phib = phi*pi4/45.d0
            theta0 = theta1 + b*xi +c*xi*xi
            dtheta = (2.d0*pi4 -theta0)/m
            thetab  = (2.d0*pi4)-(k-1)*dtheta

            xg(1,n) = r*dsin(thetab)*dcos(phib)
            xg(2,n) = r*dsin(thetab)*dsin(phib)
            xg(3,n) = r*dcos(thetab)

          else if(n.gt.j4)then                 ! bloc3
            d =  r*r - x*x - y*y
c           if(d.lt.0.d0) d=0.d0
            xg(3,n) = dsqrt(d)
          end if
        end if
      end do
c.... print new values
      mct = 0
      do n = 1, numnp
        if(mct.gt.0) go to 200
        if(prt)             write(iow,2000) ityp,cname(ityp),(i,i=1,ndm)
        if(prt.and.ior.lt.0)write(*,  2000) ityp,(i,i=1,ndm)
        mct = 100
200     if(prt)              write(iow,2001) n,(xg(i,n),i=1,ndm)
        if(prt.and.ior.lt.0) write(*,  2001) n,(xg(i,n),i=1,ndm)
        mct = mct - 1
      end do
      return
c.... formats
2000  format(/,
     1 '  New cartesian coordinates computed from CMOD',/,
     2 '  Typ ',i3,' for ',a15,/,
     3 4x,'node',3(i6,'-coord'))
2001  format(i8,3e12.4)
3000  format(' Input: Type,1. and last node, parameter  >',$)
      end
