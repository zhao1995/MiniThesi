c----------------------------------------------------------------------
c
      subroutine plotn(v,x,ix,idis,sv,nxy,k1,k2,k3,ndf,nen,nen1,ndm,
     +                 numnp,numnp2,numel,xc,yc,ck)
c----------------------------------------------------------------------
c
c     Purpose: plot displacement / stresses along a line for
c              plate and shells
c              3d problems are possible on the faces if hide is active
c
c     * line is scaled to 0-1 for intersected part of line
c       and  not between the marked points
c
c    *  intersection algorithm bases on projection of
c       coordinates in 1-2 direction,  3-d formulation is
c       possible, but note that now hidden line technic is used
c
c    * the intersection bases only on the 4 corner nodes
c    * 8/9 nodes are possible but the additional nodes are not used
c
c    * possible screen macros: move
c                              rot
c                              iso
c                              pers
c    * storage is given for numnp2=0.5*numnp points on line and
c                                     numel intersected elements
c
c.... Inputs:
c     x          = coordinate field
c     v          = field of displ./stress at given dof k1
c     ix         = field of node numbers at element
c     idis       = field of new element numbers
c     k1         = number of displacment/stress
c                < 0 on deformed mesh
c     k2         = 0 plot new   value in diagram
c                = 1 plot other value in same diagram  with    rescaling
c                = 2 plot other value in same diagram  without rescaling
c                < 0 do not plot, print coordinates
c     k3         = no. of  eigenvector(EPLO)  or layer(SPLO)
c     xc,yc      = s0(1),s0(2)
c     ck         = d(isplacement)/t(ime)/s(tress)/r(eac)/e(igv) in legend
c     sne(2,nen) = dl,dn
c     sv(2,numnp2)= coordinates for plot(x,y)
c     nxy(numel) = array with numbers of intersected elements
c     xo         = coordinates of line points (screen/feap values)
c     npt        = number of intersected elements
c     npv        = number of points on line
c
c     Outputs
c
c     open problems/ideas
c     * clip,pers
c     * disp/vs time see timeplot.for but macro tplo!
c     * macro eplo for eigv similar to displacements possible
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      USE dspos
      USE fdata
      USE iofile
      USE pdata1
      USE pdata2
      implicit double precision(a-h,o-z)
      logical flge
      character*1 ck
      real*8 xo(4),x(ndm,numnp),v(ndf,numnp),sne(60),sv(numnp)
      integer*4  nxy(numel),ix(nen1,numel),idis(*)
      save xo,npt
      if(nen.gt.30) then
          write(*,1000) nen
          return
      end if

c.....for repetition without new calculation
      if(k2.gt.0) go to 110

c.... get coordinates of line points 
      dscorval=dot(dscor,dscor,6)
      if (dscorval.eq.0.d0) then 
c....   get coordinates from screen interactive
        call plopen
        call getlcor1(xo)
        call plclos
        if(xo(1).eq.xo(3) .and. xo(2).eq.xo(4)) then
          call drawmess('No line input',7,1)
        return
        end if
      else 
c....   values are set via SPLO,SET
        call getlcor2(dscor,xo,ndm,scale,s0,sx,dx,iso)
      end if

c      if(pfr .and. ior.lt.0) write(*,2000) xo(1),xo(2),xo(3),xo(4)

      iclear=0
      call plopen
      npt    = 0
c.... loop over all elements
      do 10 n = 1,numel
          new=idis(n)
          ma = ix(nen1,new)
          if(iplma(ma).eq.0) go to 10
          flge = .false.
c....     find intersected elements (flge = .false.)
          call findxy(xo,sne,x,ix(1,n),nen,ndm,scale,s0,dx,sx,iso,flge)
          if(flge) then
            npt = npt + 1
            nxy(npt) = n
          end if
10    continue
      if(npt.gt.numel) then
          write(*,1001)
          return
      end if
c.... if no intersection points found
      if(npt.eq.0) then
          call drawmess('No points found',7,1)
          return
      end if
110   continue
      npv = 0
c.... loop over all intersected elements
      do n = 1,npt
          ne = nxy(n)
          flge = .true.
c....     work on intersected elements (flge = .true.)
          call findxy(xo,sne,x,ix(1,ne),nen,ndm,scale,s0,dx,sx,iso,flge)
c....     calculate coordinates of points
          call findpt(v,k1,ck,sne,ix(1,ne),ndf,nen,sv,npv)
      end do
      if(npv.gt.numnp2) then
          write(*,1002)
          return
      end if
c.... sort the values into increasing order
      call sortpt(sv,npv)
      if(k2.lt.0) then
c....     print the values on the line
          call princl(sv,npv,k1,k3,ck)
      else
c....     plot the values on the line
          call plotcl(sv,npv,ndf,k1,k2,k3,xc,yc,ck,ipgl)
      end if
      return
c.... formats
1000  format(' field sne to small due to nen in plotn, nen= ',i5)
1001  format(' Too many intersected elements in plotn')
1002  format(' Too many points on line in plotn')
2000  format(' Screen coordinates: x1/y1',2 e12.5,/,
     +       ' Screen coordinates: x2/y2',2 e12.5)
      end
c
c-----------------------------------------------------------------------
c

      subroutine findxy(xo,sne,x,ixl,nen,ndm,scale,s0,dx,sx,iso,flge)
c-----------------------------------------------------------------------
c
c     Purpose:
c      1) find elements which are intersected (flge=.false.)
c      2) store coordinates of intersection points (flge=.true.)
c
c     Inputs:
c     xo         = boundary coordinates of line (screen cor.)
c     x          = coordinate field
c     ixl        = field of node numbers for the actual element
c     sne(2,nen) = dl,dn
c
c     Output
c
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      USE ppers
      implicit double precision(a-h,o-z)
      logical iso,flge,plus,negn
      real*8 xo(4), x(ndm,*),sne(2,*),s0(2),dx(2),sx(2),xg(3)
      integer*4     ixl(nen)
c.... values of line
      x1  = xo(1)
      y1  = xo(2)
      x2  = xo(3)
      y2  = xo(4)
      ddx  = x2 - x1
      ddy  = y2 - y1
      dl2  = ddx*ddx + ddy*ddy
      plus = .false.
      negn = .false.
c.... loop over all element nodes
      do k = 1,nen
          n = ixl(k)
          n = abs(n)
          if(n.gt.0) then
c....       recover coordinates
            xx1 = x(1,n)
            xx2 = x(2,n)
            xx3 = 0.d0
            if(ndm.eq.3) xx3 = x(3,n)
c....       perform perspective tranformation
            if(kpers.eq.1) then
              xg(1) = xx1
              xg(2) = xx2
              xg(3) = xx3
              call perspj(xg,xg,3,3,1,errv)
              xx1 = xg(1)
              xx2 = xg(2)
              xx3 = xg(3)
            end if
c....       compute the screen coordinates of point
            s1 = scale*(xx1 + xx1 -sx(1)) + s0(1)
            s2 = scale*(xx2 + xx2 -sx(2)) + s0(2)
c....       if isometric recompute values
            if(iso) then
              xmul = 2.*max(dx(1),dx(2))/(dx(1)+dx(2))
              te   = (s0(1) + 0.5*(s1 - s2))*xmul
              s2   = (0.2885*(s1 + s2))*xmul + scale*(xx3+xx3)
              s1   = te - 0.1
            end if
c....       check if line passes through the element
            dl = (s1-x1)*ddx + (s2-y1)*ddy
            if(dl.gt.0.0 .and. dl.lt.dl2) then
              dn = (x1-s1)*ddy + (s2-y1)*ddx
              if(dn.ge.0.0) plus = .true.
              if(dn.le.0.0) negn = .true.
            end if
c.....  store points of element on line
            if(flge) then
              sne(1,k) = dl/dl2
              sne(2,k) = dn/dl2
            end if
          end if
      end do
      flge = plus .and. negn
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine findpt(v,k1,ck,sne,ixl,ndf,nen,sv,npv)
c-----------------------------------------------------------------------
c
c.... Purpose: find intersections
c
c     Inputs:
c     v          = field of displ./stress/reac at given dof k1
c     ixl        = field of node numbers for the actual element
c     sne(2,nen) = dl,dn
c     sv (1,   ) = x - values for plot (screen/feap values)
c     sv (2,   ) = y - values for plot        (real values)
c
c     Output
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      integer*4 ixl(nen)
      real*8 sne(2,*),v(ndf,*),sv(2,*)
      character*1 ck
      kp = 1
      if(ck.eq.'d'.or.ck.eq.'r'.or.ck.eq.'e') kp = k1
      do n = 1,nen
          if(ixl(n).ne.0) nel = n
      end do
c.... interpolation only for 4 nodes elements
      nel = min(nel,4)
c.... find the intersections
      n1 = nel
      ne1= abs(ixl(n1))
      do n = 1,nel
          ne = abs(ixl(n))
          if(ne.gt.0) then
            if(sne(2,n1)*sne(2,n).lt.0.0) then
              ratio   = sne(2,n)/(sne(2,n) - sne(2,n1))
              svtry   = sne(1,n) + ratio*(sne(1,n1) - sne(1,n))
              if(npv.gt.0) then
                do i = 1,npv
                      if(abs(svtry - sv(1,i)).lt.1.e-4) go to 110
                end do
              end if
              npv = npv + 1
              sv(1,npv) = svtry
              sv(2,npv) = v(kp,ne) + ratio*(v(kp,ne1) - v(kp,ne))
            end if
110         n1 = n
            ne1= ne
          end if
      end do
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sortpt(sv,npv)
c-----------------------------------------------------------------------
c
c.... Purpose: sort the values in array sv into increasing order
c
c     Inputs:
c
c     Outputs
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      real*8 sv(2,*)
c.... sort on first values
      do n = 1,npv-1
          do i = n+1,npv
            if(sv(1,i).lt.sv(1,n)) then
              sn = sv(1,n)
              sv(1,n) = sv(1,i)
              sv(1,i) = sn
              sn = sv(2,n)
              sv(2,n) = sv(2,i)
              sv(2,i) = sn
            end if
          end do
      end do
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plotcl(sv,npv,ndf,k1,k2,k3,xc,yc,ck,ipgl)
c-----------------------------------------------------------------------
c
c.... Purpose: plot contour line
c
c     Inputs:
c     sv(1,*)   = coordinate values on line (screen values)
c     sv(2,*)   = associated displ./stress value   both for point(x,y)
c     npv       = number of points on line
c     ndf       = dofs of freedom in problem
c     k1        = dof to plot
c     k2        = 0 plot new   value in diagram
c               = 1 plot other value in same diagram  with    rescaling
c               = 2 plot other value in same diagram  without rescaling
c               < 0 do not plot, print coordinates
c     k3        = no. of  eigenvector(EPLO)  or layer(SPLO)
c     xc,yc     = s0(1),s0(2)
c     ck        = d(isplacement)/t(ime)/s(tress)/r(eac)/e(igv) in legend
c     ipgl      = 1 - fill box, = 3 - line box
c
c     Outputs
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      character ck*1,yy*80
      real*8    sv(2,*)
      save vmin,vmax
      dx1 = 0.01d0
c...  calculate extremal values in x dir
      smin = sv(1,1)
      smax = smin
      do n = 1,npv
          smin = min(smin,sv(1,n))
          smax = max(smax,sv(1,n))
      end do
c...  calculate extremal values in y dir
      vminl = sv(2,1)
      vmaxl = vminl
      do n = 1,npv
          vminl = min(vminl,sv(2,n))
          vmaxl = max(vmaxl,sv(2,n))
      end do
      if(k2.lt.2) then
          vmin = vminl
          vmax = vmaxl
      end if
c.... plot box and grid
      if(k2.eq.0) call plota(vmin,vmax,xc,yc)
      if(npv.ge.2) then
c....     put extremal values on screen
          y0 = 0.76d0
          dy = 1./18.
          if(ck.eq.'d') write(yy,1004) ' DISPLACEMENT  '
          if(ck.eq.'r') write(yy,1004) ' REACTION      '
          if(ck.eq.'s') write(yy,1002) abs(k3)
          if(ck.eq.'e') write(yy,1003) abs(k3)
          call drawtxt(1,1.d0,y0,1,1,15,yy)
          write(yy,1000) ck,k1,vminl
          kp = k1
          if(kp.gt.10) kp = 10
          ycor = 0.76d0 - kp*dy
          call drawtxt(1,1.d0,ycor,k1,1,20,yy)
          ycor = 0.76d0 - (kp+0.5d0)*dy
          write(yy,1001) ck,k1,vmaxl
          call drawtxt(1,1.d0,ycor,k1,1,20,yy)
          pp   = 0.6d0
          ph   = 0.3d0
          fact = 1.0d0
          vminl= vmin
          if(xc.ne.0.5d0) fact = 0.5d0
          dx  = (smax - smin)/pp
          dy  = (vmax - vmin)/pp
          if(dy.le.0.0d0) then
            dy = 1.0d0
            vmin = vmin -  ph
          end if
          ii = 2
          if(k1.gt.0) ii = 1
          do 210 i = ii,2
            j   = 3
            do 200 n = 1,npv
              vx = xc + ((sv(1,n) - smin)/dx - ph)*fact
c.....    test if y-value is on screen
              yval=sv(2,n)
              if(yval.lt.vminl) then
                yval = vminl
                call pppcol(32)
              else if(yval.gt.vmax) then
                yval = vmax
                call pppcol(32)
              end if
              vy = yc + ((yval - vmin)/dy - ph)*fact
c....         plot line and actual point
              call dplot(vx,vy,j,0)
              call ppbox1(vx,vy,dx1,dx1,ipgl)
              if(i.eq.2) then
                j  = 2
              else
c.....          plot label
c               call plabl(k1)
              end if
c.....        reset color
              call pppcol(k1)
200         continue
210       continue
      end if
c.... plot zero line if in window in border color
      if(vmin*vmax .lt. 0.0) then
          call pppcol(8)
          vx = xc - ph*fact
          vy = yc - (vmin/dy + ph)*fact
          call dplot(vx,vy,3,0)
          vx = xc + ((smax - smin)/dx - ph)*fact
          call dplot(vx,vy,2,0)
      end if
      return
1000  format(a1,i2,1x,'min',e10.3)
1001  format(a1,i2,1x,'max',e10.3)
1002  format('STRESS-Layer',i3)
1003  format('EIGENVECTOR ',i3)
1004  format(a15)
      end
c
c-----------------------------------------------------------------------
c
      subroutine plota(vmin,vmax,xc,yc)
c-----------------------------------------------------------------------
c
c.... Purpose: plot box for splo/dplo
c              add points and plot the vertical and horizontal legends
c
c     Inputs:
c     vmin
c     vmax
c     xc,yc  = center of box
c
c     Outputs
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      character yy*10
c.... put a box around a plot
      if(xc.eq.0.5 .and. yc.eq.0.5) then
          dx = .02
          dy = .06
          x0 = 0.2
          y0 = 0.2
          xh = .8
          yh = .8
      else
          dx = .01
          dy = .04
          x0 = 0.05
          y0 = 0.05
          if(xc.gt.0.5) x0 = 0.55
          if(yc.gt.0.5) y0 = 0.55
          xh = 0.4 + x0
          yh = 0.4 + y0
      end if
c.... add labels to the lower axis to show line direction
cww   call drawtxt(1,x0-dy,y0-dy,7,1,7,'Point 1')  ! yellow
      call drawtxt(1,x0-dy,y0-dy,2,1,7,'Point 1')  ! red
cww   call drawtxt(1,xh-dy,y0-dy,7,1,7,'Point 2')  ! yellow
      call drawtxt(1,xh-dy,y0-dy,2,1,7,'Point 2')  ! red

c.... plot  grid
      call pppcol(4) ! cyan
cww   call pppcol(1) ! white
c.... plot the vertical lines
      x  = x0
      do n = 1,11
          call dplot(x,y0-dx,3,0)
          call dplot(x,yh+dx,2,0)
          x = x + dy
      end do
c.... plot the horizontal lines
      y  = y0
      do n = 1,11
          call dplot(x0-dx,y,3,0)
          call dplot(xh+dx,y,2,0)
          y = y + dy
      end do
c.... plot the vertical values
      x = 0.02
      y = y0
      valy = vmin
      dvaly = (vmax - vmin)/10.
      do n = 1,11
          write(yy,1000) valy
          call drawtxt(1,x,y,1,1,10,yy)
          y = y + dy
          valy = valy + dvaly
      end do
c.... plot the horizontal values
      y = 0.17
      x = x0-dx
      valx = 0.
      dvalx = 0.2
      do n = 1,6
          write(yy,1001) valx
          call drawtxt(1,x,y,1,1,3,yy)
          x = x + 2.*dy
          valx = valx + dvalx
      end do

c.... plot a surrounding box
      call pppcol(7) ! yellow
      call ppbox(x0,y0,0.6d0,0.6d0,3)  ! box border in yellow
cww   call ppbox(x0,y0,0.6d0,0.6d0,1)  ! box filled in yellow

1000  format(1pe10.3)
1001  format(f3.1)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine princl(sv,npv,k1,k3,ck)
c-----------------------------------------------------------------------
c
c.... Purpose: print contour line values
c
c     Inputs:
c     sv(1,npv)  = coordinate values on line (screen values)
c     sv(2,npv)  = associated displ./stress value   both for point(x,y)
c     k1         = dof to plot
c     k3         = no. of  eigenvector(EPLO)  or layer(SPLO)
c     ck         = d(isplacement)/s(tress)/r(eac)/e(igv)
c
c     Outputs
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision(a-h,o-z)
      character ck*1
      real*8    sv(2,*),dsv(2)
c.... modify values in sv(1,..) between 0 and 1!
      smin = sv(1,1)
      smax = smin
      do n = 1,npv
          smin = min(smin,sv(1,n))
          smax = max(smax,sv(1,n))
      end do
      ds = smax-smin
c
      if(ck.eq.'d') then
          if(ior.lt.0) write(*  ,1000) k1
                           write(iow,1000) k1
      else if(ck.eq.'r') then
          if(ior.lt.0) write(*  ,1004) k1
                           write(iow,1004) k1
      else if(ck.eq.'e') then
          if(ior.lt.0) write(*  ,1005) k1,k3
                           write(iow,1005) k1,k3
      else if(ck.eq.'s') then
        if(ior.lt.0) write(  *,1003) k1,k3
                     write(iow,1003) k1,k3
      end if
c
      do n = 1,npv
        dsv(1) = 1.d0-(smax-sv(1,n))/ds
        dsv(2) = sv(2,n)
          if(ior.lt.0) write(  *,1002) (dsv(i),i=1,2)
                       write(iow,1002) (dsv(i),i=1,2)
      end do
      return
1000  format(/,'  DISPLACEMENT ',i3,' at Contour Line',/,
     1         '  s(0-1)     value')
1001  format(/,'  STRESS ',i3,' for Layer ',i2,' at Contour Line',/,
     1         '  s(0-1)     value')
1002  format(2x,f8.4,3x,e12.5)
1003  format(/,'  STRESS ',i3,' for Layer ',i2,' at Contour Line',/,
     1         '  s(0-1)     value')
1004  format(/,'  REACTION     ',i3,' at Contour Line',/,
     1         '  s(0-1)     value')
1005  format(/,'  DISPLACEMENT ',i3,'of EIGENVECTOR ',i3,
     1         ' at Contour Line',/,'  s(0-1)     value')
      end
c
c-----------------------------------------------------------------------
c
      subroutine ppbox1(x,y,dx,dy,is)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot a contour or filled box with  center at x,y
c
c     Inputs
c       x
c       y
c       dx
c       dy
c       is = 1 - filled box
c       is = 3 - contour box
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c
      integer is
      real*8 x,y,dx,dy,dx2,dy2
      dx2=0.5d0*dx
      dy2=0.5d0*dy
      call dplot(x-dx2,y-dy2,is,0)
      call dplot(x+dx2,y-dy2, 2,0)
      call dplot(x+dx2,y+dy2, 2,0)
      call dplot(x-dx2,y+dy2, 2,0)
      call dplot(x-dx2,y-dy2, 2,0)
      if(is.eq.1) call clpan
c.... return to starting point
      call dplot(x,y,3,0)
      end
c
c-----------------------------------------------------------------------
c
      subroutine plotxyl1(xy,npv,ispv,xc,yc)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot xy line in a diagram
c
c     Inputs
c     xy(2,numnp2)= coordinates for plot(x,y)
c     xy(1,npv)   = x-value  (real values)
c     xy(2,npv)   = y-value  (real values)
c     npv         = number of points <=0.5*numnp*npstr=numnp*13
c     ispv        = no of stress
c     xc,yc       = s0(1),s0(2)
c
c     Outputs
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
c
      USE pdata2
      implicit double precision(a-h,o-z)
      character yy*80
      real*8 xy(2,npv)

      dx1 = 0.01d0

c...  calculate extremal values in x-dir
      xmin = xy(1,1)
      xmax = xmin
      xsum = 0.d0
      do n = 1,npv
        xmin = min(xmin,xy(1,n))
        xmax = max(xmax,xy(1,n))
        xsum = xsum + dabs(xy(1,n))
      end do

c...  check
      if(xsum.lt.1.e-8) then
        write(*,*) 'All values zero'
        return
      end if

      if(ispv.lt.1.or.ispv.gt.8) then
        write(*,*) 'Only stresses 1-8 allowed'
        return
      end if

c...  calculate extremal values in y-dir
      ymin = xy(2,1)
      ymax = ymin
      do n = 1,npv
        ymin = min(ymin,xy(2,n))
        ymax = max(ymax,xy(2,n))
      end do

c.... plot box and grid
      iclear=0
      call plopen
      call plotxyl2(xmin,xmax,ymin,ymax,xc,yc)

      if(npv.ge.2) then
c....   draw axial values

        if(ispv.eq.1) write(yy,1004) ' Stress S_xx   '
        if(ispv.eq.2) write(yy,1004) ' Stress S_yy   '
        if(ispv.eq.3) write(yy,1004) ' Stress S_zz   '
        if(ispv.eq.4) write(yy,1004) ' Stress T_xy   '
        if(ispv.eq.5) write(yy,1004) ' Stress T_xz   '
        if(ispv.eq.6) write(yy,1004) ' Stress T_yz   '
        if(ispv.eq.7) write(yy,1004) ' Warping Phi_x '
        if(ispv.eq.8) write(yy,1004) ' Warping Phi_y '
        x = 0.45d0
        y = 0.15d0
        call drawtxt(1,x,y,1,1,15,yy)
        write(yy,1004) ' Coordinate z  '
        x = 0.02d0
        y = 0.85d0
        call drawtxt(1,x,y,1,1,15,yy)

c....   draw extremal values
        y0 = 0.76d0
        dy = 1./18.

        write(yy,1000) xmin
        y = 0.76d0 - 0.5d0*dy
        call drawtxt(1,1.d0,y,1,1,15,yy)

        write(yy,1001) xmax
        y = y - 0.5d0*dy
        call drawtxt(1,1.d0,y,1,1,15,yy)

        write(yy,1002) ymin
        y = y - 0.5d0*dy
        call drawtxt(1,1.d0,y,1,1,15,yy)

        write(yy,1003) ymax
        y = y - 0.5d0*dy
        call drawtxt(1,1.d0,y,1,1,15,yy)
c
c....   plot line and points
        pp   = 0.6d0
        ph   = 0.3d0
        fact = 1.0d0
        if(xc.ne.0.5d0) fact = 0.5d0
        dx  = (xmax - xmin)/pp
        dy  = (ymax - ymin)/pp
        if(dx.le.0.0d0) then
          dx = 1.0d0
          xmin = xmin -  ph
        end if
        if(dy.le.0.0d0) then
          dy = 1.0d0
          ymin = ymin -  ph
        end if

        do 210 i = 1,2
          j   = 3
          do 200 n = 1,npv
            vx = xc + ((xy(1,n) - xmin)/dx - ph)*fact
            vy = yc + ((xy(2,n) - ymin)/dy - ph)*fact
            call pppcol(2)
            call dplot(vx,vy,j,0)
            call ppbox1(vx,vy,dx1,dx1,ipgl)
            if(i.eq.2) then
              j  = 2
            else
c.....        plot label
cww           call plabl(ispv)
            end if
200       continue
210     continue
      end if
c.... plot zero lines if in window
      call pppcol(8)

      if(xmin*xmax .lt. 0.0) then
        vy = yc - ph*fact
        vx = xc - (xmin/dx + ph)*fact
        call dplot(vx,vy,3,0)
        vy = yc + ((ymax - ymin)/dy - ph)*fact
        call dplot(vx,vy,2,0)
      end if

      if(ymin*ymax .lt. 0.0) then
        vx = xc - ph*fact
        vy = yc - (ymin/dy + ph)*fact
        call dplot(vx,vy,3,0)
        vx = xc + ((xmax - xmin)/dx - ph)*fact
        call dplot(vx,vy,2,0)
      end if
      return
1000  format(1x,'Xmin',e10.3)
1001  format(1x,'Xmax',e10.3)
1002  format(1x,'Ymax',e10.3)
1003  format(1x,'Ymax',e10.3)
1004  format(a15)
      end
c
c-----------------------------------------------------------------------
c
      subroutine plotxyl2(xmin,xmax,ymin,ymax,xc,yc)
c-----------------------------------------------------------------------
c
c.... Purpose: plot box for xy-line
c              add points and plot the vertical and horizontal legends
c
c     Inputs:
c     xmin   = min x-value
c     xmax   = max x-value
c     ymin   = min y-value
c     ymax   = max y-value
c     xc,yc  = center of box
c
c     Outputs
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      character yy*10
c.... put a box around a plot
      if(xc.eq.0.5 .and. yc.eq.0.5) then
          dx = .02
          dy = .06
          x0 = 0.2
          y0 = 0.2
          xh = .8
          yh = .8
      else
          dx = .01
          dy = .04
          x0 = 0.05
          y0 = 0.05
          if(xc.gt.0.5) x0 = 0.55
          if(yc.gt.0.5) y0 = 0.55
          xh = 0.4 + x0
          yh = 0.4 + y0
      end if
c.... plot  grid
      call pppcol(4) ! cyan

c.... plot the vertical lines
      x  = x0
      do n = 1,11
          call dplot(x,y0-dx,3,0)
          call dplot(x,yh+dx,2,0)
          x = x + dy
      end do

c.... plot the horizontal lines
      y  = y0
      do n = 1,11
          call dplot(x0-dx,y,3,0)
          call dplot(xh+dx,y,2,0)
          y = y + dy
      end do

c.... plot a surrounding box
      call pppcol(7) ! yellow
      call ppbox(x0,y0,0.6d0,0.6d0,3)  ! box border in yellow

c.... plot the vertical values
      x = 0.02
      y = y0+0.01
      valy  =  ymin
      dvaly = (ymax - ymin)/10.
      do n = 1,11
          write(yy,1000) valy
          call drawtxt(1,x,y,1,1,10,yy)
          y = y + dy
          valy = valy + dvaly
      end do
c.... plot the horizontal values (in 2 lines)
      y  = 0.17
      x  = x0-0.05
      dx = 0.06
      dy = 0.015
      valx  =  xmin
      dvalx = (xmax - xmin)/10.
      nn = 0
      do n = 1,11
          write(yy,1001) valx
          call drawtxt(1,x,y,1,1,10,yy)
          x    = x + dx
          if(nn.eq.0) then
             y    = y + dy
             nn = 1
          else
             y    = y - dy
             nn = 0
          end if
          valx = valx + dvalx
      end do
      return
1000  format(1pe10.3)
1001  format(1pe10.3)
      end
c
      subroutine getlcor2(xc,xo,ndm,scale,s0,sx,dx,iso)
c-----------------------------------------------------------------------
c
c     Purpose:
c      compute screen coordinates of line from points set by SPLO,SET
c
c     Inputs:
c     xc(3,2)    = boundary coordinates of line (cart cor.)
c
c     Output
c     xo(4)      = boundary coordinates of line (screen cor.)
c
c     W. Wagner IBS KIT 03/15
c----------------------------------------------------------------------
      USE pdata11
      USE pltran
      USE ppers
      implicit double precision(a-h,o-z)
      logical iso
      real*8 xo(4),xc(3,2),xcl(3,2),xg(3),s0(2),sx(2),dx(2)  

c.... XFAC
      do n = 1,2 
        xcl(1,n) = xc(1,n)*xfac(1)
        xcl(2,n) = xc(2,n)*xfac(2)
        xcl(3,n) = 0.d0
        if(ndm.eq.3) xcl(3,n) = xc(3,n)*xfac(3)
      end do 

c.... ROT
      call plxtrn(xcl,tra,vr,3,2)

c.... Screen Coordinates
      do n = 1,2 
        xx1 = xcl(1,n)
        xx2 = xcl(2,n)
        xx3 = xcl(3,n)

c....   PERS
        if(kpers.eq.1) then
          xg(1) = xcl(1,n)
          xg(2) = xcl(2,n)
          xg(3) = xcl(3,n)
          call perspj(xg,xg,3,3,1,errv)
          xx1 = xg(1)
          xx2 = xg(2)
          xx3 = xg(3)
        end if
c....   Screen Coordinates of point
        s1 = scale*(xx1 + xx1 -sx(1)) + s0(1)
        s2 = scale*(xx2 + xx2 -sx(2)) + s0(2)
c....   ISO
        if(iso) then
          xmul = 2.*max(dx(1),dx(2))/(dx(1)+dx(2))
          te   = (s0(1) + 0.5*(s1 - s2))*xmul
          s2   = (0.2885*(s1 + s2))*xmul + scale*(xx3+xx3)
          s1   = te - 0.1
        end if
c....   SET
        if(n.eq.1) then
          xo(1)=s1
          xo(2)=s2
        else if(n.eq.2) then
          xo(3)=s1
          xo(4)=s2
        end if 
      end do
      return
      end
c
