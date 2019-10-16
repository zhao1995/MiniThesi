c$Id:$
      subroutine plotcl(v,svx,svy,nds,npv,nd,numnp,jc,k,kp,xc,yc,ck)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output the curve for a value for dplot/splot command

c      Inputs:
c         v(nd,*)   - Values at all nodes for plot quantity
c         svx(2,*)  - Sorted positions for plot
c         svy(2,*)  - Sorted values for plot
c         nds       - Dimension of svx,svy values
c         npv       - Number of svx,svy values
c         nd        - Dimension of v array
c         numnp     - Number of nodes
c         jc        - Plot number
c         k         - Component number for plot
c         kp        - Rescale indicator
c         xc,yc     - Coordinates on where to center plot
c         ck        - Plot type character

c      Outputs:
c         none      - Plot output to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character ck*1,yyy*23
      integer   nds,npv,nd,numnp,jc,k,kp
      integer   i,ii,inc,j,n,num
      real*8    xc,yc,vx,vy,smin,smax,vmin,vmax,pp,ph,dx,dy,fact,xp,yp
      real*8    svx(nds,1),svy(nds,1),v(nd,1)

      save

      smin = svx(1,1)
      smax = smin
      vmin = svy(1,1)
      vmax = vmin
      num  = npv

      do n = 1,num
        smin = min(smin,svx(1,n))
        smax = max(smax,svx(1,n))
        vmin = min(vmin,svy(1,n))
        vmax = max(vmax,svy(1,n))
      end do ! n

      if(num.ge.2) then

c       Put label on screen

        call pppcol(max(1,k),-1)
        write(yyy,1000) jc,ck,k,vmin,vmax
        xp = 0.982d0
        yp = 0.782d0 - k*0.05455d0
        call pltext(xp,yp,-23,yyy)

c       Rescale if needed

        if(kp.le.0) then
          do n = 1,numnp
            vmin = min(vmin,v(1,n))
            vmax = max(vmax,v(1,n))
          end do ! n
        endif
        pp  = 0.8
        ph  = 0.4
        fact= 1.0
        if(xc.ne.0.5) fact = 0.5
        dx  = (smax - smin)/pp
        dy  = (vmax - vmin)/pp
        if(dy.le.0.0) then
          dy = 1.0
          vmin = vmin -  ph
        endif
        inc = max(1,num/12)
        ii = 2
        if(k.gt.0) ii = 1

        do i = ii,2
          j   = -3
          do n = 1,num
            vx = xc + ((svx(1,n) - smin)/dx - ph)*fact
            vy = yc + ((svy(1,n) - vmin)/dy - ph)*fact
            call dplot(vx,vy,j)
            if(i.eq.2) then
              j  = -2
            else
              if(mod(n,inc).eq.0) call plabl(-k)
            endif
          end do ! n
        end do ! i

      endif

c     Plot zero line if in window

      if(vmin*vmax .lt. 0.0) then
        vx = xc - ph*fact
        vy = yc - (vmin/dy + ph)*fact
        call dplot(vx,vy,-3)
        vx = xc + ((smax - smin)/dx - ph)*fact
        call dplot(vx,vy,-2)
      endif

c     Clear plot buffer to preserve line color

      call dplot(0.d0,0.d0,-3)
      call dplot(0.d0,0.d0, 3)

1000  format(i2,a1,i2,1p,2e9.2)

      end
