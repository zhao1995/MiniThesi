c$Id:$
      subroutine vblkn(nr,ns,nt,xl,x,ixl,dr,ds,dt,
     &                 ni,ndm,ctype,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'cap' option for spherical cap               11/11/2008
c          and 'xcyl' and 'ycyl' for cylindrical sectors
c       2. Add option for generating sperical cap with      13/11/2008
c          specified thickness and constant radius.
c       3. Change print to 1p,5e13.4                        16/04/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate a block of 3-d 8-node brick elements

c      Inputs:
c         nr        - Number elements in 1-local coordinate dir.
c         ns        - Number elements in 2-local coordinate dir.
c         nt        - Number elements in 3-local coordinate dir.
c         xl(ndm,*) - Block nodal coordinate array
c         ixl(*)    - Block nodal connection list
c         dr        - 1-local coordinate increment
c         ds        - 2-local coordinate increment
c         dt        - 3-local coordinate increment
c         ni        - Initial node number for block
c         ndm       - Spatial dimension of mesh
c         ctype     - Type of block coordinates
c         prt       - Output generated data if true
c         prth      - Output title/header data if true

c      Outputs:
c         x(ndm,*)  - Nodal coordinates for block
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat2.h'
      include  'iofile.h'
      include  'trdata.h'

      logical   prt,prth,phd, pcomp
      character xh*6, ctype*15
      integer   ni,ndm,nr,ns,nt,i,j,k,l,m,n,mct, ixl(*)
      real*8    dr,ds,dt, rr,sn2,cn2,sn3,cn3,afac
      real*8    ss(3),xl(3,*),x(ndm,*),xx(3)

      save

      data      xh/' coord'/

c     Check that all corners of brick are defined

      do k = 1,3
        xx(k) = 0.0d0
      end do ! k

      do k = 1,8
        if(ixl(k).ne.k) go to 900
      end do ! k
      call bcor3d(ixl,xl)
      n = ni
      mct = 0
      ss(3) = -1.0d0
      do k = 1,nt
        ss(2) = -1.0d0
        do j = 1,ns
          ss(1) = -1.0d0
          do i = 1,nr

c           Compute coordinates of node

            call xbcor3d(ss,xl, xx)

c           Convert coordinates if necessary

            if(pcomp(ctype,'pola',4) .or. pcomp(ctype,'cyli',4)) then
              call pdegree(xx(2), sn2,cn2)
              rr    = xx(1)
              xx(1) = x0(1) + rr*cn2
              xx(2) = x0(2) + rr*sn2
              xx(3) = x0(3) + xx(3)
            elseif(pcomp(ctype,'sphe',4)) then
              call pdegree(xx(2), sn2,cn2)
              call pdegree(xx(3), sn3,cn3)
              rr    = xx(1)
              xx(1) = x0(1) + rr*cn2*sn3
              xx(2) = x0(2) + rr*sn2*sn3
              xx(3) = x0(3) + rr*cn3
            elseif(pcomp(ctype,'xcyl',3)) then
              xx(3) = sqrt(abs(xx(3)**2 - xx(1)**2))
            elseif(pcomp(ctype,'ycyl',3)) then
              xx(3) = sqrt(abs(xx(3)**2 - xx(2)**2))
            elseif(pcomp(ctype,'cap',3)) then
              if(meshrad.gt.0.0d0) then
                afac  = (meshrad + xx(3))/meshrad
                xx(3) = sqrt(abs(meshrad**2 - xx(1)**2 - xx(2)**2))
                xx(1) = xx(1)*afac
                xx(2) = xx(2)*afac
                xx(3) = xx(3)*afac
              else
                xx(3) = sqrt(abs(xx(3)**2 - xx(1)**2 - xx(2)**2))
              endif
            endif

c           Transform to global coordinates

            do m = 1,ndm
              x(m,n) = xr(m)+tr(m,1)*xx(1)+tr(m,2)*xx(2)+tr(m,3)*xx(3)
            end do ! m

c           Output point

            if(prt) then
               mct = mct + 1
               phd = mod(mct,50).eq.1
               call prtitl(prth.and.phd)
               if(phd) write(iow,2000) (l,xh,l=1,ndm)
               write(iow,2001) n,(x(l,n),l=1,ndm)
               if(ior.lt.0) then
                 if(phd) write(*,2000) (l,xh,l=1,ndm)
                 write(*,2001) n,(x(l,n),l=1,ndm)
               endif
            endif
            n = n + 1
            ss(1) = ss(1) + dr
          end do ! i
          ss(2) = ss(2) + ds
        end do ! j
        ss(3) = ss(3) + dt
      end do ! k

      return

c     Error

900   write(iow,3000) k
      write(ilg,3000) k
      if(ior.lt.0) then
        write(*,3000) k
        return
      endif
      call plstop()

c     Formats

2000  format(/'  N o d a l   C o o r d i n a t e s'//6x,'Node',5(i7,a6))

2001  format(i10,1p,5e13.4)

3000  format(' *ERROR* VBLKN: Block node',i3,' is undefined')

      end
