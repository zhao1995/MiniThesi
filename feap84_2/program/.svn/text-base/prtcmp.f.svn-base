c$Id:$
      subroutine prtcmp(x,br,bi,ttim,prop,ndm,ndf,n1,n2,n3,
     &                  ndf1,ndf2,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output nodal values for complex solutions

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         br(*)     - Real part of solution
c         bi(*)     - Imaginary part of solution
c         ttim      - Value of solution time
c         prop      - Value of total proportional load
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         n1        - First node to output
c         n2        - Last noed to output
c         n3        - Increment to n1
c         ndf1      - First component to output
c         ndf2      - Last  component to output
c         prth      - Output title/header data if true

c      Outputs:
c         None      - Outputs to file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'fdata.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'xtout.h'

      logical   prth
      character cd*6,dr*6,di*6,rad*6,ang*6
      integer   ndm,ndf,n1,n2,n3, ndf1,ndf2, i,n, count, nxt1
      real*8    ttim,prop, radii,theta, x(ndm,*),br(ndf,*),bi(ndf,*)

      save

      data      cd/' coord'/, dr, di/'-Real','-Imag'/
      data      rad, ang /'-Radii','-Angle'/

      count = 0
      nxt1  = max(1,nxt)
      do n = n1,n2,n3
        if(nxt.eq.0 .or. (abs(x(nxt1,n)-xt).le.xtol) ) then
          count = count - 1
          if(count.le.0) then
            call prtitl(prth)
            write(iow,2000) ttim,prop,(i,cd,i=1,ndm),1,rad,2,ang,
     &                         (i,dr,i,di,i=ndf1,ndf2)
            if(ior.lt.0.and.pfr) then
              write(*,2000) ttim,prop,(i,cd,i=1,ndm),1,rad,2,ang,
     &                        (i,dr,i,di,i=ndf1,ndf2)
            endif
            count = 48000000
          endif
          if(ndm.lt.2 .or. (x(1,n).eq.0.0d0 .and. x(2,n).eq.0.0d0)) then
            radii = 0.0d0
            theta = 0.0d0
          else
            radii = sqrt(x(1,n)**2 + x(2,n)**2)
            theta = 180.0d0/pi*atan2(x(2,n),x(1,n))
          endif
          write(iow,2001) n,(x(i,n),i=1,ndm),radii,theta,
     &                    (br(i,n),bi(i,n),i=ndf1,ndf2)
          if(ior.lt.0.and.pfr) then
            write(*,2001) n,(x(i,n),i=1,ndm),radii,theta,
     &                    (br(i,n),bi(i,n),i=ndf1,ndf2)
          endif
        endif
      end do ! n

2000  format(/'  N o d a l   D i s p l a c e m e n t s',5x,
     &  'Time',0p,e18.5/31x,'Prop. Ld. (Eigenvalue)',1p,e13.5//
     &  '   Node',6(i6,a6):/(7x,6(i6,a6):))
2001  format( i7,1p,6e12.4:/(7x,1p,6e12.4:))

      end
