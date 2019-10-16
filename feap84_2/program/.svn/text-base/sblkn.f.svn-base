c$Id:$
      subroutine sblkn(nr,ns,xl,ixl,shp,x,dr,ds,ni,n,ndm,nodinc,ntyp,
     &                 nm,ctype,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'cap' option for spherical cap               11/11/2008
c       2. Change ss to ss(1) on bernshp                    09/01/2012
c       3. Add 'sphe'rical option for surfaces              06/06/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate nodes and elements for 2-d problems

c         nm < 4 : Generates a line of elements
c           ns = 1     2-node line
c           ns = 2     3-node line
c         nm > 3 : Generates a block of elements
c           ntyp = 1   4-node quadrilaterals
c           ntyp = 2   3-node triangles - diags ll to ur
c           ntyp = 3   3-node triangles - diags ul to lr
c           ntyp = 4   3-node triangles - diags
c           ntyp = 5   3-node triangles - diags
c           ntyp = 6   3-node triangles - diags union jack
c           ntyp = 7   6-node triangles - diags ll to ur
c           ntyp = 8   8-node quadrilaterals
c           ntyp = 9   9-node quadrilaterals

c           ntyp =-1   3-node triangles - crossed pattern
c           ntyp =-7   7-node triangles - diags ll to ur (bubble node)

c      Inputs:
c         nr        - Number nodes in 1-local coordinate dir.
c         ns        - Number nodes in 2-local coordinate dir.
c         xl(ndm,*) - Block nodal coordinate array
c         ixl(*)    - Block nodal connection list
c         shp(3,*)  - Shape functions for block
c         dr        - 1-local coordinate increment
c         ds        - 2-local coordinate increment
c         ni        - Initial node number for block
c         ndm       - Spatial dimension of mesh
c         nodinc    - Increment to node numbers at end of each line
c         ntyp      - Block type
c         nm        - Number of nodes on block
c         ctype     - Type of block coordinates
c         prt       - Output generated data if true
c         prth      - Output title/header data if true

c      Outputs:
c         n         - Final node number for block
c         x(ndm,*)  - Nodal coordinates for block
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat2.h'
      include  'iofile.h'
      include  'trdata.h'

      logical   prt,prth, pcomp, bernfl
      character xh*6, ctype*15
      integer   i,j,k,n, mct,me, nm, inc8
      integer   nr,ns,nsn,ni,ndm,nodinc,ntyp, ixl(*)
      real*8    rr,sn,cn, sn2,cn2, sn3,cn3,  drh,dsh,dr,ds, xsj
      real*8    xl(3,*),x(ndm,*),shp(3,*),shp1(2,3),xx(3),ss(2),tt(2)

      save

      data      xh /' coord'/

c     Set flag for Bernstein generations

      bernfl = pcomp(ctype,'bern',4)

c     Initialize data

      do k = 1,3
        xx(k) = 0.0d0
      end do ! k

      drh =  dr*0.5d0
      dsh =  ds*0.5d0
      n   =  ni
      mct =  0
      if(nm.lt.4 .or. bernfl) then
        nsn = 1
      else
        nsn = ns
      endif
      ss(2)   = -1.0
      do j = 1,nsn
        ss(1) = -1.0
        if(ntyp.eq.8 .and. mod(j,2).eq.0) then
          inc8 = 2
        else
          inc8 = 1
        endif
        do i = 1,nr,inc8
          if(bernfl) then
            call bernshp(ss(1),nm, xl,ndm, xx)
          elseif (nm.gt.3) then
            call shp2d(ss,xl,shp,xsj,3,nm,ixl,.true.)
            call cfunc(shp,xl,ixl,nm,ndm,xx)
          else
            call shap1d(ss(1),nm,shp1)
            do k = 1,ndm
              xx(k) = 0.0d0
              do me = 1,nm
                xx(k) = xx(k) + shp1(2,me)*xl(k,me)
              end do ! me
            end do ! k
          endif
          if(ndm.ge.2 .and.
     &      (pcomp(ctype,'pola',4) .or. pcomp(ctype,'cyli',4))) then
            call pdegree(xx(2), sn,cn)
            rr    = xx(1)
            xx(1) = x0(1) + rr*cn
            xx(2) = x0(2) + rr*sn
          elseif(ndm.eq.3 .and. pcomp(ctype,'sphe',4)) then
            call pdegree(xx(2), sn2,cn2)
            call pdegree(xx(3), sn3,cn3)
            rr    = xx(1)
            xx(1) = x0(1) + rr*cn2*sn3
            xx(2) = x0(2) + rr*sn2*sn3
            xx(3) = x0(3) + rr*cn3
          elseif(ndm.eq.3 .and. pcomp(ctype,'cap',3)) then
            xx(3) = sqrt(xx(3)**2 - xx(1)**2 - xx(2)**2)
          endif
          do k = 1,ndm
            x(k,n) = xr(k)+tr(k,1)*xx(1)+tr(k,2)*xx(2)+tr(k,3)*xx(3)
          end do ! k
          if(prt) then
            mct = mct + 1
            if(mod(mct,50).eq.1) then
              call prtitl(prth)
              write(iow,2003) (k,xh,k=1,ndm)
              if(ior.lt.0) then
                write(*,2003) (k,xh,k=1,ndm)
              endif
            endif
            write(iow,2004) n,(x(k,n),k=1,ndm)
            if(ior.lt.0) then
              write(*,2004) n,(x(k,n),k=1,ndm)
            endif
          endif
          n = n + 1
          if(nm.gt.3.and.ntyp.eq.-1.and.i.lt.nr.and.j.lt.ns) then
            tt(1) = ss(1) + drh
            tt(2) = ss(2) + dsh
            call shp2d(tt,xl,shp,xsj,3,nm,ixl,.true.)
            call cfunc(shp,xl,ixl,nm,ndm,x(1,n))
            if(ndm.ge.2 .and.
     &        (pcomp(ctype,'pola',4) .or. pcomp(ctype,'cyli',4))) then
              call pdegree(x(2,n), sn,cn)
              rr     = x(1,n)
              x(1,n) = x0(1) + rr*cn
              x(2,n) = x0(2) + rr*sn
            endif
            if(prt) then
              mct = mct + 1
              if(mod(mct,50).eq.1) then
                call prtitl(prth)
                write(iow,2003) (k,xh,k=1,ndm)
              endif
              write(iow,2004) n,(x(k,n),k=1,ndm)
              if(ior.lt.0) then
                if(mod(mct,50).eq.1) then
                  write(*,2003) (k,xh,k=1,ndm)
                endif
                write(*,2004) n,(x(k,n),k=1,ndm)
              endif
            endif
            n = n + 1
          endif
          ss(1) = ss(1) + dr*inc8
        end do ! i
        n = n + nodinc
        ss(2) = ss(2) + ds
      end do ! j

c     Formats

2003  format(/'  N o d a l   C o o r d i n a t e s'//6x,'Node',5(i7,a6))

2004  format(i10,1p,5e13.4)

      end
