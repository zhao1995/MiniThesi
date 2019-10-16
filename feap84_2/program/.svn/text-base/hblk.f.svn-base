c$Id:$
      subroutine hblk(nr,ns,xl,ixl,shp,x,ix,dr,ds,ni,ne,ndm,nel1,
     &                nodinc,ntypg,nm,ma,xlm,ctype,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate a block of (shell) elements

c      Programmed: M.S.Rifai/JCS for generation
c                  of top-surface shell mesh

c     ntypg = 20   4-node quadrilaterals
c     ntypg = 21   3-node triangles - diags ll to ur
c     ntypg = 22   3-node triangles - diags ul to lr
c     ntypg = 23   3-node triangles - diags
c     ntypg = 24   3-node triangles - diags
c     ntypg = 25   3-node triangles - diags union jack
c     ntypg = 26   3-node triangles - diags union jack
c     ntypg = 27   6-node triangles - diags ll to ur
c     ntypg = 28   8-node quadrilaterals
c     ntypg = 29   9-node quadrilaterals


c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat2.h'
      include  'iofile.h'
      include  'crotas.h'
      include  'trdata.h'

      character ctype*15, xh*6
      logical   prt,prth,ityp, pcomp
      integer   i, j, k, mct, ntyp, nm, ni, ndm, me, nodinc
      integer   nr, ns, ne, nel1, ntypg, ma, n, inc
      integer   ixl(*),ix(nel1,*)
      real*8    dr, ds, xsj, rr, sind2,cosd2, sind3,cosd3
      real*8    ss(2),xl(3,*),x(ndm,*),shp(3,*),xlm(9,6,*),xx(3)

      save

      data      xh/' coord'/

c     Set flags and options

      ntyp = ntypg - 20
      nm   = nm / 2

c     Generate nodes

      do k = 1,3
        xx(k) = 0.0d0
      end do ! k
      n     =  ni
      mct   =  0
      ss(2) = -1.0d0
      do j = 1,ns
        ss(1) = -1.0d0
        do i = 1,nr
          call shp2d(ss,xl,shp,xsj,3,nm,ixl,.true.)
          call cfunc(shp,xl,ixl,nm,ndm,xx)

c         Convert coordinates if necessary

          if(pcomp(ctype,'pola',4) .or. pcomp(ctype,'cyli',4)) then
            call pdegree(xx(2), sind2,cosd2)
            rr    = xx(1)
            xx(1) = x0(1) + rr*cosd2
            xx(2) = x0(2) + rr*sind2
          elseif(ndm.ge.3 .and. pcomp(ctype,'sphe',4)) then
            call pdegree(xx(2), sind2,cosd2)
            call pdegree(xx(3), sind3,cosd3)
            rr    = xx(1)
            xx(1) = x0(1) + rr*cosd2*sind3
            xx(2) = x0(2) + rr*sind2*sind3
            xx(3) = x0(3) + rr*cosd3
          endif

c         Transform to global coordinates

          do k = 1,ndm
            x(k,n) = xr(k) + tr(k,1)*xx(1)
     &                     + tr(k,2)*xx(2)
     &                     + tr(k,3)*xx(3)
          end do ! k

c         Output coordinates

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
          n     = n + 1
          ss(1) = ss(1) + dr
        end do ! i
        n     = n + nodinc
        ss(2) = ss(2) + ds
      end do ! j

c     Generate transformation matrices

      if (frotas) then
        do k = 1,3
          xx(k) = 0.0d0
        end do ! k
        n     =  ni
        mct   =  0
        ss(2) = -1.0
        do j = 1,ns
          ss(1) = -1.0
          do i = 1,nr
            call shp2d(ss,xl,shp,xsj,3,nm,ixl,.true.)
            call cfunc(shp,xl(1,nm+1),ixl,nm,ndm,xx)

c           Convert coordinates if necessary

            if(pcomp(ctype,'pola',4) .or. pcomp(ctype,'cyli',4)) then
              call pdegree(xx(2), sind2,cosd2)
              rr    = xx(1)
              xx(1) = x0(1) + rr*cosd2
              xx(2) = x0(2) + rr*sind2
            elseif(ndm.ge.3 .and. pcomp(ctype,'sphe',4)) then
              call pdegree(xx(2), sind2,cosd2)
              call pdegree(xx(3), sind3,cosd3)
              rr    = xx(1)
              xx(1) = x0(1) + rr*cosd2*sind3
              xx(2) = x0(2) + rr*sind2*sind3
              xx(3) = x0(3) + rr*cosd3
            endif

c           Transform to global coordinates

            do k = 1,ndm
              xlm(k+6,1,n) = xr(k) + tr(k,1)*xx(1)
     &                             + tr(k,2)*xx(2)
     &                             + tr(k,3)*xx(3)
            end do ! k

            n     = n + 1
            ss(1) = ss(1) + dr
          end do ! i
          n     = n + nodinc
          ss(2) = ss(2) + ds
        end do ! j
      endif ! frotas

c     Generate element connectivity

      if(ne.gt.0) then
        me  = ne - 1
        inc = 1
        if(ntyp.ge.7) inc = 2
        do j = 1,ns-1,inc
          n = (nr+nodinc)*(j-1) + ni
          do i = 1,nr-1,inc
            n           = n + 1
            me          = me + 1
            ix(nel1,me) = ma
            if(ntyp.eq.0) then
              ix(1,me)    = n - 1
              ix(2,me)    = n
              if(ndm.ne.1) then
                ix(3,me)  = n + nr + nodinc
                ix(4,me)  = n + nr - 1 + nodinc
              endif
            elseif(ntyp.eq.7) then
              ix(1,me)    = n-1
              ix(4,me)    = n
              ix(2,me)    = n+1
              ix(6,me)    = nr+nodinc + n
              ix(5,me)    = nr+nodinc + n+1
              ix(3,me)    = 2*(nr+nodinc) + n+1
              me          = me + 1
              ix(1,me)    = n-1
              ix(6,me)    = nr+nodinc + n-1
              ix(4,me)    = nr+nodinc + n
              ix(3,me)    = 2*(nr+nodinc) + n-1
              ix(5,me)    = 2*(nr+nodinc) + n
              ix(2,me)    = 2*(nr+nodinc) + n+1
              ix(nel1,me) = ma
              n           = n+1
            elseif(ntyp.ge.8) then
              ix(1,me) = n-1
              ix(5,me) = n
              ix(2,me) = n+1
              ix(8,me) = nr+nodinc + n-1
              if(ntyp.gt.8) ix(9,me) = nr+nodinc + n
              ix(6,me) = nr+nodinc + n+1
              ix(4,me) = 2*(nr+nodinc) + n-1
              ix(7,me) = 2*(nr+nodinc) + n
              ix(3,me) = 2*(nr+nodinc) + n+1
              n        = n+1
            else
              ityp = (ntyp.eq.1) .or. (ntyp.eq.3.and.mod(j,2).eq.1)
     &          .or. (ntyp.eq.4.and.mod(j,2).eq.0) .or. (ntyp.eq.5.and.
     &          mod(i+j,2).eq.0) .or. (ntyp.eq.6.and.mod(i+j,2).eq.1)
              if(ityp) then
                ix(1,me)    = n-1
                ix(2,me)    = n + nr + nodinc
                ix(3,me)    = n + nr + nodinc - 1
                me          = me + 1
                ix(1,me)    = n-1
                ix(2,me)    = n
                ix(3,me)    = n + nr + nodinc
                ix(nel1,me) = ma
              else
                ix(1,me)    = n-1
                ix(2,me)    = n
                ix(3,me)    = n + nr + nodinc - 1
                me          = me + 1
                ix(1,me)    = n
                ix(2,me)    = n + nr + nodinc
                ix(3,me)    = n + nr + nodinc - 1
                ix(nel1,me) = ma
              endif
            endif
          end do ! i
        end do ! j
      endif

c     Formats

2003  format('  N o d a l   C o o r d i n a t e s'//6x,'Node',5(i7,a6))
2004  format(i10,1p,5e13.4)

      end
