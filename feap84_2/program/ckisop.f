c$Id:$
      subroutine ckisop(ix,xl,shp,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check isoparametric elements for data input errors

c      Inputs:
c         ix(*)     - List of nodes connected to element
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh

c      Outputs:
c         None

c      Scratch:
c         shp(*)    - Storage for shape functions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   trifl
      integer   ndm, ineg, i,j,l, xn(9),yn(9),ic(2,16),ix(*)
      real*8    xsj, ss(2),shp(*),xl(ndm,*),jac(16)

      save

      data      xn/-1, 1,1,-1, 0,1,0,-1,0/
      data      yn/-1,-1,1, 1,-1,0,1, 0,0/

c     Check element for input errors

      ineg = 0
      do l = 1,nel
        if(ix(l).gt.0) then
          if(mr(np(190)+ix(l)-1).lt.0) then
            ineg       = ineg + 1
            ic(1,ineg) = l
            ic(2,ineg) = abs(ix(l))
          endif
        endif
      end do ! l
      if(ineg.gt.0) then
        write(iow,2000) n,(ic(1,i),ic(2,i),i=1,ineg)
        if(ior.lt.0) then
          write(*,2000) n,(ic(1,i),ic(2,i),i=1,ineg)
        endif
      else
        ineg = 0
        do l = 1,nel
          ss(1) = xn(l)
          ss(2) = yn(l)
          call  shp2d (ss,xl,shp,xsj,ndm,nel,ix,.false.)
          if(xsj.le.0.0d0) then
            ineg       = ineg + 1
            ic(1,ineg) = l
            ic(2,ineg) = abs(ix(l))
          endif
        end do ! l
        if(ineg.gt.0) then
          trifl = .false.
          if(ineg.eq.2 .and. (ic(2,1).eq.ic(2,2))) then
            trifl = .true.
          elseif(ineg.eq.1 .and. nel.eq.4) then
            i = ic(1,1)
            j = mod(i,4) + 1
            if(i.eq.1) then
              l = 4
            else
              l = i - 1
            endif
            if((ix(i).eq.ix(j)) .or. (ix(i).eq.ix(l))) then
              trifl = .true.
            endif
          endif
          if(.not.trifl) then
            write(iow,2001) n,(ic(1,i),ic(2,i),jac(i),i=1,ineg)
            if(ior.lt.0) then
              write(*,2001) n,(ic(1,i),ic(2,i),jac(i),i=1,ineg)
            endif
            call iprint(ix,1,nel,1,'IX (CHECk may change order)')
            call mprint(xl,2,nel,ndm,'XL (May be local projections)')
          endif
        endif

c       Try to fix element

        if(ineg.eq.nel) then
          if(nel.eq.3) then
            l     = ix(2)
            ix(2) = ix(3)
            ix(3) = l
            write(iow,2002) n,(ix(l),l=1,3)
            if(ior.lt.0) write(*,2002) n
          elseif(nel.eq.4) then
            trifl = .false.
            do i = 1,nel
              l = mod(i,4) + 1
              if(ix(i).eq.ix(l)) then
                trifl = .true.
                exit
              end if
            end do ! i
            if(.not.trifl) then
              l     = ix(1)
              ix(1) = ix(4)
              ix(4) = l
              l     = ix(2)
              ix(2) = ix(3)
              ix(3) = l
              write(iow,2002) n,(ix(l),l=1,4)
              if(ior.lt.0) write(*,2002) n
            endif
          endif
        elseif(ineg.eq.2 .and. nel.eq.4) then
          if((ic(2,1).ne.ic(2,2))) then
            l           = ix(ic(1,1))
            ix(ic(1,1)) = ix(ic(1,2))
            ix(ic(1,2)) = l
            write(iow,2002) n,(ix(l),l=1,4)
            if(ior.lt.0) write(*,2002) n
          endif
        endif
      endif

2000  format(' >Element',i9,' coordinates not input for nodes:'/
     &      ('                    Local =',i3,' Global =',i9))

2001  format(/' >Element',i9,' has zero or negative jacobian':,
     &        ' at nodes:'/
     &      ('     Local =',i3,' Global =',i9,' Jacobian =',1p,1e12.5))

2002  format(' >Element',i9,' Reverse numbers to fix negative jacobian':
     &      /'            New IX Order:',4(i10:))
      end
