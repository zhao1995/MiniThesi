c$Id:$
      subroutine pcormv(x,ndm,numnp,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Move coordinates to new specified locations

c      Inputs:
c         x(ndm,*) - Nodal coordinates before move
c         ndm      - Spatial dimension of mesh
c         numnp    - Number of nodes in mesh
c         prt      - Print generated data if true
c         prth     - Print title/header data if true

c      Outputs:
c         x(ndm,*) - Nodal coordinates after moves
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   prt,prth, flg, xact(18), errck, pinput
      integer   ndm,numnp, i,i1,j,n
      real*8    xmin,xmax,tolx, x(ndm,numnp),xold(6),xnew(6),td(18)

      save

      if(ior.lt.0) then
        write(*,3000)
        call pprint('  ->')
      endif
      errck = pinput(td,3*ndm)

      if(prt) then
        call prtitl(prth)
        if(ior.lt.0) write(*,2000)
        write(iow,2000)
      end if

      do i = 1,ndm
        xact(i) = .false.
      end do ! i

      do i = 1,3*ndm,3
        i1 = nint(td(i))
        if(i1.gt.0) then
          xact(i1) = .true.
          xold(i1) = td(i+1)
          xnew(i1) = td(i+2)
          if(prt) then
            if(ior.lt.0) write(*,2001) i1,xold(i1),i1,xnew(i1)
            write(iow,2001) i1,xold(i1),i1,xnew(i1)
          end if
        end if
      end do ! i

c     Loop over nodes to reposition

      if(prt) then
        if(ior.lt.0) write(*,2002) (j,j=1,ndm)
        write(iow,2002) (j,j=1,ndm)
      end if

      xmin = x(1,1)
      xmax = x(1,1)
      do i = 1,ndm
        do n = 1,numnp
          xmin = min(xmin,x(i,n))
          xmax = max(xmin,x(i,n))
        end do ! n
      end do ! i
      tolx = 1.0d-5*(xmax-xmin)

      do i = 1,ndm
        if(xact(i)) then
          do n = 1,numnp
            if(abs(x(i,n)-xold(i)).lt.tolx) then
              flg = .true.
              do j = i+1,ndm
                if(xact(j).and.abs(x(j,n)-xold(j)).gt.tolx) then
                  flg = .false.
                end if
              end do ! j
              if(flg) then
                do j = 1,ndm
                  if(xact(j)) x(j,n) = xnew(j)
                end do ! j
                if(prt) then
                  if(ior.lt.0) write(*,2003) n,(x(j,n),j=1,ndm)
                  write(iow,2003) n,(x(j,n),j=1,ndm)
                end if
              end if
            end if
          end do ! n
        end if
      end do ! i

c     Formats

2000  format(5x,'R e p o s i t i o n   C o o r d i n a t e s'/)
2001  format(' x_old(',i1,') = ',1p,e12.5,
     &       ' x_new(',i1,') = ',1p,e12.5)
2002  format(/'   Node',6(i6,'-coord'))
2003  format(i7,6e12.5)

3000  format(' Input: dir,x_old,x_new values')

      end
