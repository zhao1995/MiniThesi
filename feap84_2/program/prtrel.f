c$Id:$
      subroutine prtrel(r,x,u,ttim,ndm,ndf,nlist,list,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'dabs' to 'abs'                           17/11/2006
c       2. Add coords and compute moment sums               12/04/2013
c          Add displacements to compute current moments.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output nodal reactions based on list

c      Inputs:
c         r(*)      - Current value of reactions
c         x(*)      - Nodal coordinates
c         u(*)      - Nodal displacements
c         ttim      - Value of solution time
c         ndf       - Number dof/node
c         nlist     - Number items in list
c         list(*)   - List of node numbers to output
c         prth      - Output title/header data if true

c      Outputs:
c         None      - Outputs to file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'

      logical   prth
      integer   ndm,ndf,nlist, i,k,n,nn, count, list(*)
      real*8    ttim
      real*8    x(ndm,*), u(ndf,*), r(ndf,*)
      real*8    rsum(30),asum(30),psum(30),msum(3),usum(3)

      save

      do k = 1,ndf
        psum(k) = 0.d0
        rsum(k) = 0.d0
        asum(k) = 0.d0
      end do ! k
      msum = 0.0d0
      usum = 0.0d0
      do i = 1,numnp
        do k = 1,ndf
          rsum(k) = rsum(k) - r(k,i)
          asum(k) = asum(k) + abs(r(k,i))
        end do ! k
      end do ! i

      count = 0

      do nn = 1,nlist
        n = list(nn)
        count = count - 1
        do k = 1,ndf
          psum(k) = psum(k) - r(k,n)
        end do ! k
        if(ndm.eq.2 .and. ndf.ge.2) then
          msum(3) = msum(3) + x(2,n)*r(1,n) - x(1,n)*r(2,n)
          usum(3) = usum(3) + u(2,n)*r(1,n) - u(1,n)*r(2,n)
        elseif(ndm.eq.3 .and. ndf.ge.3) then
          msum(1) = msum(1) + x(3,n)*r(2,n) - x(2,n)*r(3,n)
          msum(2) = msum(2) + x(1,n)*r(3,n) - x(3,n)*r(1,n)
          msum(3) = msum(3) + x(2,n)*r(1,n) - x(1,n)*r(2,n)
          usum(1) = usum(1) + u(3,n)*r(2,n) - u(2,n)*r(3,n)
          usum(2) = usum(2) + u(1,n)*r(3,n) - u(3,n)*r(1,n)
          usum(3) = usum(3) + u(2,n)*r(1,n) - u(1,n)*r(2,n)
        endif
        if(count.le.0) then
          call prtitl(prth)
          write(iow,2000) ttim,(k,k=1,ndf)
          if(ior.lt.0) then
            write(*,2000) ttim,(k,k=1,ndf)
          endif
          count = 50000000
        endif
        if(ior.lt.0) then
          write(*,2001) n,(-r(k,n),k=1,ndf)
        endif
        write(iow,2001) n,(-r(k,n),k=1,ndf)
      end do ! nn

c     Print sum checks

      write(iow,2002) ' Pr.sum',(psum(k),k=1,ndf)
      write(iow,2003) '   Sum ',(rsum(k),k=1,ndf)
      write(iow,2003) '  |Sum|',(asum(k),k=1,ndf)
      write(iow,2003) '  M_ref',(msum(k),k=1,3)
      write(iow,2003) '  M_cur',(msum(k)+usum(k),k=1,3)
      if(ior.lt.0) then
        write(*,2002) ' Pr.sum',(psum(k),k=1,ndf)
        write(*,2003) '   Sum ',(rsum(k),k=1,ndf)
        write(*,2003) '  |Sum|',(asum(k),k=1,ndf)
        write(*,2003) '  M_ref',(msum(k),k=1,3)
        write(*,2003) '  M_cur',(msum(k)+usum(k),k=1,3)
      endif

c     Formats

2000  format('  N o d a l    R e a c t i o n s',12x,
     &  'Time',e18.5//'   Node',6(i8,' dof'):/(7x,6(i8,' dof')))

2001  format(i7,1p,6e12.4:/(7x,1p,6e12.4:))

2002  format(/a,1p,6e12.4:/(7x,1p,6e12.4:))

2003  format( a,1p,6e12.4:/(7x,1p,6e12.4:))

      end
