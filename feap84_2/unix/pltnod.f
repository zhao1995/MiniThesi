c$Id:$
      subroutine pltnod(x,ip,ndm,numnp,n1,n2,n3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot location of nodal points in mesh

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         ip(*)     - Plot only nodes with positive values
c         ndm       - Dimension of x array
c         numnp     - Number of nodes in mesh
c         n1        - Place number near node if .ne. 0
c         n2        - First node to plot
c         n3        - Last node to plot

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'idptr.h'
      include  'pdata1.h'
      include  'pdata4.h'
      include  'pdatay.h'
      include  'plflag.h'
      include  'ppers.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   zoom
      integer   n,ne, nsy, ndm, numnp, n1, n2, n3, ip(*)
      real*8    x1, x2, x3, dx1, shft1, shft2, shft3, x(ndm,*)

      save

c     Open plot and plot locations of all nodes: Add labels if n1 .ne. 0

      dx1 = .002d0/scale
      if(kpers.ne.0) then
        shft1 =   0.d0*dx1 ! -12
        shft2 =   0.d0*dx1 !   2
        shft3 =   0.d0*dx1 ! -12
      else
        shft1 = -12.d0*dx1
        shft2 =   2.d0*dx1
        shft3 =   0.d0*dx1
      endif
      x3 = 0.0d0
      do nsy = 1,nsym
        lsym = isym(nsy)
        call pltsym(x,ndm,numnp,lsym)
        do n = n2,n3
          if(n1.lt.0) then
            ne = mr(nprn+n-1)
          else
            ne = n
          endif
          if(ip(ne).gt.0 .and. zoom(x(1,ne),ndm)) then
            x1 = x(1,ne)
            x2 = x(2,ne)
            if(ndm.ge.3) x3 = x(3,ne)
            call plotl(x1-dx1 , x2+dx1 , x3, 3)
            call plotl(x1-dx1 , x2-dx1 , x3, 2)
            call plotl(x1+dx1 , x2-dx1 , x3, 2)
            call plotl(x1+dx1 , x2+dx1 , x3, 2)
            call plotl(x1-dx1 , x2+dx1 , x3, 2)
            if(n1.ne.0) then
              call plotl(x1+shft1, x2+shft2, x3+shft3, 3)
              if(clip) call plabl(n)
            endif
          endif
        end do ! n
        call pltsym(x,ndm,numnp,lsym)
      end do ! nsy

      end
