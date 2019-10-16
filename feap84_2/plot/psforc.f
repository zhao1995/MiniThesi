c$Id:$
      subroutine psforc(f,nty,x,idir,val,ndf,ndm,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set specificed values of force components using mouse

c      Inputs:
c         nty(*)    - Node type
c         x(ndm,*)  - Nodal coordinates in deformed state
c         idir      - Number of boundary condition component to set
c         val       - Value of component set
c         ndf       - Number dof/node
c         ndm       - Dimension of x array
c         numnp     - Number of nodes in mesh

c      Outputs:
c         f(ndf,*)  - Modified force values
c                     Results saved in file finp.for
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'pdata1.h'

      logical   noerr, exst,outf
      character button*1, fnamr*132, fext*4
      integer   idir, n,ndf,ndm,numnp,nz
      integer   nty(numnp)
      real*8    val, x1,y1,x3,dx1, xm,ym, f(ndf,numnp), x(ndm,numnp)

      save

c     Set default increment

      x3  =  0.00d0
      dx1 = .006d0/scale

c     Pick point from screen

      write(*,2001) val,idir,idir

      button = 'l'
100   call gin(x1,y1,noerr,button)

      if(noerr .and. button.ne.'m') then

        x1 = 0.5d0*(sx(1) + (x1 - s0(1))/scale)
        y1 = 0.5d0*(sx(2) + (y1 - s0(2))/scale)

c       Find closest node

        xm = abs(x(1,1) - x1)**2 + (x(2,1) - y1)**2
        nz = 1
        do n = 2,numnp
          ym = (x(1,n)-x1)**2+(x(2,n)-y1)**2
          if(ym.lt.xm) then
            xm = ym
            nz = n
          endif
        end do ! n

        if(button.eq.'l' .or. button.eq.' ') then
          f(idir,nz) = val
          call pppcol(1,1)
        elseif(button.eq.'r') then
          f(idir,nz) = 0.0d0
          call pppcol(0,1)
        endif

c       Plot boundary forces as slash (lines = val)

        if(ndm.ge.3) x3 = x(3,nz)
        call plotl(x(1,nz)+dx1, x(2,nz)-dx1, x3, 3)
        call plotl(x(1,nz)-dx1, x(2,nz)+dx1, x3, 2)
        call plclos
        call plopen

        go to 100

      endif

c     Put boundary conditions in file called 'finp.frc'

      fnamr =  finp
      fext  =  'frc'
      call addext(fnamr,fext,128,3)
      call opnfil(fext,fnamr,-1,ios,exst)
      rewind ios

      write(ios,2003)
      do nz = 1,numnp
        if(nty(nz).ge.0) then
          outf = .false.
          do n = 1,ndf
            if(f(n,nz).ne.0.0d0) outf = .true.
          end do ! n
          if(outf) write(ios,2002) nz,(f(n,nz),n=1,ndf)
          end if
      end do ! nz

      write(ios,2004)

      close(ios)

c     Formats

2001  format('  ->Place cursor over NODE:'/
     &       '    Forced value =',1p,e13.5/
     &       '    Add a  force to',i2,' dof with LEFT  button'/
     &       '    Delete force to',i2,' dof with RIGHT button'/
     &       '    EXIT with MIDDLE button (or SHIFT+button).')

2002  format(i6,' 0 ',1p,6E13.5:/(1p,8e13.5))
2003  format('forced conditions')
2004  format(' '/'incl,end')

      end
