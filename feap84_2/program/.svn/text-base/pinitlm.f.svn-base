c$Id:$
      subroutine pinitlm(u,ud,ix,ip,nen,nen1,ndf,numel,numnp,
     &                   prt,prth,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change dimension on 'ix' from 'numnp' to 'numel' 17/04/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input constant initial condition for material set
c               [init,mate] - set initial velocity for material sets
c                  or
c               [init,regi] - set initial velocity for regions

c               Data records: <disp,rate> nn (u(i),i=1,ndf)
c                 - nn   = material or region number
c                 - u(i) = values to set on each dof
c                 - Terminate on blank record

c      Inputs:
c        ix(nen1,numel)  - Element connection array
c        nen             - Maximum nodes/element
c        nen1            - First dimension of 'ix'
c        ndf             - Degree of freedoms/node
c        numel           - Number of elements in mesh
c        numnp           - Number of nodes in mesh
c        prt             - Output if true
c        prth            - Header print flag
c        isw             - Switch: 1 = material, 2 = regions

c      Working array:
c        ip(numnp)       - Store active nodes

c      Outputs:
c        u(ndf,*)        - Initial displacements
c        ud(ndf,numnp,*) - Initial displacements
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'
      include   'prflag.h'

      character  text*15, type(2)*13
      logical    prt,prth, errck,tinput, pcomp
      integer    nen,nen1,ndf,numel,numnp,isw, nc,nn, n,i,j
      integer    ix(nen1,numel),ip(numnp)
      real*8     u(ndf,*), ud(ndf,numnp,*),td(15)

      save

      data       type /' Velocity    ',' Displacement' /

c     Check isw value

      if(isw.eq.1) then  ! based on material set
        nc = nen1
      elseif(isw.eq.2) then ! based on region number
        nc = nen1 - 1
      else
        write(  *,*) '  INITIAL condition error: isw = ',isw
        write(iow,*) '  INITIAL condition error: isw = ',isw
        write(ilg,*) '  INITIAL condition error: isw = ',isw
        call plstop()
      endif

c     Do inputs

      text = 'start'
      do while (.not.pcomp(text,'    ',4))
        errck = tinput(text,1,td,15)
        nn = nint(td(1))

c       Set velocity condition

        if(pcomp(text,'rate',4)) then
          ivelfl = .true.
          do n = 1,numnp
            if(ix(nc,n).eq.nn) then
              do i = 1,nen
                ip(ix(i,n)) = 1
                if(ix(i,n).gt.0) then
                  do j = 1,ndf
                    ud(j,ix(i,n),1) = td(j+1)
                  end do ! j
                endif
              end do ! i
            endif
          end do ! n

c       Set displacement condition

        elseif(pcomp(text,'disp',4)) then
          idisfl = .true.
          do n = 1,numnp
            if(ix(nc,n).eq.nn) then
              do i = 1,nen
                ip(ix(i,n)) = 2
                if(ix(i,n).gt.0) then
                  do j = 1,ndf
                    u(j,ix(i,n)) = td(j+1)
                  end do ! j
                endif
              end do ! i
            endif
          end do ! n
        endif
      end do ! while

c     Output values input

      if(prt) then
        call prtitl(prth)
        write(iow,2000) type(isw),(j,type(isw),j=1,ndf)
        do n = 1,numnp
          if(ip(n).eq.1) then
            write(iow,2001) n,(ud(j,n,1),j=1,ndf)
          elseif(ip(n).eq.2) then
            write(iow,2001) n,(u(j,n),j=1,ndf)
          endif
        end do ! n
      endif

c     Formats

2000  format(5x,'Nodal',a//6x,'node',6(i5,a6):/(10x,6(i5,a6)))

2001  format(i10,1p,6e11.3:/(10x,1p,6e11.3))

      end
