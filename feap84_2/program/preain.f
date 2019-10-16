c$Id:$
      subroutine preain(id,f,ul,ndf,numnp, fname,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change sign on 'ul' added to f(*,*)              29/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Reaction to force inputs from file

c      Inputs:
c         id(*)      - Equation numbers for active dof
c         ndf        - Number dof/node
c         numnp      - Number of nodes in mesh
c         fname      - File name to use for reads
c         prt        - Output results if true
c         prth       - Output title/header data if true

c      Scratch:
c         ul(*)      - Local array of data

c      Outputs:
c         f(*)       - Nodal force/displacement
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iodata.h'
      include  'iofile.h'

      character fname*(*),fnamr*18,fext*8
      logical   exst,nonz,prt,prth
      integer   i,j,n, ndf,numnp
      integer   id(ndf,numnp)
      real*8    f(ndf,numnp),ul(*)

      save

c     Input reactions from 'name' file

      fnamr = fname
      fext  = 'reac'
      call opnfil(fext,fnamr,-2,ios,exst)
      if(exst) then

        call prtitl(prth)
        if(prt) write(iow,2000) (i,i=1,ndf)

        do i = 1,numnp

c         Input a record

          read(ios,*,end=100) n,(ul(j),j=1,ndf)

c         Check if current dof is a force: add reaction

          nonz = .false.
          do j = 1,ndf
            if(id(j,n).eq.0) then
              f(j,n) = f(j,n) + ul(j)
              nonz = .true.
            else
              ul(j) = 0.0d0
            endif
          end do ! j
          if(nonz.and.prt) write(iow,2001) n,(ul(j),j=1,ndf)

        end do ! i
100     close(ios)
        if(ior.lt.0) write(*,2002) fnamr
        write(iow,2002) fnamr

      else
        write(*,*) ' *ERROR* PREAIN: No file with reactions'
      endif

c     Formats

2000  format(/'  R e a c t i o n   N o d a l    F o r c e s'/
     &       /4x,'Node',6(i6,'-Force'):/(8x,6(i6,'-Force')))

2001  format(i8,1p,6e12.4:/(8x,1p,6e12.4))

2002  format(/'   Reaction forces loaded from file: ',a)

      end
