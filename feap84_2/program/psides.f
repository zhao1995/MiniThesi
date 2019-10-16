c$Id:$
      subroutine psides(is,side,nside,prt,prth,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'nurb' side type                             03/12/2008
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Input sides for blending functions

c     Input:
c       nside   - Total number of sides
c       prt     - Output if true
c       prth    - Output title if true
c       isw     - 1=sides; 2=faces

c     Output:
c       is(*,*) - List of side nodes/face sides
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblend.h'
      include  'iofile.h'

      logical   prt,prth, pcomp, setvar,tinput
      character type*15
      integer   i,n,side, isw,nside,is(nside,*)
      real*8    td(16)

      save

      if(prt) then
        call prtitl(prth)
        if(isw.eq.1) write(iow,2000)
        if(isw.eq.2) write(iow,2001)
        if(ior.lt.0) then
          if(isw.eq.1) write(*,2000)
          if(isw.eq.2) write(*,2001)
        endif
      endif

c     Input Sides for blending inputs

      type = 'initial'
      do while(.not.pcomp(type,'    ',4))

        if(ior.lt.0) then
          if(isw.eq.1) then
            write(*,3000)
          elseif(isw.eq.2) then
            write(*,3001)
          endif
          call pprint('   >')
        endif

        setvar = tinput(type,1,td,nside-1)

        if(.not. pcomp(type,'    ',4)) then
          side   = side + 1
          if(pcomp(type,'pola',4)) then
           is(1,side) = 1
          elseif(pcomp(type,'cart',4) .or. pcomp(type,'lagr',4)) then
           is(1,side) = 0
          elseif(pcomp(type,'segm',4)) then
           is(1,side) = 2
          elseif(pcomp(type,'elip',4)) then
           is(1,side) = 3
          elseif(pcomp(type,'bern',4)) then
           is(1,side) = 4
          elseif(pcomp(type,'nurb',4)) then
           is(1,side) = 5
          else
           is(1,side) = -1
          endif

          if(is(1,side).ge.0) then
            do n = 2,nside
              is(n,side) = nint(td(n-1))
            end do ! n
            do n = nside,1,-1
              if(is(n,side).ne.0) go to 100
            end do ! n

100         if(prt) then
              write(iow,2002) side,type,(is(i,side),i=2,n)
              if(ior.lt.0) then
                write(*,2002) side,type,(is(i,side),i=2,n)
              endif
            endif
          else
            side = side - 1
          endif

        endif

      end do ! while

      numsd = max(numsd,side)

3000  format('   Input: side, is(side,snode),snode=2,nside')

3001  format('   Input: face, is(face,side),side=1,nface')

2000  format('   S i d e   S u p e r N o d e   C o n n e c t i o n s'//
     &'       Side Type 1-Nd 2-Nd 3-Nd 4-Nd 5-Nd 6-Nd 7-Nd 8-Nd 9-Nd'/)

2001  format('   F a c e    S i d e    C o n n e c t i o n s'//
     &'       Face 1-Sd 2-Sd 3-Sd 4-Sd 5-Sd 6-Sd 7-Sd 8-Sd 9-Sd'/)

2002  format(5x,i4,1x,a4,1x,9i5:/(15x,9i5))

      end
