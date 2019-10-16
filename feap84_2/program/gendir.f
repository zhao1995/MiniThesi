c$Id:$
      subroutine gendir(x,ctr,prt,prth,err,prtz)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate real data arrays by linear interpolation

c         Corrected: I. Romero                Date: (February  3, 2000)
c                    R. Taylor                      (February 21, 2000)

c      Inputs:
c         ctr(*)   - Header type
c         prt      - Output generated data if true
c         prth     - Output title/header data if true
c         prtz     - Do not print zero entries if true
c      Outputs:
c         x(54,*)  - Generated data (used later as x(9,6,*)
c         err      - Error flag
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'dstars.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      character cd*12, ctr*(*)
      logical   prt,prth,prtz,err,errck,pinput,nozero
      integer   i,j,mct,n,ng,nn,l,lg
      real*8    xli, x(54,*),xl(3),td(5)

      save

      cd  = ctr
      mct = 0
      n   = 0
      ng  = 0

c     Input next record (save old numbers into l and lg)

100   l   = n
      lg  = ng

c     Call input routine - values returned in td and then moved

101   if(ior.lt.0) then
        write(*,3000)
        call pprint('   >')
      endif
      errck = pinput(td,5)
      if(errck) go to 101

      nn    = nint(td(1))
      n     = nn + starnd
      ng    = nint(td(2))
      xl(1) = td(3)
      xl(2) = td(4)
      xl(3) = td(5)

c     Exit if node zero or greater than numnp

      if(n.gt.numnp) then
        write(ilg,4001) n,cd
        write(iow,4001) n,cd
        if(ior.lt.0) then
          write(*,4001) n,cd
        endif
      endif
      if(nn.le.0.or.n.gt.numnp) go to 109

c     Assign data to node 'n'

      do i = 1,3
        x(i+6,n) = xl(i)
      end do ! I

c     Fill missing data

      if(lg.ne.0) then
        lg  = sign(lg,n-l)

c       compute increment in data

        xli = abs(lg)/(abs(n-l+lg)-1)
        do i = 1,3
          xl(i) = (x(i+6,n) - x(i+6,l))*xli
        end do ! I

c       Loop 'l+lg' to 'n-' to fill

106     l = l + lg
        if((n-l)*lg.le.0) go to 100

        if(l.le.0.or.l.gt.numnp) go to 108
        do i = 1,3
          x(i+6,l) = x(i+6,l-lg) + xl(i)
        end do ! I
        go to 106

c       Print error message -- then continue for checking only

108     write(ilg,4000) l,cd
        write(iow,4000) l,cd
        if(ior.lt.0) then
          write(*,4000) l,cd
        endif
        err = .true.

      endif

      go to 100

c     Output quantities

109   if(prt) then
        do j = 1,numnp
          nozero = prtz
          if(.not.prtz) then
            do l = 7,9
              if(x(l,j).ne.0.0d+0) nozero = .true.
            end do ! l
            if(nozero) then
              mct = mct - 1
              if(mct.le.0) then
                mct = 50
                call prtitl(prth)
                write(iow,2000) cd,(l,cd,l=1,3)
                if(ior.lt.0) then
                  write(*,2000) cd,(l,cd,l=1,3)
                endif
              endif
              if(mr(np(190)+j-1).ge.0) then
                write(iow,2001) j,(x(l+6,j),l=1,3)
                if(ior.lt.0) then
                  write(*,2001) j,(x(l+6,j),l=1,3)
                endif
              endif
            endif
          endif
        end do ! j
      endif ! prt

c     Formats

2000  format(5x,'Nodal',a//6x,'Node',6(i7,a6)/(10x,6(i7,a6)))

2001  format(i10,1p,6e13.4:/(10x,1p,6e13.4))

3000  format(' Input: node#, inc., values')

4000  format(' *ERROR* GENDIR: Attempt to generate node',i5,' in ',a)

4001  format(' *ERROR* GENDIR: Attempt to input node',i5,', terminate',
     &       ' input of nodes in ',a)

      end
