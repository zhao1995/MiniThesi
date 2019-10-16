c$Id:$
      subroutine genint(nd,ix,nix,num,cname,carg,prt,prth,err,type)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate integer data arrays

c      Inputs:
c         nd        - Number of integer items/node to generate
c         nix       - Dimension of integer array ix
c         num       - Number of generated items
c         cname     - Header identification for outputs
c         carg      - Value identifier for outputs
c         prt       - Output generated data if true
c         prth      - Output title/header if true
c         type      - Generation type: 1=node, 2=element

c      Outputs:
c         ix(nix,*) - Integer data generated
c         err       - True if error occurred
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'dstars.h'
      include  'iofile.h'

      logical   prt,prth,err,errck, pinput, oflg
      character cd*12, cname*(*), carg*(*)
      integer   i,j,l,n,lg,ng,nd,nn,nix,num,mct,type
      integer   ix(nix,*),ixl(14)
      real*8    td(16)

      save

      err = .false.
      cd  = carg
      mct = 0
      n   = 0
      ng  = 0
100   l   = n
      lg  = ng

c     Call input routine - values returned in td and then moved

101   if(ior.lt.0) then
        write(*,2010)
        call pprint('   >')
      endif
      errck = pinput(td,2+nd)
      if(errck) go to 101
      nn = nint(td(1))
      if(type.eq.1) then
        n = nn + starnd
      elseif(type.eq.2) then
        n = nn + starel
      else
        write(*,*) ' *ERROR* Wrong call to GENINT'
        call plstop()
      endif
      if(n.gt.num) then
        write(ilg,3001) n,cname
        write(iow,3001) n,cname
        if(ior.lt.0) write(*,3001) n,cname
      endif
      if(nn.le.0.or.n.gt.num) go to 104
      ng  = nint(td(2))
      do i = 1,nd
        ixl(i)  = nint(td(i+2))
        ix(i,n) = ixl(i)
      end do ! i
      if(lg.ne.0) then
        lg = sign(lg,n-l)
102     l = l + lg
        if((n-l)*lg.le.0) go to 100
        if(l.le.0.or.l.gt.num) go to 103
        do i = 1,nd
          ix(i,l) = ix(i,l-lg)
        end do ! i
        go to 102
103     write(ilg,3000) l,cname
        write(iow,3000) l,cname
        if(ior.lt.0) write(  *,3000) l,cname
        err = .true.
      endif
      go to 100

104   if(.not.prt) return

      do j = 1,num

c       Output if non-zero

        oflg = .false.
        do l = 1,nd
          if(ix(l,j).ne.0) oflg = .true.
        enddo ! l

        if(oflg) then

c         Output header

          mct = mct - 1
          if(mct.le.0) then
            mct = 50
            call prtitl(prth)
            if(type.eq.1) then
              write(iow,2001) cname,(l,cd,l=1,nd)
              if(ior.lt.0) then
                write(*,2001) cname,(l,cd,l=1,nd)
              endif
            else
              write(iow,2002) cname,(l,cd,l=1,nd)
              if(ior.lt.0) then
                write(*,2002) cname,(l,cd,l=1,nd)
              endif
            endif
          endif

c         Output values

          write(iow,2009) j,(ix(l,j),l=1,nd)
          if(ior.lt.0) then
            write(*,2009) j,(ix(l,j),l=1,nd)
          endif
        endif

      end do ! j

c     Formats

2001  format(5x,a//'      Node',6(i5,a6):/(10x,6(i5,a6):))
2002  format(5x,a//'   Element',6(i5,a6):/(10x,6(i5,a6):))

2009  format(i10,6i11:/(10x,6i11):)

2010  format(' Input: item#, inc., values')

3000  format(' *ERROR* GENINT: Attempt to generate item',i5,' in ',a)

3001  format(' *ERROR* GENINT: Attempt to input item',i5,', terminate',
     & ' input in ',a)

      end
