c$Id:$
      subroutine gapnd(d,ul,xl,ix,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set plot type for element to 0 (no plot)         21/08/2008
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Penalty or Lagrange multiplier form on position
c               constraint.

c     Penalty - Augmented: 3-node element

c            o--------------o    -----> +idir
c            1              2

c          Node 1: Body 1 connection
c          Node 2: Body 2 connection

c     Lagrange Multiplier: 3-node element

c            o-------o------o    -----> +idir
c            1       3      2

c          Node 1: Body 1 connection
c          Node 2: Body 2 connection
c          Node 3: Lagrange multiplier
c-----[--.----+----.----+----.-----------------------------------------]
c     Input data: (isw = 1)

c       Record ('dire',D)
c         D - Direction of contact

c       Record ('pena',K)
c         K - Penalty parameter

c     Outputs: (isw = 4)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'comblk.h'

      logical   errck, tinput, pcomp, elast, doflg, conflg
      character type*15
      integer   i, idir, idof,isgn, j, ndf, ndm, nst, isw, ix(*)
      real*8    dg,fr,fc,td(1)
      real*8    d(*),ul(ndf,*),xl(ndm,*),s(nst,*),r(ndf,*)

      save

c     INPUT MATERIAL PROPERTIES

      if(isw.eq.1) then

c       Record: 'dire',D

        elast = .true.
        doflg = .true.

        type = 'start'
        do while(.not.pcomp(type,'    ',4))

          if(ior.lt.0) then
            write(*,3000)
            call pprint('     >')
          endif
          errck = tinput(type,1,td,1)

c         Directions

          if    (pcomp(type,'dire',4)) then

            d(1)  = max(-ndm,min(ndm,int(td(1))))
            elast = .false.

c         Directions

          elseif(pcomp(type,'degr',4)) then

            d(2)  = max(1,min(ndf,int(td(1))))
            doflg = .false.

c         Penalty parameter

          elseif(pcomp(type,'pena',4)) then

            d(3) = td(1)

          end if

        end do ! while

c       Output material properties

        if(elast) then
          write(iow,4002) ma
          if(ior.lt.0) then
            write(*,4002) ma
            return
          end if
          call plstop()
        end if

c       Check if dof set: Default is same as |direction|

        if(doflg) then
          d(2) = abs(nint(d(1)))
        endif

        write(iow,2000) nint(d(1)),nint(d(2)),d(3)
        if(ior.lt.0) then
          write(*,2000) nint(d(1)),nint(d(2)),d(3)
        end if

c       Set history terms

        nh1 = 1
        nh3 = 2

c       Deactivate all dof in element

        do i = 1,ndf
          ix(i) = 0
        end do ! i

c       Activate gap degree-of-freedom only

        i     = max(1,min(ndf,nint(d(2))))
        ix(i) = 1

c       Plot type

        pstyp = 0  ! No plot

c     CHECK ELEMENTS

      elseif(isw.eq.2) then

        if(ix(1).eq.0 .or. ix(2).eq.0) then
          write(iow,4000) n,ix(1),ix(2)
          if(ior.lt.0) then
            write(*,4000) n,ix(1),ix(2)
          endif
        elseif(nel.gt.2 .and. ix(3).lt.max(ix(1),ix(2))) then
          write(iow,4001) n,ix(3),ix(1),ix(2)
          if(ior.lt.0) then
            write(*,4001) n,ix(3),ix(1),ix(2)
          endif
        endif

c     COMPUTE ELEMENT STIFFNESS AND RESIDUAL ARRAYS

      elseif(isw.eq.3  .or. isw.eq.4 .or.
     &       isw.eq.6  .or. isw.eq.10 ) then

        isgn   = nint(d(1))
        idir   = abs(isgn)
        idof   = nint(d(2))
        fc     = 0.0d0
        fr     = 0.0d0
        i      = ndf + idof
        conflg = .false.

c       Compute the gap function

        if(idir.eq.idof) then
          dg  = (ul(idof,2) - ul(idof,1)) + (xl(idir,2) - xl(idir,1))
        else
          dg  =  ul(idof,2) - ul(idof,1)
        endif

c       Reverse for negative directions

        if(isgn.lt.0) then
          dg  = -dg
        endif

c       Lagrange multiplier formulation

        if(nel.gt.2) then

          j      = ndf + i

c         Extract contact force

          fc = ul(idof,3)

c         Compute contact state

          if(fc.gt.0.0d0 .or. dg.lt.0.0d0) then
            if(isw .eq. 3) then
              s(idof,j  ) = -ctan(1)
              s(j  ,idof) = -ctan(1)
              s(i  ,j  )  =  ctan(1)
              s(j  ,i  )  =  ctan(1)
              s(j  ,j  )  =  0.d0
            endif

c           Contact element residual

            if(isgn.gt.0) then
              r(idof,1) = -fc
              r(idof,2) =  fc
              r(idof,3) = -dg
            else
              r(idof,1) =  fc
              r(idof,2) = -fc
              r(idof,3) =  dg
            endif

c           Set flag for regularized term

            conflg = .true.
            fr     =  d(3)*dg

c         Release state

          else

            s(j,j)    =  ctan(1)
            r(idof,3) = -ul(idof,3)
            fc        =  0.0d0

          endif

c       Penalty/Augmented formulation (nel.eq.2)

        else

c         Extract contact force

          fc =  hr(nh2) + augf*d(3)*dg
          if(fc .lt. 0.0d0) then
            conflg = .true.
            fr     =  fc
          endif

        endif

c       Contact state

        if(conflg) then
          s(idof,idof) =  d(3)*ctan(1)
          s(idof,i   ) = -d(3)*ctan(1)
          s(i   ,idof) = -d(3)*ctan(1)
          s(i   ,i   ) =  d(3)*ctan(1)
          if(isgn.gt.0) then
            r(idof,1) = r(idof,1) + fr
            r(idof,2) = r(idof,2) - fr
          else
            r(idof,1) = r(idof,1) - fr
            r(idof,2) = r(idof,2) + fr
          endif

c       Release state

        else
          hr(nh2) = 0.0d0
          fc      = 0.0d0
          fr      = 0.0d0
        endif

c       Tplot state

        tt(1) = fc
        tt(2) = dg

c       Output contact force/element

        if(isw.eq.4) then
          write(iow,2001) n,fc,dg
          if(ior.lt.0) then
            write(*,2001) n,fc,dg
          endif
        endif

c       Augmented update

        if(isw.eq.10 .and. nel.eq.2) then
          hr(nh2) =  fr
        endif

c     Initialize element strains for activation

      elseif(isw.eq.17) then

        idof      = nint(d(2))
        hr(nh3  ) = ul(idof,1)
        hr(nh3+1) = ul(idof,2)

c     Initialize element strains for deactivation

      elseif(isw.eq.18) then

        hr(nh3  ) = 0.0d0
        hr(nh3+1) = 0.0d0

c     External node check

      elseif(isw.eq.26) then

      endif

c     I/O Formats

2000  format(9x,'C o n t a c t   G a p   E l e m e n t'//
     &      10x,'Direction         =',i5/
     &      10x,'Degree of freedom =',i5/
     &      10x,'Penalty value     =',1p,e12.5/)

2001  format(5x,'Contact element',i5,', Force =',1p,1e13.5,
     &       ' Gap =',1p,1e13.5)

3000  format(' Input: DIREction, value '/'        DEGRee of freedom'/
     &       '        PENAlty,   value '/)

4000  format(' *ERROR* Element',i7,' has nodes',2i8)

4001  format(' *ERROR* Element',i7,' has node-3:',i8,
     &       'less than nodes 1,2:',2i8)

4002  format(' *ERROR* Material',i3,' No direction specified')

      end
