c$Id:$
      subroutine aropk(n,u,v,w,mode)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Author: D.S. Bindel 6/21/2005
c       Original version                                    01/11/2006
c       1. Add computation of A*v for mode = -1             03/09/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Apply K^-1 for Arnoldi or Lanczos

c      Inputs:
c         n    - Number of equations active
c         u(*) - Current vector
c         mode - Solution mode value

c      Temporary array
c         w(*) - Temporary vector

c      Outputs:
c         mode:1
c           v(*) - Solution: L*K*L*u = v; L = sqrt(M)-inverse
c         mode:2
c           v(*) - Solution: inv[M]*K*u = v
c         mode:3
c           v(*) - Solution: K * v = u
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'compas.h'
      include   'evdata.h'
      include   'fdata.h'
      include   'ndata.h'
      include   'part0.h'
      include   'p_int.h'
      include   'pointer.h'
      include   'comblk.h'

      integer    mode, i, n
      real*8     u(n),v(n),w(n)

c     Mode 1:

      if(mode.eq.1) then
        call vecpmult(hr(np(157)),u,v,n)
        if(ittyp.eq.-3) then
          do i = 1,n
            w(i) = 0.0d0
          end do ! i
          call promul(hr(nal),hr(nau),hr(na),v,w, mr(np(20+npart)),
     &                n,.true.)
        elseif(ittyp.ge.1) then
          call caprod(hr(na),hr(nau),v,w, mr(np(93)),mr(np(94)), n)
        endif
        call vecpmult(hr(np(157)),w,v,n)
      elseif(mode.eq.2) then

c     Mode 3:

      elseif(mode.eq.3) then
        fp(1)  = na
        fp(2)  = nau
        fp(3)  = nal
        fp(4)  = np(20+npart)
        do i = 1,n
          v(i) = u(i)
        end do ! i
        call psolve(ittyp,v,fp,.false.,.true.,.true.,.false.)

      elseif(mode.eq.-1) then
        fp(1) = na
        fp(2) = nau
        fp(3) = nal
        fp(4) = np(20+npart)
        call pamult(ittyp,u,v,fp)
      endif

      end

      subroutine vecpmult(m,u,v,neq)

      implicit   none

      integer    neq,n
      real*8     m(*),u(*),v(*)

      do n = 1,neq
        v(n) = m(n)*u(n)
      end do ! n

      end
