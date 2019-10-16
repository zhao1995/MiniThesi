c$Id:$
      subroutine elmlib(d,u,x,ix,t,s,p,i,j,k,jel,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Element library driver routine
c               N.B. Must set library flags in Subroutine PELNUM
c                    for new program modules

c      Inputs:
c         d(*)    - Material parameters
c         u(*)    - Element solution parameters
c         x(*)    - Element nodal coordinates
c         ix(*)   - Element nodal numbers
c         t(*)    - Element temperatures
c         i       - Number dof/node           (ndf)
c         j       - Spatial dimension of mesh (ndm)
c         k       - Size of element arrays    (nst)
c         jel     - Element type number

c      Outputs:
c         d(*)    - Material parameters (isw = 1 only)
c         s(*,*)  - Element array
c         p(*)    - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'iofile.h'

      include  'setups.h'

      integer   i,j,k,jel,isw, ix(*)
      real*8    p(*),s(*),d(*),u(*),x(*),t(*)

      save

      if(isw.ge.3 .and. k.gt.0) then
        call pzero(s,k*k)
        call pzero(p,k  )
      endif

c     User element routines

      if(jel.gt.0) then

        if(    jel.eq. 1) then
          call elmt01(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq. 2) then
          call elmt02(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq. 3) then
          call elmt03(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq. 4) then
          call elmt04(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq. 5) then
          call elmt05(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq. 6) then
          call elmt06(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq. 7) then
          call elmt07(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq. 8) then
          call elmt08(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq. 9) then
          call elmt09(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.10) then
          call elmt10(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.11) then
          call elmt11(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.12) then
          call elmt12(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.13) then
          call elmt13(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.14) then
          call elmt14(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.15) then
          call elmt15(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.16) then
          call elmt16(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.17) then
          call elmt17(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.18) then
          call elmt18(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.19) then
          call elmt19(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.20) then
          call elmt20(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.21) then
          call elmt21(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.22) then
          call elmt22(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.23) then
          call elmt23(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.24) then
          call elmt24(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.25) then
          call elmt25(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.26) then
          call elmt26(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.27) then
          call elmt27(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.28) then
          call elmt28(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.29) then
          call elmt29(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.30) then
          call elmt30(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.31) then
          call elmt31(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.32) then
          call elmt32(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.33) then
          call elmt33(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.34) then
          call elmt34(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.35) then
          call elmt35(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.36) then
          call elmt36(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.37) then
          call elmt37(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.38) then
          call elmt38(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.39) then
          call elmt39(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.40) then
          call elmt40(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.41) then
          call elmt41(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.42) then
          call elmt42(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.43) then
          call elmt43(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.44) then
          call elmt44(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.45) then
          call elmt45(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.46) then
          call elmt46(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.47) then
          call elmt47(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.48) then
          call elmt48(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.49) then
          call elmt49(d,u,x,ix,t,s,p,i,j,k,isw)
        elseif(jel.eq.50) then
          call elmt50(d,u,x,ix,t,s,p,i,j,k,isw)
        else
          go to 400
        endif

c     Program element library routines

      else

c       1-D Element library

        if(j.eq.1) then

          if(    jel.eq. -1) then
            call solid1d(d,u,x,ix,t,s,p,i,j,k,isw)
          elseif(jel.eq. -2) then
            call trussnd(d,u,x,ix,t,s,p,i,j,k,isw)
          elseif(jel.eq. -3) then
            write(iow,*) ' No 1-d frame element available. n =',n
            call plstop()
          elseif(jel.eq. -4) then
            write(iow,*) ' No 1-d plate element available. n =',n
            call plstop()
          elseif(jel.eq. -5) then
            write(iow,*) ' No 1-d shell element available. n =',n
            call plstop()
          elseif(jel.eq. -6) then
            write(iow,*) ' No 1-d membrane element available. n =',n
            call plstop()
          elseif(jel.eq. -7) then
            call gapnd(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -8) then
            call therm1d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -9) then
            call convec1d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq.-10) then
            call point(d,u,ix,s,p,i,j,k,isw)
          elseif(jel.eq.-11) then
            write(iow,*) ' No 1-d pressure load available. n =',n
            call plstop()
          elseif(jel.eq.-12) then
            write(iow,*) ' No 1-d St.Venant torsion available. n =',n
            call plstop()
          elseif(jel.eq.-13) then
            write(iow,*) ' No 1-d coupled thermo-mechanical',
     &                   '  available. n =',n
            call plstop()
          else
            go to 400
          endif

c       2-D element library

        elseif(j.eq.2) then

          if(    jel.eq. -1) then
            call solid2d(d,u,x,ix,t,s,p,i,j,k,isw)
          elseif(jel.eq. -2) then
            call trussnd(d,u,x,ix,t,s,p,i,j,k,isw)
          elseif(jel.eq. -3) then
            call frame2d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -4) then
            call plate2d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -5) then
            call shell2d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -6) then
            write(iow,*) ' No 2-d membrane element available. n =',n
            call plstop()
          elseif(jel.eq. -7) then
            call gapnd(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -8) then
            call therm2d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -9) then
            call convec2d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq.-10) then
            call point(d,u,ix,s,p,i,j,k,isw)
          elseif(jel.eq.-11) then
            call presld(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq.-12) then
            call stvtor(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq.-13) then
            write(iow,*) ' No 2-d coupled thermo-mechanical',
     &                   '  available. n =',n
            call plstop()
          else
            go to 400
          endif

c       3-D Element library

        elseif(j.eq.3) then

          if(    jel.eq. -1) then
            call solid3d(d,u,x,ix,t,s,p,i,j,k,isw)
          elseif(jel.eq. -2) then
            call trussnd(d,u,x,ix,t,s,p,i,j,k,isw)
          elseif(jel.eq. -3) then
            call frame3d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -4) then
            write(iow,*) ' No 3-d plate element available: Use SHELL.',n
            call plstop()
          elseif(jel.eq. -5) then
            call shell3d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -6) then
            call membr3d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -7) then
            call gapnd(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -8) then
            call therm3d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq. -9) then
            call convec3d(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq.-10) then
            call point(d,u,ix,s,p,i,j,k,isw)
          elseif(jel.eq.-11) then
            call presld(d,u,x,ix,s,p,i,j,k,isw)
          elseif(jel.eq.-12) then
            write(iow,*) ' No 3-d St.Venant torsion available. n =',n
            call plstop()
          elseif(jel.eq.-13) then
            call solid3dtm(d,u,x,ix,t,s,p,i,j,k,isw)
          else
            go to 400
          endif

        endif
      endif

      return

c     Error

400   write(ilg,4000) n,jel,isw, rank
      if(ior.gt.0) write(iow,4000) n,jel,isw, rank
      if(ior.lt.0) write(  *,4000) n,jel,isw, rank
      call plstop()

c     Format

4000  format('  *ERROR* ELMLIB: Element:',i6,', type number',i3,
     &       ' input, isw =', i3/' RANK =',i4)

      end
