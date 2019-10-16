c$Id:$
      subroutine ielmlib(d1,u1,x1,ix1,t1,d2,u2,x2,ix2,t2,
     &                   intnod,s,p,isw)

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
c         di(*)  - Material parameters
c         ui(*)  - Element solution parameters
c         xi(*)  - Element nodal coordinates
c         ixi(*) - Element nodal numbers
c         ti(*)  - Element temperatures
c         isw    - Element type number

c      Outputs:
c         d1(*)  - Material parameters (isw = 1 only)
c         s(*,*) - Element array
c         p(*)   - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'ieldat.h'

      integer    isw, i,j
      integer    ix1(*),ix2(*), intnod(*)
      real*8     d1(*),u1(*),x1(*),t1(*),d2(*),u2(*),x2(*),t2(*)
      real*8     s(nsts,nsts),p(nsts)

      save

      if(isw.ge.3 .and. nsts.gt.0) then
        do i = 1,nsts
          p(i) = 0.0d0
          do j = 1,nsts
            s(j,i) = 0.0d0
          end do ! j
        end do ! i
      endif

      if(iel1.eq. 1) then
        call ielmt01(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq. 2) then
        call ielmt02(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq. 3) then
        call ielmt03(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq. 4) then
        call ielmt04(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq. 5) then
        call ielmt05(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq. 6) then
        call ielmt06(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq. 7) then
        call ielmt07(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq. 8) then
        call ielmt08(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq. 9) then
        call ielmt09(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.10) then
        call ielmt10(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.11) then
        call ielmt11(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.12) then
        call ielmt12(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.13) then
        call ielmt13(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.14) then
        call ielmt14(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.15) then
        call ielmt15(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.16) then
        call ielmt16(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.17) then
        call ielmt17(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.18) then
        call ielmt18(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.19) then
        call ielmt19(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.20) then
        call ielmt20(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.21) then
        call ielmt21(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.22) then
        call ielmt22(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.23) then
        call ielmt23(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.24) then
        call ielmt24(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.25) then
        call ielmt25(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.26) then
        call ielmt26(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.27) then
        call ielmt27(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.28) then
        call ielmt28(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.29) then
        call ielmt29(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.30) then
        call ielmt30(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.31) then
        call ielmt31(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.32) then
        call ielmt32(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.33) then
        call ielmt33(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.34) then
        call ielmt34(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.35) then
        call ielmt35(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.36) then
        call ielmt36(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.37) then
        call ielmt37(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.38) then
        call ielmt38(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.39) then
        call ielmt39(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.40) then
        call ielmt40(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.41) then
        call ielmt41(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.42) then
        call ielmt42(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.43) then
        call ielmt43(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.44) then
        call ielmt44(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.45) then
        call ielmt45(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.46) then
        call ielmt46(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.47) then
        call ielmt47(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.48) then
        call ielmt48(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.49) then
        call ielmt49(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      elseif(iel1.eq.50) then
        call ielmt50(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,p,isw)
      else
        write(*,*) ' INTERFACE NOT AVAILABLE'
      endif

      end
