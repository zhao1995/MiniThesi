c$Id:$
      subroutine ctied2d (ndm,ndf,x,u,csw,npair,cs02,cp0,
     &                    ix1,ix2,ch1,ch2,ch3,ww1,ww3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Introduce 'setix' to ensure all surfaces set     15/01/2008
c       2. Remove 'slavfl' and use xi1-xi2 to check facet   18/07/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor             6 March 2003            1.0

c      Acronym: Contact DRIVER for 2D Tied

c      Purpose: Management of specific contact formulation

c      Inputs :
c         ndm     - Space dimension of mesh
c         ndf     - Number dof/node
c         x(*)    - Nodal coordinates
c         u(*)    - Current nodal solution vectors
c         csw     - Contact switch
c         npair   - # of current pair
c         cs02(*) - Contact surface 2 control data
c         cp0(*)  - Contactpair control data
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2

c      Outputs:
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)
c         w1(*)   - Dictionary of variables for CH1 & CH2
c         w3(*)   - Dictionary of variables for CH3
c                 - Data exchange with main program subroutine calls
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include   'cdata.h'
      include   'c_0.h'
      include   'c_comnd.h'
      include   'c_contac.h'
      include   'c_geom.h'
      include   'c_keyh.h'
      include   'c_mate.h'
      include   'c_pair.h'
      include   'c_tole.h'

      include   'compac.h'
      include   'compas.h'
      include   'counts.h'
      include   'debugs.h'
      include   'eqsym.h'
      include   'iofile.h'
      include   'ndata.h'
      include   'print.h'
      include   'umac1.h'

      include   'pointer.h'
      include   'comblk.h'

      character ww1(*)*(*),ww3(*)*(*)
      integer   ndm,ndf
      real*8    x(ndm,*),u(ndf,*)
      real*8    cs02(nr0,n0c1:*), cp0(nr0,n0c3:*)
      real*8    ch1(lh1,*),ch2(lh1,*),ch3(lh3,*)

      integer   csw,npair,ix1(dnope1,*),ix2(dnope2,*), ida(6)
      integer   ixl(30),is(10),js(3),isz(2)
      integer   fel,lel,nel,mel,nlm,lint, ke,kn, mm,me,mn,i,j,k,l
      integer   me1,me2, setix
      real*8    sg(2,10),shp(2,4),shpm(2,4),shph(2,4), zlen
      real*8    xi1,xi2, sgc,xjac,denom, xi0,gn0, zmax,tolz
      real*8    xs(2,4),zs(2,4),xm(2,4),zm(2,4),lm(2,4)
      real*8    cxm(2),dxm(2),clm(2),cxs(2),nn(2)
      real*8    tanm(36,36), resn(36), aa(4,4),bb(4)

      logical   ifprt,once,oneck, masfl,cmast2,cmasp2
      logical   slavt2

      save

c     Set active dof and dof order (idl length = ndf)

      data      ida     /1,2,3,4,5,6/
      data      tolz    / 1.d-8 /

      call cdebug0 ('    ctied2d',csw)

c-----[--.----+----.----+----.-----------------------------------------]
c     Automatic call: Set once

      if (csw.eq.0) then

        once  = .true.

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PCONTR for activation of requested history variables

      elseif (csw.eq.1) then

        if (once) then
          once  = .false.

c         Number of quadrature points per segment

          lint  = nint(cp0(4,0))
          if(lint.eq.0) then
            lint = 4
          endif
          if(lint.gt.6) then
            write(iow,3000) lint
            write(ilg,3000) lint
            call plstop()
          endif
          write(iow,*) ' LINT =',lint,cp0(4,0)
          write(  *,*) ' LINT =',lint,cp0(4,0)
          call int1d(lint,sg) ! Gauss-Legendre

c         Define history variables

          call defhvt2(ww1,ww3)
          do i = 1,4
            is(i) = (i-1)*ndf
          end do ! i
          js(1) = 0

        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from FORMFE to compute stiffness and residual

      elseif ((csw.eq.3).or.(csw.eq.6).or.(csw.eq.206)) then

        do ke = 1,neps1

c         Loop over slave surfaces

          nel = 0
          nlm = 0
          zmax = 0.0d0
          do kn = 1,nope1
            if(ix1(kn,ke).gt.0) then
              nel       = nel + 1
              ixl(nel)  = ix1(kn,ke)
              zs(1,nel) = x(1,ixl(nel))
              zs(2,nel) = x(2,ixl(nel))
              xs(1,nel) = x(1,ixl(nel)) + u(1,ixl(nel))
              xs(2,nel) = x(2,ixl(nel)) + u(2,ixl(nel))
              lm(1,nel) = ch2(p1(1)+nlm  ,ke)
              lm(2,nel) = ch2(p1(1)+nlm+1,ke)
              zmax      = max(zmax,abs(zs(1,nel)),abs(zs(2,nel)))
              nlm       = nlm + ndf
            endif
          end do ! kn
          zmax = zmax*tolz

c         Loop over Master surfaces

          me1 = nint(ch3(p3(1)  ,ke))
          me2 = nint(ch3(p3(1)+1,ke))

          if(me1.gt.0) then

c           Compute dual shape functions at Gauss points

            do i = 1,nel
              bb(i) = 0.0d0
              do j = 1,nel
                aa(j,i) = 0.0d0
              end do ! j
            end do ! i
            zlen = 0.0d0

            do l = 1,lint
              call shp1dn(sg(1,l),shp,nel)
              do j = 1,2
                dxm(j) = 0.0d0
                do i = 1,nel
                  dxm(j) = dxm(j) + shp(1,i)*zs(j,i)
                end do ! i
              end do ! j
              xjac =  sqrt(dxm(1)**2 + dxm(2)**2)*sg(2,l)
              zlen = zlen + xjac
              do i = 1,nel
                sgc   = shp(2,i)*xjac
                bb(i) = bb(i) + sgc
                do j = 1,nel
                  aa(j,i) = aa(j,i) + shp(2,j)*sgc
                end do ! j
              end do ! i
            end do ! l

            if(zlen.lt.zmax) then
              go to 100
            endif

c           Invert and multiply by integral of shape function

            call invert(aa,nel,4)
            do i = 1,nel
              do j = 1,nel
                aa(j,i) = aa(j,i)*bb(i)
              end do ! j
            end do ! i

c           Compute tangent and residual

            me = me1
            do mm = me1,me2,-1

c             Zero tangent and residual

              do j = 1,36
                resn(j) = 0.0d0
                do i = 1,36
                  tanm(i,j) = 0.0d0
                end do ! i
              end do ! j

c             Set master reference point

              do i = 1,2
                cxm(i) = x(i,ix2(1,me))
              end do ! i
              masfl = slavt2( cxm,zs,nel,xi1 )
              do i = 1,2
                cxm(i) = x(i,ix2(2,me))
              end do ! i
              masfl = slavt2( cxm,zs,nel,xi2 )
              xi1 = min( 1.d0,xi1)
              xi2 = max(-1.d0,xi2)

c             If projects on positive segment of slave facet

              if(abs(xi1-xi2).gt.0.0d0) then
                js(2) = nel*ndf
                do l = 1,lint
                  sgc   = 0.5d0*((xi1 - xi2)*sg(1,l) + xi1 + xi2)
                  masfl = cmasp2(sgc,x,zs,xm,ix2,ndm,nel,
     &                           me,xi0,gn0,xjac,nn)
                  if(masfl) then

c                   Do quadrature for element facet

                    call shp1dn(sgc,shp,nel)

c                   Form Dual shape functions for multiplier

                    do i = 1,nel
                      shph(1,i) = 0.0d0
                      shph(2,i) = 0.0d0
                      do j = 1,nel
                        shph(1,i) = shph(1,i) + aa(j,i)*shp(1,j)
                        shph(2,i) = shph(2,i) + aa(j,i)*shp(2,j)
                      end do ! j
                    end do ! i

c                   Form slave and lagrange multiplier parts

                    cxs(1) = 0.0d0
                    cxs(2) = 0.0d0
                    clm(1) = 0.0d0
                    clm(2) = 0.0d0
                    do i = 1,nel
                      cxs(1) = cxs(1) + shp (2,i)*xs(1,i)
                      cxs(2) = cxs(2) + shp (2,i)*xs(2,i)
                      clm(1) = clm(1) + shph(2,i)*lm(1,i)
                      clm(2) = clm(2) + shph(2,i)*lm(2,i)
                    end do ! i

c                   Set nodes on master surface

                    mel   = 0
                    do mn = 1,nope2
                      if(ix2(mn,me).gt.0) then
                        mel          = mel + 1
                        ixl(mel+nel) = ix2(mn,me)
                        zm(1,mel)    = x(1,ixl(mel+nel))
                        zm(2,mel)    = x(2,ixl(mel+nel))
                        xm(1,mel)    = x(1,ixl(mel+nel))
     &                               + u(1,ixl(mel+nel))
                        xm(2,mel)    = x(2,ixl(mel+nel))
     &                               + u(2,ixl(mel+nel))
                      endif
                    end do ! mn

c                   Compute the master point

                    cxm(1) = gn0*nn(1)
                    cxm(2) = gn0*nn(2)
                    xjac   = xjac*sg(2,l)*(xi1-xi2)*0.5d0
                    call shp1dn(xi0,shpm,mel)
                    do i = 1,mel
                      cxm(1) = cxm(1) + shpm(2,i)*xm(1,i)
                      cxm(2) = cxm(2) + shpm(2,i)*xm(2,i)
                    end do ! i

c                   Compute residual

                    js(3) = js(2) + mel*ndf
                    do j = 1,nel
                      do i = 1,2
                        resn(js(1)+is(j)+i) = resn(js(1)+is(j)+i)
     &                                      - shp(2,j)*clm(i)*xjac
                        resn(js(3)+is(j)+i) = resn(js(3)+is(j)+i)
     &                                + shph(2,j)*(cxm(i)-cxs(i))*xjac
                      end do ! i
                    end do ! j
                    do j = 1,mel
                      do i = 1,2
                        resn(js(2)+is(j)+i) = resn(js(2)+is(j)+i)
     &                                      + shpm(2,j)*clm(i)*xjac
                      end do ! i
                    end do ! j

c                   Compute tangent

                    if(csw.eq.3) then
                      do j = 1,nel
                        do k = 1,nel
                          denom = shp(2,k)*shph(2,j)*xjac
                          do i = 1,2
                            tanm(js(1)+is(k)+i,js(3)+is(j)+i) =
     &                      tanm(js(1)+is(k)+i,js(3)+is(j)+i) + denom
                            tanm(js(3)+is(j)+i,js(1)+is(k)+i) =
     &                      tanm(js(3)+is(j)+i,js(1)+is(k)+i) + denom
                          end do ! i
                        end do ! k
                        do k = 1,mel
                          denom = shpm(2,k)*shph(2,j)*xjac
                          do i = 1,2
                            tanm(js(2)+is(k)+i,js(3)+is(j)+i) =
     &                      tanm(js(2)+is(k)+i,js(3)+is(j)+i) - denom
                            tanm(js(3)+is(j)+i,js(2)+is(k)+i) =
     &                      tanm(js(3)+is(j)+i,js(2)+is(k)+i) - denom
                          end do ! i
                        end do ! k
                      end do ! j
                    endif
                  endif
                end do ! l

c               Modify residual/tangent for fixed b.c.

                isz(1) = 0
                isz(2) = 0
                if(    (abs(zm(1,1) - zs(1,2)).lt.zmax  .and.
     &                  abs(zm(2,1) - zs(2,2)).lt.zmax) ) then
                  isz(2) = ixl(nel+1)
                endif

                if(    (abs(zm(1,2) - zs(1,1)).lt.zmax  .and.
     &                  abs(zm(2,2) - zs(2,1)).lt.zmax) ) then
                  isz(1) = ixl(nel+2)
                endif

                call cmodify(isz,ixl,mr(np(31)+ndf*numnp),
     &                       nel,mel,ndf, 36, tanm,resn)

c               Save stiffness/residual

                i = nel + mel
                call constass(ixl,ida,i,ndf,ixl,nel,ndf,36,tanm,resn)

              endif ! abs(xi1-xi2) .gt. 0.0d0
              me = ix2(dnope2-1,me) ! set next adjacent element
            end do ! mm
100         continue
          endif ! me1 > 0


        end do ! ke

c-----[--.----+----.----+----.-----------------------------------------]
c     History variables initialization

      elseif (csw.eq.14) then

c-----[--.----+----.----+----.-----------------------------------------]
c     CSW = 103: Called from PMACR1 to check geometry of contacts
c     CSW = 304: Called from PMACR3 to check geometry of contacts

      elseif (csw.eq.103 .or. csw.eq.304) then

c       Only do once

        setix = max(setix,npair)
        if(oneck) then
          if(npair.lt.setix) then
            oneck = .false.
          endif
          do ke = 1,neps1
            nel = 0
            do kn = 1,nope1
              if(ix1(kn,ke).gt.0) then
                nel       = nel + 1
                ixl(nel)  = ix1(kn,ke)
                xs(1,nel) = x(1,ixl(nel))
                xs(2,nel) = x(2,ixl(nel))
              endif
            end do ! kn

c           Check how many master elements intersect slave surface

            masfl = cmast2(-1.d0,x,xs,xm,ix2,ndm,nel,
     &                     me,xi0,gn0,xjac,nn)
            if(masfl .and. abs(xi0).le.1.d0+tolz) then
              me1 =  me
              xi1 =  xi0
            else
              me1 = -me
              xi1 = -1.d0
            endif
            masfl = cmast2( 1.d0,x,xs,xm,ix2,ndm,nel,
     &                     me,xi0,gn0,xjac,nn)
            if(masfl .and. abs(xi0).le.1.d0+tolz) then
              me2 =  me
              xi2 =  xi0
            else
              me2 = -me
              xi2 =  1.d0
            endif
            if(me1.lt.0 .and. me2.lt.0) then
            elseif(me1.lt.0) then
              me1 = me2
            elseif(me2.lt.0) then
              me2 = me1
            endif

            ch3(p3(1)  ,ke) = me1
            ch3(p3(1)+1,ke) = me2
          end do ! ke
        endif ! oneck

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR5 to show element information

      elseif (csw.eq.200) then

         write (*,2000)

c-----[--.----+----.----+----.-----------------------------------------]
c     Print contact status

      elseif(csw.eq.204) then

c       Get print flag and range

        call setcprt (ifprt,fel,lel)

c       Print title

        if (ifprt) then
          write (iow,2001) npair
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PPLOTF for plot of contact geometry

      elseif (csw.eq.305) then
        if(nint(cs02(1,0)).eq.6) then
          call c2rigplt (ix1,cs02(2,0),2,8)
        else
          call c2geoplt (ix1,ix2,2,8)
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PPLOTF to set profile and range to plot variable
c     Called from CONTACT for plot contours of a contact variable

      elseif (csw.eq.308 .or. csw.eq.408) then
         call c2varplt (ix1,ch1,ch2,ch3,npair,csw)

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PCONTR for initialization of history variables

      elseif (csw.eq.313) then

        nset = neps1
        call acthvt2 (nset)

c-----[--.----+----.----+----.-----------------------------------------]
c     Updates for Lagrange multipliers

      elseif (csw.eq.314) then

        if(ifsolm.eq.2) then
          do ke = 1,neps1
            nel = 0
            do kn = 1,nope1
              if(ix1(kn,ke).gt.0) then
                nel       = nel + 1
                ixl(nel)  = ix1(kn,ke)
              endif
            end do ! kn
            if(nel.gt.0) then
              call getlagm(ixl,nel,ndf,ch2(p1(1),ke))
            endif
          end do ! ke
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PCONTR resets flag for "once action" for next problem

      elseif (csw.eq.400) then
        oneck = .true.
        setix = 0

c-----[--.----+----.----+----.-----------------------------------------]
c     Reset profile for contact

      elseif (csw.eq.403) then

c       Check each slave element

        do ke = 1,neps1
          nel = 0
          do kn = 1,nope1
            if(ix1(kn,ke).gt.0) then
              nel       = nel + 1
              ixl(nel)  = ix1(kn,ke)
            endif
          end do ! kn

c         Do master elements

          me1   = nint(ch3(p3(1)  ,ke))
          if(me1.gt.0) then
            me2 = nint(ch3(p3(1)+1,ke))
            do me = me1,me2,-1

c             Set master segment

              mel   = nel
              do mn = 1,nope2
                if(ix2(mn,me).gt.0) then
                  mel       = mel + 1
                  ixl(mel)  = ix2(mn,me)
                endif
              end do ! mn
              call modprofl(ixl,ida,mel,ndf,ixl,nel,ndf)
            end do ! me
          else
            write(*,*) ' ME1 = 0: No Master'
          endif
        end do ! ke

      endif

c     Formats

2000  format(//10x,'2-D Tied Interface - Dual form')
2001  format (/'     I n t e r f a c e   O u t p u t   f o r   P a i r',
     &         i5)

3000  format(/' **ERROR** Too many quadrature points for segment'/
     &        '           Number requested =',i4/
     &        '           Maximum number   =  10'/)

      end
