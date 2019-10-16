c$Id:$
      subroutine fld3d2(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'pconstant.h'                 14/11/2006
c       2. Revise statement for augmenting                  26/01/2007
c       3. Modify augmenting for higher order elements      13/03/2007
c          Use d(185) for augmenting factor
c       4. Correct computation of 4-node tet intertia;      29/03/2007
c          Add definition of 'cfac','lfac'.
c       5. Add 14/15 node tet computation                   06/09/2007
c          Set number of element multipliers for 14/15 tet  07/09/2007
c       6. Introduce 'nvn' to compute center coordinate     14/09/2007
c       7. Move 'xu' to isw.eq.4 area, redefine press(l),   25/09/2007
c          Eliminate need for 'mpress' and 'dsigtr'
c       8. Change 'ord' to 'nel' on 'tetshp' calls          05/11/2007
c       9. Modify to compute inerta and internal force      21/12/2007
c          separately.
c      10. Revise quadrature order for tetrahedra           28/12/2007
c          Modify form for computing intertia
c      11. Correct prints of logarithmic lambdas            28/08/2008
c      12. Add direct call to quadr3d and interp3d          11/11/2008
c      13. Remove 'nel' from call to 'quadr3d'              23/01/2009
c      14. Increase arrays to store 64 node brick           03/02/2009
c          Add quadratic terms to phi array.
c      15. Revise augment: augfp -> augf; nhi - > fp(1)     14/04/2009
c      16. Add 'defgrd.h' for deformation gradient values   27/04/2009
c      17. Dimension theta(3,125)                           05/06/2009
c      18. Set 'npm' for augment by call to quadr3d         15/08/2009
c      19. Add prints for Almansi strains                   19/10/2009
c      20. Correct phi( 6,*) to phi( 8,*); also phi(10,*)   21/11/2010
c      21. Increase arrays to store 125 node brick          20/12/2010
c      23. Increase dimension on hh,rr,hsig to 10           22/12/2010
c      24. Move inertial computation after stiffness        01/01/2011
c          Add isw to iner3d to control tangent computation
c      25. Add computation of 'dvol' for stress computation 29/03/2011
c      26. Convert to computations using the displacement   21/11/2011
c          gradient. dimension theta(4,*), etc.
c      27. Add average of density for multiscale            10/05/2012
c      28. Add eps to slcn3d call                           01/01/2013
c      29. Initialize sigm(1) for outputs                   05/03/2013
c      30. Use pengy3d to compute momenta and energies      06/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Finite deformation mixed model: u-p-theta

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         r(ndf,*)  - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'defgrd.h'
      include  'elcoor.h'
      include  'eldata.h'
      include  'elengy.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'oelmt.h'
      include  'pconstant.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'prstrs.h'
      include  'ptdat6.h'
      include  'p_int.h'
      include  'qudshp.h'
      include  'rdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   bflg
      integer   ndf,ndm,nst,isw,i,ii,i1,j,jj,j1,l,nhv,nn, istrt

      real*8    xlamd,ha, d1, epp, thlog, dtheta,qfact, ta
      real*8    dpress,bdb
      real*8    d(*),  ul(ndf,nen,*),  xl(ndm,*), s(nst,*),r(ndf,*)

      real*8    bbd(3,7),  bei(6), ad(7,7,5,125), dd(7,7), sigm(10)
      real*8    dvol(125), xxm(3), xr(3,125),ur(3,125), epsl(6,125)
      real*8    sigl(10,125),shpbar(3,125,125),weng(125),theta(4,125)
      real*8    xu(3,125),ru(3,125),bpra(3),body(3), bf(3),bt(3,3)
      real*8    hh(10,10),rr(10),hsig(10),phi(10,125),press(125),x0(3)

      save

c     TEMPORARY TEMPERATURE

      data      ta    / 0.0d0 /

c     Adjust data storage for augmenting

      if(isw.eq.1) then

        if(nen.le.10) then
          nh1 = nh1 + 2
        elseif(nen.le.27) then
          nh1 = nh1 + 8
        else
          nh1 = nh1 + 20
        endif

c     Augmented Lagrangian update for nested iteration

      elseif(isw.eq.10) then

        call quadr3d(d,.true.)

        d1  = augf*d(185)
        fp(1) = nh2 - 1
        fp(2) = fp(1) + npm
        do i = 1,npm
          hr(fp(1)+i) = hr(fp(1)+i) + d1*hr(fp(2)+i)
        end do ! i

c     Compute tangent stiffness and residual force vector

      elseif(isw.eq.3 .or. isw.eq.4  .or. isw.eq.6  .or.
     &       isw.eq.8 .or. isw.eq.13 .or. isw.eq.14 .or.
     &       isw.eq.16) then

       estore = 0.0d0

c      Set quadrature and order

       call quadr3d(d,.true.)

c      Mean coordinate of vertex nodes

       if(npm.gt.1) then
         do i = 1,3
           x0(i) = 0.0d0
           do j = 1,nvn
             x0(i) = x0(i) + xl(i,j)
           end do ! j
           x0(i) = x0(i)/dble(nvn)
         end do ! i
       endif

c      Get shape functions and derivatives in geometry at time t_n+1

       do l = 1,lint

         call interp3d(l, xl, ndm,nel)

         do i = 1,3
           xr(i,l) = 0.0d0
           ur(i,l) = 0.0d0
           do j = 1,nel
             xr(i,l) = xr(i,l) + shp3(4,j,l)*xl(i,j)
             ur(i,l) = ur(i,l) + shp3(4,j,l)*ul(i,j,1)
           end do ! j
c          ur(i,l) = ur(i,l) + xr(i,l)
         end do ! j

         phi(1,l) = 1.d0
         if(npm.gt.1) then
           phi(2,l) = xr(1,l) - x0(1)
           phi(3,l) = xr(2,l) - x0(2)
           phi(4,l) = xr(3,l) - x0(3)
           if(npm.gt.4) then
             phi( 5,l) = phi(2,l)**2
             phi( 6,l) = phi(2,l)*phi(3,l)
             phi( 7,l) = phi(3,l)**2
             phi( 8,l) = phi(3,l)*phi(4,l)
             phi( 9,l) = phi(4,l)**2
             phi(10,l) = phi(4,l)*phi(2,l)
           endif
         endif
       end do ! l

c      Set number of history terms / quadradure point

       nhv   = nint(d(15))
       istrt = nint(d(84))

c      MECHANICAL ELEMENT

       if(isw.eq. 3 .or. isw.eq. 6 .or.
     &    isw.eq.13 .or. isw.eq.14) then

c        Compute f, finv, df and det(fei) at conf t-n+1

         if(isw.eq.14) then

           call pfinit(f,df,finv,detf, lint)

         else

           call kine3d(shp3,ul,f,finv,df,detf,ndf,nel,nen,lint)

         endif

c        Mixed model for volumetric response

         call bbar3m(phi,shp3,jac,detf,lint,nel,npm,hh,theta,shpbar)

c        Compute mixed model deformation gradient

         call fbarm(f,detf,theta,lint)

c        Initialize storage for augmented function

         if(npm.gt.1) then
           do i = 1,npm
             rr(i) = 0.0d0
           end do ! i
         endif

c        Compute Cauchy stresses and spatial tangent tensor at t-n+1

         nn = 2*npm
         do l = 1,lint

           do i = 1,3
             xref(i) = xr(i,l)
             xcur(i) = xr(i,l) + ur(i,l)
           end do ! i

c          Set augmented multipler

           xlamd = hr(nh2)
           do i = 1,npm-1
             xlamd = xlamd + phi(i+1,l)*hr(nh2+i)
           end do ! i

           call modlfd(l,d,f(1,1,l),finv(1,l),df(1,l),theta(1,l),ta,
     &                 hr(nn+nh1),hr(nn+nh2),nhv,istrt,ad(1,1,1,l),
     &                 sigl(1,l),bei,xlamd,ha,.true.,isw)
           weng(l) = estore
           if(npm.gt.1) then
             do i = 1,npm
               rr(i) = rr(i) + phi(i,l)*ha*jac(l)
             end do ! i
           endif
           nn = nn + nhv
         end do ! l

         if(isw.eq.14) return

c        Set augmented function

         if(npm.eq.1) then
           hr(nh2+1) = ha
         else
           fp(1) = nh2 + npm  - 1
           do i = 1,npm
             hr(fp(1)+i) = 0.0d0
             do j = 1,npm
               hr(fp(1)+i) = hr(fp(1)+i) + hh(i,j)*rr(j)
             end do ! j
           end do ! i
         endif

c        ENERGY COMPUTATIONS

         if(isw.eq.13) then

           call pengy3d(ndm,ndf,nel, weng, d,ul,xl)

c        Compute mixed pressure

         elseif(isw.eq.3 .or. isw.eq.6) then

c          Set body forc values

           call sbodyf(d, body)

           if(npm.eq.1) then

             press(1) = 0.0d0
             do l = 1,lint

c             Modify volume element and integrate pressure
c             over reference volume

              press(1) = press(1) + one3*(sigl(1,l) + sigl(2,l)
     &                                  + sigl(3,l))*jac(l)
              dvol(l)  = jac(l) * theta(1,l)

             end do ! l

c            Divide pressure by reference volume

             press(1) = press(1) * hh(1,1)
             do l = 2,lint
               press(l) = press(1)
             end do ! l
           else
             do i = 1,npm
               sigm(i) = 0.0d0
             end do ! i

             do l = 1,lint

c              Modify volume element and integrate pressure
c              over reference volume

               dpress  = one3*(sigl(1,l)+sigl(2,l)+sigl(3,l))*jac(l)
               sigm(1) = sigm(1) + dpress
               do j = 2,npm
                 sigm(j) = sigm(j) + dpress*phi(j,l)
               end do ! j
               dvol(l) = jac(l) * theta(1,l)

             end do ! l

c            Divide pressure by reference volume

             do i = 1,npm
               hsig(i) = 0.0d0
               do j = 1,npm
                 hsig(i) = hsig(i) + hh(i,j)*sigm(j)
               end do ! j
             end do ! i
             do l = 1,lint
               press(l) = hsig(1)
               do j = 2,npm
                 press(l) = press(l) + hsig(j)*phi(j,l)
               end do ! j
             end do ! l

           endif ! npm
           bflg  = d(4).gt.0.0d0 .and. d(65).gt.0.0d0

c          Compute final residual and tangent arrays

           do l = 1,lint

c            Compute parameters for multiscale plane strain

             v_avg = v_avg + jac(l)
             v_rho = v_rho + jac(l)*d(4)

c            Angular velocity: d(4) = rho; d(65) = omega

             do i = 1,3
               bf(i) = 0.0d0
             end do ! i
             if(bflg) then
               call sbodyw(d(4),d(65),ur(1,l), bf,bt, .true.)
             endif

c            Compute mixed stress and multiply by volume element

             press(l) = press(l)*detf(1,l)/theta(1,l)
             dpress   = press(l) - (sigl(1,l)+sigl(2,l)+sigl(3,l))*one3

             sigm(1)  =  sigl(1,l) + dpress
             sigm(2)  =  sigl(2,l) + dpress
             sigm(3)  =  sigl(3,l) + dpress
             sigm(4)  =  sigl(4,l)
             sigm(5)  =  sigl(5,l)
             sigm(6)  =  sigl(6,l)

c            Store time history plot data for element

             i = 6*(l-1)
             do j = 1,6
               tt(j+i) = sigm(j)
               sigm(j) = sigm(j)*dvol(l)
             end do ! j

c            Compute residual

             do j = 1,nel

               ru(1,j) = shp3(1,j,l)*sigm(1)
     &                 + shp3(2,j,l)*sigm(4)
     &                 + shp3(3,j,l)*sigm(6)

               ru(2,j) = shp3(1,j,l)*sigm(4)
     &                 + shp3(2,j,l)*sigm(2)
     &                 + shp3(3,j,l)*sigm(5)

               ru(3,j) = shp3(1,j,l)*sigm(6)
     &                 + shp3(2,j,l)*sigm(5)
     &                 + shp3(3,j,l)*sigm(3)

               do i = 1,3
                 r(i,j)  = r(i,j) - ru(i,j)
     &                   + shp3(4,j,l)*jac(l)*(body(i) + bf(i))
               end do ! i
             end do ! j

c            Compute mixed tangent stiffness matrix

             if(isw.eq.3) then

c              Part 1: Geometric tangent matrix

               if(gflag) then

                 i1 = 0
                 do i = 1,nel
                   j1 = 0
                   do j = 1,nel
                     bdb = (shp3(1,i,l)*ru(1,j)
     &                   +  shp3(2,i,l)*ru(2,j)
     &                   +  shp3(3,i,l)*ru(3,j))*ctan(1)
                     do jj = 1,3
                       s(i1+jj,j1+jj) = s(i1+jj,j1+jj) + bdb
                     end do ! jj
                     j1 = j1 + ndf
                   end do ! j
                   i1 = i1 + ndf
                 end do ! i

               endif ! gflag

c              Part 2: Material tangent matrix

c              Modify tangent moduli for stress factors

               dpress = one3*(sigl(1,l) + sigl(2,l) + sigl(3,l))
               call dmatdx(ad(1,1,1,l),sigl(1,l),dpress,press(l))

c              Multiply tangent moduli by volume element

               d1 = dvol(l)*ctan(1)
               do j = 1,7
                 do i = 1,7
                   dd(i,j) = ad(i,j,1,l)*d1
                 end do ! i
               end do ! j

c              Compute row terms

               i1    = 0
               do i = 1,nel

c                Angular velocity tangent

                 if(bflg) then
                   do jj = 1,3
                     do ii = 1,3
                       bdb = shp3(4,i,l)*bt(ii,jj)
                       j1  = 0
                       do j = 1,nel
                         s(i1+ii,j1+jj) = s(i1+ii,j1+jj)
     &                                  + bdb*shp3(4,j,l)
                         j1             = j1 + ndf
                       end do ! j
                     end do ! ii
                   end do ! jj
                 endif ! bflg

c                Compute bmat-t * dd * dvol

                 do jj = 1,7

                   bbd(1,jj) =   shp3(1,i,l)*dd(1,jj)
     &                       +   shp3(2,i,l)*dd(4,jj)
     &                       +   shp3(3,i,l)*dd(6,jj)
     &                       + shpbar(1,i,l)*dd(7,jj)

                   bbd(2,jj) =   shp3(2,i,l)*dd(2,jj)
     &                       +   shp3(1,i,l)*dd(4,jj)
     &                       +   shp3(3,i,l)*dd(5,jj)
     &                       + shpbar(2,i,l)*dd(7,jj)

                   bbd(3,jj) =   shp3(3,i,l)*dd(3,jj)
     &                       +   shp3(2,i,l)*dd(5,jj)
     &                       +   shp3(1,i,l)*dd(6,jj)
     &                       + shpbar(3,i,l)*dd(7,jj)
                 end do ! jj

                 j1 = 0
                 do j = 1,nel

c                  Compute mechanics part of tangent stiffness

                   do jj = 1,3

                     s(i1+jj,j1+1) = s(i1+jj,j1+1)
     &                             + bbd(jj,1)*shp3(1,j,l)
     &                             + bbd(jj,4)*shp3(2,j,l)
     &                             + bbd(jj,6)*shp3(3,j,l)
     &                             + bbd(jj,7)*shpbar(1,j,l)

                     s(i1+jj,j1+2) = s(i1+jj,j1+2)
     &                             + bbd(jj,2)*shp3(2,j,l)
     &                             + bbd(jj,4)*shp3(1,j,l)
     &                             + bbd(jj,5)*shp3(3,j,l)
     &                             + bbd(jj,7)*shpbar(2,j,l)

                     s(i1+jj,j1+3) = s(i1+jj,j1+3)
     &                             + bbd(jj,3)*shp3(3,j,l)
     &                             + bbd(jj,5)*shp3(2,j,l)
     &                             + bbd(jj,6)*shp3(1,j,l)
     &                             + bbd(jj,7)*shpbar(3,j,l)
                   end do ! jj
                   j1 = j1 + ndf
                 end do ! j
                 i1 = i1 + ndf
               end do ! i
             endif ! isw = 3
           end do ! l

c          Compute inertia effects: shflg = .true. for eigen shifts

           if(ctan(3).ne.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
             call iner3d(d,xl,ul(1,1,4),ul(1,1,5), s,r,
     &                   nel,ndf,ndm,nst, isw)
           endif ! ctan(3) test

         endif ! isw = 3 or 6

c      Output stresses and computre fracture force

       elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16) then

c        Compute current geometry

         do i = 1,3
          do j = 1,nel
           xu(i,j) = xl(i,j) + ul(i,j,1)
          end do ! j
         end do ! i

c        Zero arrays

         do i = 1,10
           sigm(i) = 0.0d0
         end do ! i
         do i = 1,6
           ebig(i) = 0.0d0
           esml(i) = 0.0d0
         end do ! i
         do i = 1,3
           bpra(i) = 0.0d0
           xxm(i)  = 0.0d0
         end do ! i
         sigm(10) = 0.0d0
         epp      = 0.0d0
         dtheta   = 0.0d0
         qfact    = 1.d0/dble(lint)

c        Compute f, finv, df and det(fei) at conf t-n+1

         call kine3d(shp3,ul,f,finv,df,detf,ndf,nel,nen,lint)

c        Mixed model for volumetric response

         call bbar3m(phi,shp3,jac,detf,lint,nel,npm,hh,theta,shpbar)

c        Compute mixed model deformation gradient

         call fbarm(f,detf,theta,lint)

c        Second loop over Gauss points

         nn  = 2*npm
         do l = 1,lint

c          Compute Cauchy stresses and spatial tangent tensor at t-n+1

           do i = 1,3
             xref(i) = xr(i,l)
             xcur(i) = xr(i,l) + ur(i,l)
           end do ! i

c          Set augmented function

           xlamd = hr(nh2)
           do i = 1,npm-1
             xlamd = xlamd + phi(i+1,l)*hr(nh2+i)
           end do ! i

           call modlfd(l,d,f(1,1,l),finv(1,l),df(1,l),theta(1,l),ta,
     &                 hr(nn+nh1),hr(nn+nh2),nhv,istrt,ad,sigl(1,l),
     &                 bei,xlamd,ha,.true.,isw)
           weng(l)   = estore

c          Compute Green & Almansi strains

           call fstrain(f(1,1,l),finv(1,l), egreen,ealmansi)

           do i = 1,6
             epsl(i,l) = ealmansi(i)
           end do ! i

c          Compute principal stretches

           call pstr3d(bei, bpr)

c          Average stresses and stretches for printing

           do i = 1,3
             bpra(i) = bpra(i) + 0.5d0*qfact*log(bpr(i))
             xxm(i)  = xxm(i)  + qfact*xu(i,l)
           end do ! i
           do i = 1,6
             sigm(i) = sigm(i) + qfact*sigl(i,l)
             ebig(i) = ebig(i) + qfact*egreen(i)
             esml(i) = esml(i) + qfact*ealmansi(i)
           end do ! i
           sigm(10) = sigm(10) + qfact*sigl(10,l)
           epp      = epp      + qfact*sigl( 9,l)
           dtheta   = dtheta   + qfact*theta(1,l)

           nn = nn + nhv
         end do ! l

c        Output stresses

         if (isw .eq. 4) then


           mct = mct - 2
           if(mct.le.0) then
             write(iow,2001) o,head
             if(ior.lt.0) write(*,2001) o,head
             mct = 50
           endif

c          Output stress at quadrature point

           if(qoutfl) then

             do i = 1,3
               bpra(i) = 0.5d0*log(bpr(i))
               xxm(i)  = xu(i,l)
             end do ! i

             call pstr3d(sigl(1,l),sigm(7))
             thlog = log(abs(theta(1,l)))
             write(iow,2002) n,ma,(sigl(i,l),i=1,6),(sigm(i),i=7,9),
     &                       bpra,xxm,thlog,sigl(9,l),sigl(10,l),
     &                       (epsl(i,l),i=1,6)
             if(ior.lt.0) then
               write(*,2002) n,ma,(sigl(i,l),i=1,9),(sigm(i),i=7,9),
     &                       bpra,xxm,thlog,sigl(9,l),sigl(10,l),
     &                       (epsl(i,l),i=1,6)
             endif

c          Output averaged stress only

           else

             call pstr3d(sigm,sigm(7))

             thlog = log(abs(dtheta))
             write(iow,2002) n,ma,(sigm(i),i=1,9),bpra,
     &                       xxm,thlog,epp,sigm(10),esml
             if(ior.lt.0) then
               write(*,2002) n,ma,(sigm(i),i=1,9),bpra,
     &                       xxm,thlog,epp,sigm(10),esml
             endif

           endif

c        Project stresses onto nodes

         elseif(isw.eq.8) then

           call slcn3d(sigl,epsl, r,s, nel)

c        Compute fracture indices

         elseif(isw.eq.16) then

           call pfrac3f(f,detf,sigl,weng, dvol, r, ndf,ndm,3)

         endif ! isw = 4 or 8 or 16

       endif ! isw = 3 or 6 or 4 or 8 or 14 or 16

      endif ! isw tests

c     Formats for input-output

2001  format(a1,20a4//5x,'Element Stresses'//
     &   '  Elmt  Matl  11-stress  22-stress  33-stress',
     &   '  12-stress  23-stress  13-stress'/4x,'Cauchy  ',
     &   '   1-stress   2-stress   3-stress',
     &   '  log(lam1)  log(lam2)  log(lam3)'/12x,
     &    '   1-coord    2-coord    3-coord',
     &   '       log-J     eff-ep      Yield'/4x,'Almansi '
     &   '  11-Strain  22-Strain  33-Strain  12-Strain',
     &   '  23-Strain  31-Strain')

2002  format(/2i6,1p,6e11.3/(12x,1p,6e11.3:))

      end
