c$Id:$
      subroutine fld3d2u(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    17/01/2008
c       1. Revise augment: augfp -> augf                    14/04/2009
c       2. Add prints for Almansi strains                   19/10/2009
c       3. Increase arrays to store 125 node brick          20/12/2010
c       4. Add 'l' to modlfd call                           05/01/2012
c       5. Add average of density for multiscale            10/05/2012
c       6. Add eps to slcn3du call                          01/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose  : Finite uniform deformation gradient model: u-theta

c      Reference: J. Bonet and P. Bhargava, "A uniform deformation
c                 gradient hexahedron element with artificial hourglass
c                 control", IJNME, Vol. 38, pp 2809-2828 (1995)

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
      include  'counts.h'
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
      include  'pointer.h'
      include  'prstrs.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   bflg, setvar,palloc
      integer   ndf,ndm,nst,isw,i,ii,j,jj,l,nhv,nn, istrt
      integer   k,m

      real*8    xlamd,ha, d1, epp, thlog, dtheta, xsj,ta
      real*8    press,dpress,dmass,cfac,lfac,cmshp,lmshp,bdb,vol0,vol1
      real*8    pvol0, coni,conj
      real*8    d(*),  ul(ndf,nen,*),  xl(ndm,*), s(nst,*),r(ndf,*)

      integer   is(8)
      real*8    sg(4,8),fi(9,2),finv(9),df(9),theta(2)
      real*8    bbd(3,7),  bei(6), ad(7,7,5), dd(7,7), sigm(10)
      real*8    dvol0(8), fdet(8),detf(2),xxm(3),egreen(6),ealmansi(6)
      real*8    sigl(10),shp0(3,8),shp(3,8),shptil(3,8),weng(8)
      real*8    xr(3,8),ur(3,8),xu(3,8),ru(3,8)
      real*8    acc(3),bpra(3),body(3), bf(3),bt(3,3)

c     Stabilized modes
      real*8     alpha, kmax
      real*8     khg(8,8),hgmode(8,4),hgshp(4,8),xhg(3)

      save

c     TEMPORARY TEMPERATURE

      data      ta    / 0.0d0 /

c     Hourglass normal modes:
c       (1) xi_1*xi_2; (2) xi_2*xi_3; (3) xi_3*xi*1; (4) xi_1*xi_2*xi_3

      data       hgmode/ 1.d0,-1.d0, 1.d0,-1.d0, 1.d0,-1.d0, 1.d0,-1.d0,
     &                   1.d0, 1.d0,-1.d0,-1.d0,-1.d0,-1.d0, 1.d0, 1.d0,
     &                   1.d0,-1.d0,-1.d0, 1.d0,-1.d0, 1.d0, 1.d0,-1.d0,
     &                  -1.d0, 1.d0,-1.d0, 1.d0, 1.d0,-1.d0, 1.d0,-1.d0/

c     Input data adjustments

      if(isw.eq.1) then

        nh1 = nh1 + 2
        nh3 = nh3 + 1

c     Augmented Lagrangian update for nested iteration

      elseif(isw.eq.10) then

        hr(nh2) = hr(nh2) + augf*d(185)*hr(nh2+1)

c     Compute tangent stiffness and residual force vector

      elseif(isw.eq.3 .or. isw.eq.4  .or. isw.eq.6 .or.
     &       isw.eq.8 .or. isw.eq.14 .or. isw.eq.16) then

       estore = 0.0d0

c      Set quadrature and order

       if(nel.eq.8) then
         if(nint(d(182)).gt.0) then
           call int3dn(nel,lint,sg)
         else
           l = 2
           call int3d(l,lint,sg)
         endif
         npm = 1
         nvn = 8
       else
         write(ilg,*) ' FLD3D2U ERROR: Must be 8-node element'
         write(iow,*) ' FLD3D2U ERROR: Must be 8-node element'
         call plstop()
       endif

c      Compute shape functions and derivatives at time t_n

       do l = 1,lint
         call shp3d(sg(1,l),xsj,shp3(1,1,l),xl,ndm,nel)
         dvol0(l) = xsj*sg(4,l)
         do i = 1,3
           xr(i,l) = 0.0d0
           ur(i,l) = 0.0d0
           do j = 1,nel
             xr(i,l) = xr(i,l) + shp3(4,j,l)*xl(i,j)
             ur(i,l) = ur(i,l) + shp3(4,j,l)*ul(i,j,1)
           end do ! j
         end do ! j

       end do ! l

c      Set number of history terms / quadradure point

       nhv   = nint(d(15))
       istrt = nint(d(84))

c      MECHANICAL ELEMENT

       if(isw.eq.3 .or. isw.eq. 6 .or. isw.eq.14) then

c        Compute f, finv, df and det(f) at t_n and t-n+1

         if(isw.eq.14) then
           call pfinit(fi,df,finv,detf, 1)
         else
           call kine3d2u(shp3,dvol0,ul,fi,finv,df,fdet,detf,theta,
     &                   shp0,shp,shptil,vol0,ndf,nel,nen,lint)
         endif

c        Compute volume at t_n+1

         vol1 = vol0*theta(1)

c        Set reference and t_n+1 coordinates

         do i = 1,3
           xref(i) = xr(i,l)
           xcur(i) = xr(i,l) + ur(i,l)
         end do ! i

c        Set augmented multipler

         xlamd = hr(nh2)

c        Compute Cauchy stresses and spatial tangent tensor at t-n+1

         nn = 2
         l  = 1
         call modlfd(l,d,fi,finv,df,theta,ta,hr(nn+nh1),hr(nn+nh2),
     &               nhv,istrt,ad,sigl,bei,xlamd,ha,.true.,isw)

c        Return if initialization step

         if(isw.eq.14) return

c        Set augmented function

         hr(nh2+1) = ha

c        Compute tangent and/or residual arrays

         if(isw.eq.3 .or. isw.eq.6) then

c          Compute assembly pointer

           is(1) = 0
           do j = 2,nel
             is(j) = is(j-1) + ndf
           end do ! j

           call sbodyf(d, body) ! Compute body force values

c          Compute pressure

           press = one3*(sigl(1) + sigl(2) + sigl(3))

c          Set flag for angular spin body force computation

           bflg  = d(4).gt.0.0d0 .and. d(65).gt.0.0d0

c          Compute inertia factors

           if((d(7).ge.0.0 .or. d(183).ne.0.0d0) .and.
     &                (ndfo(1).gt.0 .or. shflg)) then
             cfac = d(7)
             lfac = 1.d0 - cfac
           else
             cfac  = 0.0d0
             lfac  = 0.0d0
           endif ! d(7) test

c          Compute integrated residual and tangent arrays

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

c            Compute acceleration

             dmass = d(4)*dvol0(l)
             cmshp = dmass*cfac      ! Consistent mass factor
             lmshp = dmass*lfac      ! Lumped     mass factor
             do i = 1,3
               acc(i) = 0.0d0
               do j = 1,nel
                 acc(i) = acc(i) + shp3(4,j,l)*ul(i,j,5)
               end do ! j
               acc(i) = acc(i)*cmshp
               xxm(i) = (body(i) + bf(i))*dvol0(l)
             end do ! i

c            Compute residual for body, spin and inertia effects

             do j = 1,nel
               do i = 1,3
                 r(i,j)  = r(i,j) + shp3(4,j,l)*(xxm(i)
     &                   - (acc(i) + ul(i,j,5)*lmshp))
               end do ! i
             end do ! j

c            Geometric, spin and Inertial tangent

             if(isw.eq.3) then

c              Geometric stiffness of mixed pressure term

               if(gflag) then

                 pvol0 = press*fdet(l)*dvol0(l)
                 do i = 1,nel
                   do j = 1,nel
                     do k = 1,3
                       coni = shp3(k,i,l)*pvol0
                       conj = shp3(k,j,l)*pvol0
                       do m = 1,3
                         s(is(i)+k,is(j)+m) = s(is(i)+k,is(j)+m)
     &                                      + coni*shp3(m,j,l)
     &                                      - shp3(m,i,l)*conj
                       end do ! m
                     end do ! k
                   end do ! j
                 end do ! i

               endif ! gflag

c              Angular velocity tangent

               if(bflg) then

                 do i = 1,nel
                   do jj = 1,3
                     do ii = 1,3
                       bdb = shp3(4,i,l)*bt(ii,jj)
                       do j = 1,nel
                         s(is(i)+ii,is(j)+jj) = s(is(i)+ii,is(j)+jj)
     &                                        + bdb*shp3(4,j,l)
                       end do ! j
                     end do ! ii
                   end do ! jj
                 end do ! i
               endif ! bflg

c              Inertial tangent

               if(ctan(3).ne.0.0d0) then

                 dmass = ctan(3)*dmass
                 do i = 1,nel
                   lmshp = shp3(4,i,l)*dmass
                   cmshp = lmshp*cfac             ! Consistent mass
                   lmshp = lmshp*lfac             ! Lumped     mass

c                  Diagonal mass contribution

                   do jj = 1,3
                     s(is(i)+jj,is(i)+jj) = s(is(i)+jj,is(i)+jj) + lmshp
                   end do ! jj

c                  Consistent mass contribution

                   do j = 1,nel
                     do jj = 1,3
                       s(is(i)+jj,is(j)+jj) = s(is(i)+jj,is(j)+jj)
     &                                      + cmshp*shp3(4,j,l)
                     end do ! jj
                   end do ! j
                 end do ! i
               endif ! ctan(3) > 0

             endif ! isw.eq.3
           end do ! l

c          Store time history data and multiply by current volume

           do j = 1,6
             tt(j)   = sigl(j)
             sigl(j) = sigl(j)*vol1
           end do ! j
           dpress = press*vol1

c          Compute deviator stress times current volume

           sigm(1)  = sigl(1) - dpress
           sigm(2)  = sigl(2) - dpress
           sigm(3)  = sigl(3) - dpress
           sigm(4)  = sigl(4)
           sigm(5)  = sigl(5)
           sigm(6)  = sigl(6)

c          Compute residual

           do j = 1,nel

             ru(1,j) = shp(1,j)*sigm(1)
     &               + shp(2,j)*sigm(4)
     &               + shp(3,j)*sigm(6)

             ru(2,j) = shp(1,j)*sigm(4)
     &               + shp(2,j)*sigm(2)
     &               + shp(3,j)*sigm(5)

             ru(3,j) = shp(1,j)*sigm(6)
     &               + shp(2,j)*sigm(5)
     &               + shp(3,j)*sigm(3)

             r(1,j) = r(1,j) - ru(1,j) - shptil(1,j)*dpress

             r(2,j) = r(2,j) - ru(2,j) - shptil(2,j)*dpress

             r(3,j) = r(3,j) - ru(3,j) - shptil(3,j)*dpress
           end do ! j

c          Compute tangent stiffness matrix

           if(isw.eq.3) then

c            Part 1: Geometric tangent matrix (uniform part)

             if(gflag) then

               do i = 1,nel
                 do j = 1,nel
                   bdb = (ru(1,i)*shp(1,j)
     &                 +  ru(2,i)*shp(2,j)
     &                 +  ru(3,i)*shp(3,j))*ctan(1)
                   do jj = 1,3
                     s(is(i)+jj,is(j)+jj) = s(is(i)+jj,is(j)+jj) + bdb
                   end do ! jj
                 end do ! j
               end do ! i

c              Modify tangent moduli for geometric stress factors

               dpress   = 0.0d0
               call dmatdx(ad,sigl,press,dpress)

             endif ! gflag

c            Part 2: Material tangent matrix

c            Multiply tangent moduli by volume element

             d1 = vol1*ctan(1)
             do j = 1,7
               do i = 1,7
                 dd(i,j) = ad(i,j,1)*d1
               end do ! i
             end do ! j

c            Compute stiffness matrix

             do i = 1,nel

               do jj = 1,7
                 bbd(1,jj) =   shp(1,i)*dd(1,jj)
     &                     +   shp(2,i)*dd(4,jj)
     &                     +   shp(3,i)*dd(6,jj)
     &                     + shptil(1,i)*dd(7,jj)

                 bbd(2,jj) =   shp(2,i)*dd(2,jj)
     &                     +   shp(1,i)*dd(4,jj)
     &                     +   shp(3,i)*dd(5,jj)
     &                     + shptil(2,i)*dd(7,jj)

                 bbd(3,jj) =   shp(3,i)*dd(3,jj)
     &                     +   shp(2,i)*dd(5,jj)
     &                     +   shp(1,i)*dd(6,jj)
     &                     + shptil(3,i)*dd(7,jj)
               end do ! jj

               do j = 1,nel
                 do jj = 1,3

                   s(is(i)+jj,is(j)+1) = s(is(i)+jj,is(j)+1)
     &                                 + bbd(jj,1)*shp(1,j)
     &                                 + bbd(jj,4)*shp(2,j)
     &                                 + bbd(jj,6)*shp(3,j)
     &                                 + bbd(jj,7)*shptil(1,j)

                   s(is(i)+jj,is(j)+2) = s(is(i)+jj,is(j)+2)
     &                                 + bbd(jj,2)*shp(2,j)
     &                                 + bbd(jj,4)*shp(1,j)
     &                                 + bbd(jj,5)*shp(3,j)
     &                                 + bbd(jj,7)*shptil(2,j)

                   s(is(i)+jj,is(j)+3) = s(is(i)+jj,is(j)+3)
     &                                 + bbd(jj,3)*shp(3,j)
     &                                 + bbd(jj,5)*shp(2,j)
     &                                 + bbd(jj,6)*shp(1,j)
     &                                 + bbd(jj,7)*shptil(3,j)
                 end do ! jj
               end do ! j
             end do ! i
           endif ! isw = 3

c          Stabilize tangent/residual (Ref: Bonet and Bhargava)

           if(nstep.le.1 .and. niter.eq.0) then

c            For block use 0.01, for tube use 0.003

             alpha = d(60)
             if(alpha.eq.0.0d0) then
               alpha = 0.0001d0
             endif
             i      = ndf*nel
             setvar = palloc(111,'TEMP1',2*i*i+5*i, 2)
             call seigen(s,i,nst,hr(np(111)),kmax)
             setvar = palloc(111,'TEMP1',        0, 2)
             alpha = alpha*kmax
c            write(*,*) ' Alpha = ',alpha
             hr(nh3) = alpha
           else
             alpha = hr(nh3)
           end if

           do i = 1,4
             do j = 1,8
               hgshp(i,j) = hgmode(j,i)
             end do ! j
             do k = 1,3
               xhg(k) = 0.0d0
               do m = 1,8
                 xhg(k) = xhg(k) + xl(k,m)*hgmode(m,i)
               end do ! m
             end do ! k
             do j = 1,8
               do k = 1,3
                 hgshp(i,j) = hgshp(i,j) - shp0(k,j)*xhg(k)
               end do ! k
             end do ! j
           end do ! i

c          Compute hourglass stiffness matrix khg

           do j = 1,8
             do i = 1,8
               khg(i,j) = 0.0d0
               do k = 1,4
                 khg(i,j) = khg(i,j) + hgshp(k,i)*hgshp(k,j)
               end do ! k
               khg(i,j) = alpha*khg(i,j)
             end do ! i
           end do ! j

c          Add hourglass force and stiffness contribution

           do j = 1,8
             do i = 1,8
               do m = 1,3
                 r(m,i)             = r(m,i) - khg(i,j)*ul(m,j,1)
                 s(is(i)+m,is(j)+m) = s(is(i)+m,is(j)+m) + khg(i,j)
               end do ! m
             end do ! i
           end do ! j

         endif ! isw = 3 or 6

c      Output stresses and computre fracture force

       elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16) then

c        Compute current geometry

         do i = 1,3
          do j = 1,nel
           xu(i,j) = xl(i,j) + ul(i,j,1)
          end do ! j
         end do ! i

c        Compute f, finv, df and det(fei) at conf t-n+1

         call kine3d2u(shp3,dvol0,ul,fi,finv,df,fdet,detf,theta,
     &                 shp0,shp,shptil,vol0,ndf,nel,nen,lint)

c        Compute Cauchy stresses and spatial tangent tensor at t-n+1

         do i = 1,3
           xref(i) = xr(i,l)
           xcur(i) = xr(i,l) + ur(i,l)
         end do ! i

c        Set augmented function

         xlamd = hr(nh2)

         nn = 2
         l  = 1
         call modlfd(l,d,fi,finv,df,theta,ta,hr(nn+nh1),hr(nn+nh2),
     &               nhv,istrt,ad,sigl,bei,xlamd,ha,.true.,isw)
         weng(1)   = estore

c        Compute Green & Almansi strains

         call fstrain(fi,finv, egreen,ealmansi)

c        Compute principal stretches

         call pstr3d(bei, bpr)

c        Average stresses and stretches for printing

         do i = 1,3
           bpra(i) = 0.5d0*log(bpr(i))
           xxm(i)  = xu(i,l)
         end do ! i
         do i = 1,6
           sigm(i) = sigl(i)
         end do ! i
         sigm(10) = sigl(10)
         epp      = sigl( 9)
         dtheta   = theta(1)

c        Output stresses

         if (isw .eq. 4) then

           call pstr3d(sigm,sigm(7))

           mct = mct - 2
           if(mct.le.0) then
             write(iow,2001) o,head
             if(ior.lt.0) write(*,2001) o,head
             mct = 50
           endif

           thlog = log(abs(dtheta))
           write(iow,2002) n,ma,(sigm(i),i=1,9),bpra,
     &                     xxm,thlog,epp,sigm(10),ealmansi
           if(ior.lt.0) then
             write(*,2002) n,ma,(sigm(i),i=1,9),bpra,
     &                     xxm,thlog,epp,sigm(10),ealmansi
           endif
         elseif(isw.eq.8) then

c          Project stresses onto nodes

           call slcn3du(sigl,ealmansi,shp3,dvol0, r,s, lint,nel,125)

c        Compute fracture indices

         elseif(isw.eq.16) then

           call pfrac3f(fi,detf,sigl,weng, dvol0, r, ndf,ndm,3)

         endif ! isw = 4 or 8 or 16

       endif ! isw = 3 or 6 or 4 or 8 or 14 or 16

      endif ! isw tests

c     Formats for input-output

2001  format(a1,20a4//5x,'Element Stresses'//
     &   '  Elmt  Matl  11-stress  22-stress  33-stress',
     &   '  12-stress  23-stress  13-stress'/4x,'Cauchy  ',
     &   '   1-stress   2-stress   3-stress',
     &   '  log(lam1)  log(lam2)  log(lam3)'/12x,
     &   '    1-coord    2-coord    3-coord',
     &   '       log-J     eff-ep      Yield'/4x,'Almansi '
     &   '  11-Strain  22-Strain  33-Strain  12-Strain',
     &   '  23-Strain  31-Strain')

2002  format(/2i6,1p,6e11.3/(12x,1p,6e11.3:))

      end
