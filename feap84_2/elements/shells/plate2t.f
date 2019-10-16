c$Id:$
      subroutine plate2t(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'pconstant.h'                 14/11/2006
c       2. Add rotational lumped masses                     05/03/2010
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Programmed by: Ferdinando Auricchio

c     Triangular plate: 3 dofs per node (w, theta-x, theta-y)
c                       2 bubble modes for rotations
c                       2 shear parameters

c     Mixed approach for shear stiffness.
c        Step 1: Condensation of bubble terms
c        Step 2: Condensation of shear terms

c     Three integration points are used.

c     Arguments:
c        d(*)       - specified parameter array
c        ul(ndf,*)  - local nodal solution values
c        xl(ndm,*)  - local nodal coordinate values
c        ix(*)      - node numbers
c        s(nst,nst) - finite element array (stiffness, mass, geometric
c                                           stiffness)
c        p(nst)     - finite element array (residual, lumped mass)
c        ndf        - number of degree of freedoms at node ( > or = 3 )
c        ndm        - spatial dimension of element         ( > or = 2 )
c        nst        - size of finite element arrays        ( > or = 9 )
c        isw        - solution option
c                   = 1: Input values and store in d(*) array
c                   = 2: Check mesh coordinate and connection inputs
c                        for errors
c                   = 3: Compute element residual (p) and stiffness (s)
c                   = 4: Output element results
c                   = 5: Compute mass (p,s)/geometric stiffness array(s)
c                   = 6: Compute element residual (p)
c                   = 8: Compute nodal projections

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Input parameters set as follows:

c         ndm = 2 (x,y cartesian coordinates at nodes)
c         ndf = 3 (w,theta-x,theta-y, at nodes)
c         nen = 3 nodes (counterclockwise around element)

c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'pconstant.h'
      include  'part0.h'
      include  'prld1.h'
      include  'ptdat2.h'
      include  'strnum.h'
      include  'rdata.h'
      include  'comblk.h'

      integer   ndf,ndm,nst,isw
      integer   i,j,k,l, lint, i1,j1, i2,j2, ns

      real*8    den, area, ar3, ar24, xx, yy, thk, thk3, psi, qt
      real*8    ctan1,ctan3

      logical   shflag
      integer   ix(*)
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),p(ndf,*)
      real*8    b(3),c(3),co(3),si(3),el(4,3),ss(9,9)
      real*8    bb(2,3,3), bbd(2,3,4), bs(3,2,3), hh(4,4), gamm(2)
      real*8    dbg(3,3), dsg(2,2), dd(3,3), dv(3), uhat(2), shear(2)
      real*8    alphm(3), epsc(3),eps(3,3),sigc(7),sig(3,3)

      save

c     Go to correct array processor

      go to (1,2,3,3,5,3,7,3,3,7,3), isw

1     return

c     Check element for errors

2     call cktris(ix,xl,bbd,ndm)
      return

c     Compute stiffness and internal force quantities

c     Compute geometric factors

3     call geompt(xl,ndm,co,si,b,c,area)

c     Set shear flag: true includes shear effects

      shflag = d(79).eq.0.0d0

c     Moduli multiplied by thickness factors

      ctan1 = ctan(1) + ctan(2)*d(78)
      if(ndfo(1).gt.0 .or. shflg) then
        ctan3 = ctan(3) + ctan(2)*d(77)
      else
        ctan3 = 0.0d0
      endif

      psi   = d(31)

      call dmatpl(d,psi,dbg,dsg,alphm)

      thk  = d(14)
      thk3 = thk**3/12.d0
      ar3  = thk3*area
      do j = 1,3
        do i = 1,3
          dd(i,j) = dbg(i,j)*ar3
        end do ! i
      end do ! j

c     Compute Bs-bar matrix

      den = 0.5d0*one3*area
      do i=1,3
        j  = mod(i,3) + 1
        k  = mod(j,3) + 1

        bs(1,1,i) =  b(i) * area
        bs(2,1,i) = (b(j)*co(j) - b(k)*co(k)) * den
        bs(3,1,i) = (b(j)*si(j) - b(k)*si(k) + 2.0d0) * den
        bs(1,2,i) =  c(i) * area
        bs(2,2,i) = (c(j)*co(j) - c(k)*co(k) - 2.0d0) * den
        bs(3,2,i) = (c(j)*si(j) - c(k)*si(k)) * den

      end do ! i

c     Numerical integration for bubble

      hh(3,3) = 0.0d0
      hh(3,4) = 0.0d0
      hh(4,4) = 0.0d0

      ar3     = 5.0625d0*one3

      do l = 1,3  ! { begin integration

c     Multiply by bending stiffness * quadrature weight

        do j = 1,3
          bbd(1,j,4) =  (c(l)*dd(2,j) + b(l)*dd(3,j))*ar3
          bbd(2,j,4) = -(b(l)*dd(1,j) + c(l)*dd(3,j))*ar3
        end do ! j

c     Bubble bending stiffness

        hh(3,3) = hh(3,3) + bbd(1,2,4)*c(l) + bbd(1,3,4)*b(l)
        hh(3,4) = hh(3,4) - bbd(1,1,4)*b(l) - bbd(1,3,4)*c(l)
        hh(4,4) = hh(4,4) - bbd(2,1,4)*b(l) - bbd(2,3,4)*c(l)

      end do  ! end integration }

c     Bubble shear mode and symmetric parts: N.B. Quadrature approx to
c                                            integral gives 0.5 not 0.45

      hh(1,4) =  0.5d0*area
      hh(2,3) = -hh(1,4)
      hh(3,2) =  hh(2,3)
      hh(4,1) =  hh(1,4)
      hh(4,3) =  hh(3,4)

c     Compute shear compliances * area

      if(shflag) then

        den = area/((dsg(1,1)*dsg(2,2) - dsg(1,2)*dsg(2,1))*thk)

        hh(1,1) = -dsg(2,2)*den
        hh(1,2) =  dsg(1,2)*den
        hh(2,1) =  dsg(2,1)*den
        hh(2,2) = -dsg(1,1)*den

c     Thin limit case

      else

        hh(1,1) = 0.0d0
        hh(1,2) = 0.0d0
        hh(2,1) = 0.0d0
        hh(2,2) = 0.0d0

      endif

      if ( isw.eq.3 .or. isw.eq.6 ) then


c     Set time varying pressure load

        if(int(d(74)).gt.0) then
          qt = d(10) + prldv(nint(d(74)))*d(71)
        else
          qt = d(10)*dm
        endif

c       Compute bending stiffness

        ar3 =  qt*area*one3
        ar24 = ar3*0.125d0

        do i = 1,3
          j          = mod(i,3) + 1
          k          = mod(j,3) + 1

          bb(1,1,i) =  0.0d0
          bb(1,2,i) = -c(i)
          bb(1,3,i) = -b(i)
          bb(2,1,i) =  b(i)
          bb(2,2,i) =  0.0d0
          bb(2,3,i) =  c(i)

          p(1,i)    =  ar3
          p(2,i)    =  ar24*(co(k)-co(j))
          p(3,i)    =  ar24*(si(k)-si(j))

        end do ! i

c       Strain-displacement times moduli

        do i = 1,3
          do j = 1,3
            bbd(1,j,i) = bb(1,2,i)*dd(2,j) + bb(1,3,i)*dd(3,j)
            bbd(2,j,i) = bb(2,1,i)*dd(1,j) + bb(2,3,i)*dd(3,j)
          end do ! j
        end do ! i

c       Bending stiffness (constant part)

        j1 = 2
        do j=1,3
          i1 = 2
          do i=1,j
            s(i1  ,j1  ) = bbd(1,2,i)*bb(1,2,j) + bbd(1,3,i)*bb(1,3,j)
            s(i1+1,j1  ) = bbd(2,2,i)*bb(1,2,j) + bbd(2,3,i)*bb(1,3,j)
            s(i1  ,j1+1) = bbd(1,1,i)*bb(2,1,j) + bbd(1,3,i)*bb(2,3,j)
            s(i1+1,j1+1) = bbd(2,1,i)*bb(2,1,j) + bbd(2,3,i)*bb(2,3,j)

c           Make symmetric part

            s(j1  ,i1  ) = s(i1  ,j1  )
            s(j1  ,i1+1) = s(i1+1,j1  )
            s(j1+1,i1  ) = s(i1  ,j1+1)
            s(j1+1,i1+1) = s(i1+1,j1+1)
            i1 = i1 + ndf
          end do ! i
          j1 = j1 + ndf
        end do ! j

c       Static condensation and load vector

        call sconpt(hh,bs,s,ndf,nst,1)

c       Compute stress residual

        do i = 1,nst
          j1 = 0
          do j = 1,3
            do k = 1,ndf
              p(i,1)    = p(i,1) - s(i,j1+k)*(ul(k,j,1)
     &                                + d(78)*ul(k,j,4))
              s(i,j1+k) = s(i,j1+k)*ctan1
            end do ! k
            j1 = j1 + ndf
          end do ! j
        end do ! i

c       Add inertial parts if necessary

        if(ctan3.ne.0.0d0) then

          call masspl(d,xl,ndm,3,9, hh,ss)

          i1 = 0
          i2 = 0
          do i = 1,3

            j1 = 0
            j2 = 0
            do j = 1,3

              do k = 1,3
                do l = 1,3
                  p(k,i)       = p(k,i) - ss(i2+k,j2+l)*(ul(l,j,5)
     &                                           + d(77)*ul(l,j,4))
                  s(i1+k,j1+l) = s(i1+k,j1+l) + ctan3*ss(i2+k,j2+l)
                end do ! l
              end do ! k

              j1 = j1 + ndf
              j2 = j2 + 3
            end do ! j

            i1 = i1 + ndf
            i2 = i2 + 3
          end do ! i

        endif

c       Tplot saves

        if(nsplts.gt.0) then

c         Compute  bending strains and moments

          epsc(1) =   b(1)*ul(3,1,1) + b(2)*ul(3,2,1) + b(3)*ul(3,3,1)
          epsc(2) = - c(1)*ul(2,1,1) - c(2)*ul(2,2,1) - c(3)*ul(2,3,1)
          epsc(3) = - b(1)*ul(2,1,1) - b(2)*ul(2,2,1) - b(3)*ul(2,3,1)
     &              + c(1)*ul(3,1,1) + c(2)*ul(3,2,1) + c(3)*ul(3,3,1)

          tt(1) = dd(1,1)*epsc(1) + dd(1,2)*epsc(2) + dd(1,3)*epsc(3)
          tt(2) = dd(2,1)*epsc(1) + dd(2,2)*epsc(2) + dd(2,3)*epsc(3)
          tt(3) = dd(3,1)*epsc(1) + dd(3,2)*epsc(2) + dd(3,3)*epsc(3)

c         Compute shear stress/strain

          gamm(1) = 0.0d0
          gamm(2) = 0.0d0
          do i = 1,3
            gamm(1) = gamm(1) - bs(1,1,i)*ul(1,i,1)
     &                        - bs(2,1,i)*ul(2,i,1)
     &                        - bs(3,1,i)*ul(3,i,1)
            gamm(2) = gamm(2) - bs(1,2,i)*ul(1,i,1)
     &                        - bs(2,2,i)*ul(2,i,1)
     &                        - bs(3,2,i)*ul(3,i,1)
          end do ! i

          tt(1) =  hh(1,1)*gamm(1) + hh(1,2)*gamm(2)
          tt(2) =  hh(2,1)*gamm(1) + hh(2,2)*gamm(2)

        endif

      else if ((isw.eq.4).or.(isw.eq.8).or.(isw.eq.11)) then

c       Constant part of bending/shear strains

        epsc(1) =   b(1)*ul(3,1,1) + b(2)*ul(3,2,1) + b(3)*ul(3,3,1)
        epsc(2) = - c(1)*ul(2,1,1) - c(2)*ul(2,2,1) - c(3)*ul(2,3,1)
        epsc(3) = - b(1)*ul(2,1,1) - b(2)*ul(2,2,1) - b(3)*ul(2,3,1)
     &            + c(1)*ul(3,1,1) + c(2)*ul(3,2,1) + c(3)*ul(3,3,1)

c       Recover bubble and shear parameters

        call sconpt(hh,bs,s,ndf,nst,0)

        gamm(1) = 0.0d0
        gamm(2) = 0.0d0
        do i = 1,3
          gamm(1) = gamm(1) - bs(1,1,i)*ul(1,i,1)
     &                      - bs(2,1,i)*ul(2,i,1)
     &                      - bs(3,1,i)*ul(3,i,1)
          gamm(2) = gamm(2) - bs(1,2,i)*ul(1,i,1)
     &                      - bs(2,2,i)*ul(2,i,1)
     &                      - bs(3,2,i)*ul(3,i,1)
        end do ! i

c       Compute shear stress/strain and bubble displacement

        shear(1) =  hh(1,1)*gamm(1) + hh(1,2)*gamm(2)
        shear(2) =  hh(2,1)*gamm(1) + hh(2,2)*gamm(2)

        dv(1)    = -hh(3,2)*shear(2)
        dv(2)    = -hh(4,1)*shear(1)

        uhat(1)  =  hh(3,3)*dv(1) + hh(3,4)*dv(2)
        uhat(2)  =  hh(4,3)*dv(1) + hh(4,4)*dv(2)

c       Get elastic material properties

        psi   = d(31)

        call dmatpl(d,psi,dbg,dsg,alphm)

        thk   = d(14)
        thk3  = thk**3/12.d0
        do i = 1,3 ! {
          do j = 1,3 ! {
            dd(i,j) = dbg(i,j)*thk3
          end do ! j   }
        end do ! i   }

        den = 1.d0/((dsg(1,1)*dsg(2,2) - dsg(1,2)*dsg(2,1))*thk)

        gamm(1) = ( dsg(2,2)*shear(1) - dsg(1,2)*shear(2))*den
        gamm(2) = (-dsg(1,2)*shear(1) + dsg(1,1)*shear(2))*den

        if(isw.eq.4) then

          xx    = (xl(1,1) + xl(1,2) + xl(1,3))/3.d0
          yy    = (xl(2,1) + xl(2,2) + xl(2,3))/3.d0

          sigc(1) = dd(1,1)*epsc(1) + dd(1,2)*epsc(2) + dd(1,3)*epsc(3)
          sigc(2) = dd(2,1)*epsc(1) + dd(2,2)*epsc(2) + dd(2,3)*epsc(3)
          sigc(3) = 0.0d0
          sigc(4) = dd(3,1)*epsc(1) + dd(3,2)*epsc(2) + dd(3,3)*epsc(3)

          call pstr2d(sigc(1),sigc(5))

c         Output moments and curvatures

          mct = mct -2
          if(mct.le.0) then
            write(iow,2001) o,head
            if(ior.lt.0) then
              write (*,2001) o,head
            endif
            mct = 50
          endif
          write(iow,2002) n,ma,(sigc(j),j=1,2),(sigc(j),j=4,6),
     &                    xx,yy,epsc,sigc(7),shear
          if(ior.lt.0) then
            write(*,2002) n,ma,(sigc(j),j=1,2),(sigc(j),j=4,6),
     &                    xx,yy,epsc,sigc(7),shear
          endif

        elseif(isw.eq.8 .or.isw.eq.11) then

c         Set quadrature

          l = -3
          call tint2d(l,lint,el)

          do l = 1, lint

c           Area weight for projection

            dv(l) = area

c           Bubble mode bending functions

            eps(1,l) = epsc(1) - 2.25d0*b(l)*uhat(2)
            eps(2,l) = epsc(2) + 2.25d0*c(l)*uhat(1)
            eps(3,l) = epsc(3) + 2.25d0*(b(l)*uhat(1) - c(l)*uhat(2))

c           Compute moments

            sig(1,l) = dd(1,1)*eps(1,l) + dd(1,2)*eps(2,l)
     &               + dd(1,3)*eps(3,l)
            sig(2,l) = dd(2,1)*eps(1,l) + dd(2,2)*eps(2,l)
     &               + dd(2,3)*eps(3,l)
            sig(3,l) = dd(3,1)*eps(1,l) + dd(3,2)*eps(2,l)
     &               + dd(3,3)*eps(3,l)

          end do ! l

c         Compute nodal stress values

          if(isw.eq.8) then

            call stcnpt(dv,el,lint,sig,eps,shear,gamm,p,s,p(nen+1,1))

c         Compute stress/energy errors

          elseif(isw.eq.11) then

             if(d(10).eq.0.0d0) then
               ns = 3
             else
               ns = 5
             endif

             call sterpt(dv,el,lint,sig,shear,eps,gamm,thk,
     &                       s,nen,d(50),ns)

          end if

        end if

      end if

      return

c     Compute consistent and lumped mass matrix

5     call masspl(d,xl,ndm,ndf,nst, p,s)
      return

c     Compute surface tractions

7     continue
      return

c     Formats for input-output

2001  format(a1,20a4//5x,'Element Moments'//'  Element Material',
     1   3x,'11-Moment',3x,'22-Moment',3x,'12-Moment',4x,
     2   '1-Moment',4x,'2-Moment'/2x,'1-Coord',2x,'2-Coord',3x,
     3   '11-Strain',3x,'22-Strain',3x,'12-Strain',12x,'Angle'/
     4   21x,'Shear-x  ',3x,'Shear-y')

2002  format(2i9,1p,5e12.3/0p,2f9.3,1p,3e12.3,0p,f18.2/18x,1p,2e12.3/1x)

      end

      subroutine stcnpt(dv,el,lint,sig,eps,shear,gamm,dt,st,ser)

      implicit  none

      include  'cdata.h'
      include  'prstrs.h'
      include  'strnum.h'

      integer   l, j, lint
      real*8    xg
      real*8    dt(*),st(nen,*),ser(*),sig(3,*),eps(3,*),dv(*)
      real*8    shear(2), gamm(2), el(4,*)

      save

c     Lumped projection routine

      do l = 1,lint
        do j = 1,3
          xg       = dv(l)*el(j,l)
          dt(j)    = dt(j)     + xg
          st(j,1 ) = st(j,1)   + sig(1,l)*xg
          st(j,2 ) = st(j,2)   + sig(2,l)*xg
          st(j,4 ) = st(j,4)   + sig(3,l)*xg
          st(j,5 ) = st(j,5)   + shear(1)*xg
          st(j,6 ) = st(j,6)   + shear(2)*xg
          st(j,7 ) = st(j,7)   + eps(1,l)*xg
          st(j,8 ) = st(j,8)   + eps(2,l)*xg
          st(j,10) = st(j,10) + eps(3,l)*xg
          st(j,11) = st(j,11) + gamm(1)*xg
          st(j,12) = st(j,12) + gamm(2)*xg

          ser(j)   = ser(j)   + erav*xg

        end do ! j
      end do ! l

      iste = 12

      end

      subroutine geompt(xl,ndm,co,si,b,c,area)

c     Compute geometric quantities and shape function derivatives
c     Linear shape functions ==> constant derivatives [b(i) and c(i)]

      implicit  none

      integer   i, ndm

      real*8    det,area
      real*8    b(3),c(3),co(3),si(3)
      real*8    xl(ndm,*)

      save

      b(1)  =  xl(2,2) - xl(2,3)
      b(2)  =  xl(2,3) - xl(2,1)
      b(3)  =  xl(2,1) - xl(2,2)

      c(1)  =  xl(1,3) - xl(1,2)
      c(2)  =  xl(1,1) - xl(1,3)
      c(3)  =  xl(1,2) - xl(1,1)

      det  = c(2)*b(1) - c(1)*b(2)
      area = 1.d0/det

      do i=1,3
        co(i) = -b(i)
        si(i) = -c(i)
        b(i)  =  b(i)*area
        c(i)  =  c(i)*area
      end do ! i

      area = det*0.5d0

      end

      subroutine sconpt(hh,bs,s,ndf,nst,isw)

c     Static condensation of internal dofs

      implicit  none

      include  'eltran.h'

      integer   isw, i, j, ii, jj, i1, j1, ndf,nst
      real*8    hh(4,4),bs(3,2,3), s(nst,nst), tem(2), det

      save

c     Solve bubble and shear modes

      det     = 1.d0/(hh(3,3)*hh(4,4) - hh(3,4)*hh(4,3))
      tem(1)  =  hh(4,4)*det
      hh(3,4) = -hh(3,4)*det
      hh(4,3) = -hh(4,3)*det
      hh(4,4) =  hh(3,3)*det
      hh(3,3) =  tem(1)

      hh(1,1) = hh(1,1) - hh(1,4)*hh(4,4)*hh(4,1)
      hh(2,1) = hh(2,1) - hh(2,3)*hh(3,4)*hh(4,1)
      hh(1,2) = hh(1,2) - hh(1,4)*hh(4,3)*hh(3,2)
      hh(2,2) = hh(2,2) - hh(2,3)*hh(3,3)*hh(3,2)

      det     =  1.0d0/(hh(1,1)*hh(2,2) - hh(1,2)*hh(2,1))
      tem(1)  =  hh(2,2)*det
      hh(1,2) = -hh(1,2)*det
      hh(2,1) = -hh(2,1)*det
      hh(2,2) =  hh(1,1)*det
      hh(1,1) =  tem(1)

c     Condense stiffness matrix

      if(isw.eq.1) then

        j1 = 0
        do jj = 1,3
          do j = 1,3

            tem(1) = hh(1,1)*bs(j,1,jj) + hh(1,2)*bs(j,2,jj)
            tem(2) = hh(2,1)*bs(j,1,jj) + hh(2,2)*bs(j,2,jj)

            i1 = 0
            do ii = 1,3
              do i = 1,3
                s(i+i1,j+j1) = s(i+i1,j+j1) - bs(i,1,ii)*tem(1)
     &                                      - bs(i,2,ii)*tem(2)
              end do ! i
            i1 = i1 + ndf
            end do ! ii
          end do ! j
          j1 = j1 + ndf
        end do ! jj

      end if

      end

      subroutine dmatpl(d,psi,dmg,dsg,alphg)

c     Rotate material arrays from principal to local element directions


c     Inputs: d        - Array with material properties
c             psi      - Angle of y1-axis (local) to 1-axis (principal)
c     Output: dmg(3,3) - Plane stress modulus matrix
c             dsg(2,2) - Transverse shear modulus matrix
c             alphg(3) - global thermal strain/temperature

c     Variables used in subroutine

c             qm(3,3)  - Transformation matrix for plane stresses

c             sig_glob = qm * sig_princ

c             dml(3,3) - Local (orthotropic ) plane modulus matrix
c             dmlqj(3) - intermediate matrix for triple product
c             alphl(3) - local  thermal strain/temperature

c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer   i, j

      real*8    psi, si, co, s2, c2, cs

      real*8    d(*)
      real*8    dml(3,3), dmg(3,3), dsg(2,2), qm(3,3), dmlqj(3)
      real*8    alphg(3)

      save

c     Set up constants for transformation

      si = sin(psi)
      co = cos(psi)
      s2 = si*si
      c2 = co*co
      cs = co*si

c     Transformation matrix for plane stress

      qm(1,1) =  c2
      qm(1,2) =  s2
      qm(1,3) =  cs
      qm(2,1) =  s2
      qm(2,2) =  c2
      qm(2,3) = -cs
      qm(3,1) = -2.d0 * cs
      qm(3,2) =  2.d0 * cs
      qm(3,3) =  c2 - s2

c     Set up local (orthotropic) plane stress matrix

      dml(1,1) = d(21)
      dml(2,2) = d(22)

      dml(1,2) = d(24)
      dml(2,1) = d(24)

      dml(3,3) = d(27)

c     Convert plane stress local to global matrix

      do j = 1,3 ! {

        dmlqj(1) = dml(1,1)*qm(1,j) + dml(1,2)*qm(2,j)
        dmlqj(2) = dml(2,1)*qm(1,j) + dml(2,2)*qm(2,j)
        dmlqj(3) = dml(3,3)*qm(3,j)

        do i = 1,3 ! {
          dmg(i,j) = qm(1,i)*dmlqj(1) + qm(2,i)*dmlqj(2)
     +             + qm(3,i)*dmlqj(3)
        end do ! i   }

      end do ! j   }

c     Set up global shear matrix

      dsg(1,1) = (c2 * d(29) + s2 * d(28))*d(37)
      dsg(2,2) = (s2 * d(29) + c2 * d(28))*d(37)
      dsg(1,2) =  cs * ( d(29) - d(28) )  *d(37)
      dsg(2,1) = dsg(1,2)

c     Convert to global thermal stiffness vector

      alphg(1) = c2 * d(47) + s2 * d(48)
      alphg(2) = s2 * d(47) + c2 * d(48)
      alphg(3) = cs * ( d(47) - d(48) )

      end

      subroutine sterpt(dv,el,lint,sig,shear,eps,gamm,thk,
     &                  st,nen,error,ns)

      implicit  none

      include  'adapt1.h'
      include  'adapt2.h'
      include  'errind.h'

      integer   i,j,l,lint,nen,ns
      real*8    error,thk,thk2

      real*8    st(nen,*),sig(3,3), sigp(5), epsp(5), eps(3,3)
      real*8    shear(2), gamm(2), dv(3),el(3,3)

      save

c     Error projection routine

      vfem   = 0.0d0
      vproj  = 0.0d0
      verror = 0.0d0
      vener  = 0.0d0
      venere = 0.0d0
      heta   = 0.0d0
      thk2   = thk*thk

c     Quadrature loop

      do l = 1,lint

c       Compute projected stresses at integration points

        do i = 1,ns
          sigp(i)= 0.0d0
          epsp(i)= 0.0d0
        end do ! i
        do i = 1,3
          do j = 1,ns
            sigp(j) = sigp(j) + el(i,l)*st(i,j)
            epsp(j) = epsp(j) + el(i,l)*st(i,j+15)
          end do ! j
        end do ! i

c       Compute integral of stress squares for error indicator use

        heta = heta + dv(l)

        do i = 1,3
          vfem   = vfem   + sig(i,l)*sig(i,l)*dv(l)
          vproj  = vproj  + sigp(i)*sigp(i)*dv(l)
          verror = verror + ((sigp(i)-sig(i,l))**2)*dv(l)
          vener  = vener  + sig(i,l)*eps(i,l)*dv(l)
          venere = venere + (sigp(i)-sig(i,l))*(epsp(i)-eps(i,l))*dv(l)
        end do ! i

        do i = 4,ns
          vfem   = vfem   + shear(i-3)*shear(i-3)*thk2*dv(l)
          vproj  = vproj  + sigp(i)*sigp(i)*thk2*dv(l)
          verror = verror + ((sigp(i)-shear(i-3))**2)*thk2*dv(l)
          vener  = vener  + shear(i-3)*gamm(i-3)*dv(l)
          venere = venere + (sigp(i)-shear(i-3))
     &                    * (epsp(i)-gamm(i-3))*dv(l)
        end do ! i

      end do ! l

c     Set error indicators

      arsq   = arsq  + heta
      efem   = efem  + vfem
      eproj  = eproj + vproj
      eerror = eerror+ verror
      eener  = eener + vener
      eenere = eenere+ venere

      areai  = heta

c     Weight for triangles

      heta  =  error*sqrt(heta)

      end

      subroutine masspl(d,xl,ndm,ndf,nst, p,s)

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   ndm,ndf,nst, i,i1, j,k, l,lint
      real*8    d(*),xl(ndm,*), p(ndf,*),s(nst,nst)
      real*8    co(3),si(3),b(3),c(3),el(4,3),an(24)
      real*8    area, ar3,ar24, ccm,clm, den, h12

      save

c     Initialize mass matrix

      do i = 1,nst
        do j = 1,nst
          s(j,i) = 0.0d0
        end do ! j
      end do ! i

c     Compute consistent and lumped mass matrix

      l = -3
      call tint2d(l,lint,el)

c     Compute geometric factors

      call geompt(xl,ndm,co,si,b,c,area)

      ar3 = d(4)*d(14)*area
      ccm = d(7)*ar3
      clm = ar3 - ccm

c     Lumped mass matrix part

      i1 = 1
      do i=1,3
        p(1,i)   = ar3*one3
        p(2,i)   = 0.0d0
        p(3,i)   = 0.0d0
        s(i1,i1) = s(i1,i1) + clm*one3
        i1       = i1 + ndf
      end do ! i

c     Consistent mass matrix part

      do l=1,lint

        i1 = 1
        do i=1,3
          j        = mod(i,3) + 1
          k        = mod(j,3) + 1
          an(i1  ) = el(i,l)
          an(i1+1) = el(i,l)*(el(j,l)*(xl(2,j) - xl(2,i))
     &                      + el(k,l)*(xl(2,k) - xl(2,i)))*0.5d0
          an(i1+2) = el(i,l)*(el(j,l)*(xl(1,i) - xl(1,j))
     &                      + el(k,l)*(xl(1,i) - xl(1,k)))*0.5d0
          i1 = i1 + ndf
        end do ! i

        den = ccm*el(4,l)
        do j=1,3*ndf
          ar24 = an(j)*den
          do k=1,j
            s(j,k) = s(j,k) + ar24*an(k)
            s(k,j) = s(j,k)
          end do ! k
        end do ! j
      end do ! l

c     Add rotational mass

      h12 = d(8)*d(14)*d(14)/12.0d0
      do i = 1,3
        p(2,i) = p(1,i)*h12
        p(3,i) = p(2,i)
      end do ! i

      end
