c$Id:$
      subroutine frans2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c     2. Set d(21) to d(1) in isw = 1 to get E              14/01/2010
c     3. Multiply u geometric stiffness by tolerance        10/10/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Two dimensional Euler-Bernoulli frame element

c     N.B. ELASTIC ONLY

c     Control Information:

c       ndm  - Spatial dimension of problem       = 2
c       ndf  - Number degree-of-freedoms at node >= 3
c              ( 1 = u_1 ; 2 = u_2 ; 3 = theta )
c       nen  - Two node element (nel = 2)        >= 2
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include   'bdata.h'
      include   'cdata.h'
      include   'eldata.h'
      include   'eltran.h'
      include   'evdata.h'
      include   'iofile.h'
      include   'part0.h'
      include   'rdata.h'

      integer    i,j,ii,jj,i1,j1,i2,j2,ndf,ndm,nst,isw
      real*8     cs,sn,le,xx,yy,xn,xm,vv,eps,chi,gam,EA,EI,RA,RI,ctan3
      real*8     d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),r(ndf,*)
      real*8     sm(6,6),rm(6), fi(3,2)

      save

c     Compute element parameters

      if(isw.gt.2) then
        cs = xl(1,2) - xl(1,1)
        sn = xl(2,2) - xl(2,1)
        le = sqrt(cs*cs+sn*sn)
        cs = cs/le
        sn = sn/le
      endif

c     Material parameter adjustment

      if(isw.eq.1) then

c       Set effective modulus to E

        d(21) = d(1)

c     Compute elment stiffness/residual arrays

      elseif(isw.eq.3 .or. isw.eq.6 .or. isw.eq.8) then

c       Set body loads

        call fbody2d(d,xl,ul, r,s, ndm,ndf,nst, isw)

c       Compute momentum residual and tangent matrix

        EA = d(21)*d(32)
        EI = d(21)*d(33)
        RA = d(4)*d(32)
        if(d(8).gt.0.0d0) then
          RI = d(4)*d(33)*d(8)
        else
          RI = 0.0d0
        endif
        call beam2d(s ,EA,EI,le,cs,sn,nst,ndf)

        if(isw.eq.3 .or. isw.eq. 6) then
          call massf2(sm,rm,d(7),RA,RI,le,cs,sn,6,3)

          if(ndfo(1).gt.0 .or. shflg) then
            ctan3 = ctan(3)
          else
            ctan3 = 0.0d0
          endif

c         Stress and mass modification to residual and stiffness

          i1 = 0
          i2 = 0
          do ii = 1,2
            do i = 1,ndf
              j1 = 0
              j2 = 0
              do jj = 1,2
                do j = 1,ndf

c                 Residual

                  r(i,ii)      = r(i,ii) -  s(i+i1,j+j1)*ul(j,jj,1)
     &                                   - sm(i+i2,j+j2)*ul(j,jj,5)

c                 Tangent

                  s(i+i1,j+j1) =  s(i+i1,j+j1)*ctan(1)
     &                         + sm(i+i2,j+j2)*ctan3

                end do ! j
                j1 = j1 + ndf
                j2 = j2 + 3
              end do ! jj

            end do ! i
            i1 = i1 + ndf
            i2 = i2 + 3
          end do ! ii

c       Stress projection

        elseif(isw.eq.8) then

          i1 = 0
          do ii = 1,2
            do i = 1,ndf
              fi(i,ii) = 0.0d0
              j1 = 0
              do jj = 1,2
                do j = 1,ndf
                  fi(i,ii) = fi(i,ii) - s(i+i1,j+j1)*ul(j,jj,1)
                end do ! jj
                j1 = j1 + ndf
              end do ! jj
            end do ! i
            i1 = i1 + ndf
          end do ! ii

          do i = 1,ndf
            fi(i,2) = -fi(i,2)
          end do ! i
          do ii = 1,2
            rm(1)    =  fi(1,ii)*cs + fi(2,ii)*sn
            fi(2,ii) = -fi(1,ii)*sn + fi(2,ii)*cs
            fi(1,ii) =  rm(1)
          end do ! ii

          do i = 1,nst
            do j = 1,nst
              s(j,i) = 0.0d0
            end do ! j
          end do ! i

          call frcn2d(fi,r,s)

        endif

c     Compute element output quantities

      elseif(isw.eq.4) then

        xx  = 0.5d0*(xl(1,1)+xl(1,2))
        yy  = 0.5d0*(xl(2,1)+xl(2,2))
        eps = (cs*(ul(1,2,1)-ul(1,1,1))
     &       + sn*(ul(2,2,1)-ul(2,1,1)))/le
        chi = (ul(3,2,1)-ul(3,1,1))/le
        gam =-12.d0*(-sn*(ul(1,1,1)-ul(1,2,1))
     &              + cs*(ul(2,1,1)-ul(2,2,1)))/le**3
     &            - 6.d0*(ul(3,1,1)+ul(3,2,1))/le/le
        xn  = d(21)*d(32)*eps
        xm  = d(21)*d(33)*chi
        vv  = d(21)*d(33)*gam

        mct = mct - 1
        if(mct.le.0) then
          write(iow,2000) o,head
          if(ior.lt.0) write(*,2000) o,head
          mct = 50
        endif

        write(iow,2001) n,ma,xx,yy,xn,xm,vv,eps,chi,gam
        if(ior.lt.0) write(*,2001) n,ma,xx,yy,xn,xm,vv,eps,chi,gam

c     Compute element mass or geometric stiffness arrays

      elseif(isw.eq.5) then

c       Mass term

        if(imtyp.eq.1) then

          RA = d(4)*d(32)
          if(d(8).gt.0.0d0) then
            RI = d(4)*d(33)*d(8)
          else
            RI = 0.0d0
          endif
          call massf2(sm,rm,d(7),RA,RI,le,cs,sn,nst,ndf)

          i1 = 0
          i2 = 0
          do ii = 1,2
            do i = 1,ndf
              j1 = 0
              j2 = 0
              do jj = 1,2
                do j = 1,ndf
                  s(i+i1,j+j1) = sm(i+i2,j+j2)
                end do ! j
                j1 = j1 + ndf
                j2 = j2 + 3
              end do ! jj
              r(i,ii) = rm(i+i2)
            end do ! i
            i1 = i1 + ndf
            i2 = i2 + 3
          end do ! ii

c       Geometric stiffness term

        elseif(imtyp.eq.2) then

          EA  = d(21)*d(32)
          eps = (cs*(ul(1,2,1)-ul(1,1,1))
     &        +  sn*(ul(2,2,1)-ul(2,1,1)))/le
          xn  = -EA*eps*ctan(1)

          s(1    ,1    ) =  xn/le*1.d0-2  ! Non-zero to avoid failure
          s(2    ,2    ) =  1.2d0*xn/le
          s(2    ,3    ) =  0.1d0*xn
          s(3    ,2    ) =  s(2,3)
          s(3    ,3    ) =  2.d0*xn*le/15.d0

          s(1    ,1+ndf) = -s(1,1)
          s(2    ,2+ndf) = -s(2,2)
          s(2    ,3+ndf) =  s(2,3)
          s(3    ,2+ndf) = -s(2,3)
          s(3    ,3+ndf) = -xn*le/30.d0

          s(1+ndf,1    ) = -s(1,1)
          s(2+ndf,2    ) = -s(2,2)
          s(2+ndf,3    ) = -s(2,3)
          s(3+ndf,2    ) =  s(2,3)
          s(3+ndf,3    ) =  s(3,3+ndf)

          s(1+ndf,1+ndf) =  s(1,1)
          s(2+ndf,2+ndf) =  s(2,2)
          s(2+ndf,3+ndf) = -s(2,3)
          s(3+ndf,2+ndf) = -s(2,3)
          s(3+ndf,3+ndf) =  s(3,3)

c         Transform to global coordinates

          call bm2trn (s,cs,sn,nst,ndf,1)

        endif
      end if

c     Formats

2000  format(a1,20a4//5x,'Beam Element Stresses'//
     &     '    Elmt  Matl     x-coord     y-coord     ',
     &     ' Force      Moment       Shear'/
     & 43x,'Strain    Curvature      Gamma')

2001  format(i8,i6,1p,5e12.4/38x,1p,3e12.4)

      end

      subroutine beam2d(s,ea,ei,le,cs,sn,nst,ndf)

c     Stiffness for frame element

      implicit none

      integer nst,ndf,i,j,k
      real*8  ea,ei,le,cs,sn,t, s(nst,nst)

      i = ndf + 1
      j = ndf + 2
      k = ndf + 3

      t = ea/le
      s(1,1) = t
      s(i,i) = t
      s(1,i) =-t
      s(i,1) =-t

      t = 12.d0*ei/le**3
      s(2,2) = t
      s(j,j) = t
      s(2,j) =-t
      s(j,2) =-t
      t = (ei+ei)/le
      s(3,3) = t + t
      s(k,k) = t + t
      s(3,k) = t
      s(k,3) = t
      t = 6.d0*ei/le**2
      s(2,3) = t
      s(3,2) = t
      s(2,k) = t
      s(k,2) = t
      s(3,j) =-t
      s(j,3) =-t
      s(j,k) =-t
      s(k,j) =-t

      call rotaf2(s,cs,sn,nst,ndf)

      end

      subroutine massf2(sm,rm,cfac,ra,ri,le,cs,sn,nst,ndf)

c     Frame mass matrix

      implicit   none

      include   'pconstant.h'

      integer    nst,ndf,i,j,l,ii(4)
      real*8     cfac,lfac,ra,ri,le,cs,sn,ta,ti,dv,s1
      real*8     rm(nst),sm(nst,nst),sg(2,4),bb(4)
      real*8     shpw(4,2), shpt(4,2)

      data       ii /2,3,5,6/

      ii(3)    = ndf + 2
      ii(4)    = ndf + 3

c     Lumped mass matrix

      ta   = 0.5d0*ra*le
      ti   = 0.5d0*ri*le
      rm(1) = ta
      rm(2) = ta
      rm(3) = ti
      rm(4) = ta
      rm(5) = ta
      rm(6) = ti

c     Consistent mass matrix

      do i = 1,6
        do j = 1,6
          sm(j,i) = 0.0d0
        end do ! j
      end do ! i

      sm(1,1) = two3*ta
      sm(1,4) = one3*ta
      sm(4,1) = sm(1,4)
      sm(4,4) = sm(1,1)

      call int1d (4,sg)
      do l = 1,4
        call shp1dh(sg(1,l),le, shpw, shpt)

c       Translational mass

        bb(1)  = shpw(4,1)
        bb(2)  = shpt(4,1)
        bb(3)  = shpw(4,2)
        bb(4)  = shpt(4,2)
        dv     = ta*sg(2,l)
        do i = 1,4
          s1 = bb(i)*dv
          do j = 1,4
            sm(ii(i),ii(j)) = sm(ii(i),ii(j)) + s1*bb(j)
          end do ! j
        end do ! i

c       Rotatory mass

        bb(1)  = shpw(1,1)
        bb(2)  = shpt(1,1)
        bb(3)  = shpw(1,2)
        bb(4)  = shpt(1,2)
        dv     = ti*sg(2,l)
        do i = 1,4
          s1 = bb(i)*dv
          do j = 1,4
            sm(ii(i),ii(j)) = sm(ii(i),ii(j)) + s1*bb(j)
          end do ! j
        end do ! i

      end do ! l

      lfac = 1.d0 - cfac
      do i = 1,6
        do j = 1,6
          sm(i,j) = cfac*sm(i,j)
        end do ! j
        sm(i,i) = sm(i,i) + lfac*rm(i)
      end do ! i

      call rotaf2(sm,cs,sn,nst,ndf)

      end

      subroutine rotaf2(s,cs,sn,nst,ndf)

c     Rotate arrays from local to global dof's

      implicit none

      integer nst,ndf,i,j,n
      real*8  cs,sn,t, s(nst,nst)

c     Check angle (perform rotation if necessary)

      if(cs.lt.0.99999999d0) then

        do i = 1,nst,ndf
          j = i + 1
          do n = 1,nst
            t      = s(n,i)*cs - s(n,j)*sn
            s(n,j) = s(n,i)*sn + s(n,j)*cs
            s(n,i) = t
          end do ! n
        end do ! i
        do i = 1,nst,ndf
          j = i + 1
          do n = 1,nst
            t      = cs*s(i,n) - sn*s(j,n)
            s(j,n) = sn*s(i,n) + cs*s(j,n)
            s(i,n) = t
          end do ! n
        end do ! i

      end if

      end
