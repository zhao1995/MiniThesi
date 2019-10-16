      subroutine mate3d10(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c
c     Purpose: calculate S and C for
c               a concrete material law
c
c     Inputs:
c         h1(nh)       - history array h1
c         d(md)        - local d-array
c         Eps          - strains
c         isw          - solution option from element
c
c     Input material parameter:
c         E            - Young's modulus
c         v            - Poisson's ratio
c         ipla         - plasticity on=1 /off=0
c         f_ctm        - tensile strength
c         f_cm         - compressive strength
c         G_f          - energy release during cracking
c         G_c          - energy release rate
c         gam_1        - fitting parameter
c         gam_2        - fitting parameter
c
c     Outputs:
c         md = 9       - number of used data for control of d-array
c         nh = 8       - number of history parameter at Gauss-Point
c         sig          - stresses
c         cmat         - tangent modulus
c         plout(10)    - plot data
c
c     Allocation of d-array:
c         d(1) = E     - Young's modulus
c         d(2) = v     - Poisson's ratio
c         d(3) = ipla  - plasticity on=1 /off=0
c         d(4) = f_ctm - tensile strength
c         d(5) = f_cm  - compressive strength
c         d(6) = G_f   - energy release during cracking
c         d(7) = G_c   - energy release rate
c         d(8) = gam_1 - fitting parameter
c         d(9) = gam_2 - fitting parameter
c
c
c     References:
c     [1] J. Schütt; Ein inelastisches 3D-Versagensmodell für Beton und
c         seine Finite-Element-Implementierung, Dissertation,
c         Jan Schütt, 2005, ISBN 978- 3-935322-08-9, Bericht 9,
c         Institut für Baustatik, Universität Karlsruhe (TH)
c
c         implemented in Matelib3d             SK 12/2010
c         modified to new version of Matelib3d WW 01/2011
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),Eps(*),Sig(*),Cmat(*),plout(10),E_p(6)
c
      if(isw.eq.1) then
c....   input data
        call mati3d10(d,md,nh)
c
      else
c....   Store history variables at time n (h1-array) in plastic variables
        do i = 1,6
           E_p(i) = h1(i)
        end do
        a = h1(7)
        b = h1(8)
c
c       compute elastic-plastic tangent modulus and stresses
        call matm3d10(d,sig,cmat,Eps,E_p,a,b,lgp,dvp)
c
c....   store history variables at time n+1 in h2-array
        do i = 1,6
           h2(i) = E_p(i)
        end do
        h2(7) = a
        h2(8) = b

        plout(1) = a
        plout(2) = b
      end if
c
      return
      end
c
      subroutine mati3d10(d,md,nh)
c-----------------------------------------------------------------------
c
c     mate3d10 concrete model
c     input material parameter
c
c-----------------------------------------------------------------------
c
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(9)

      md=9
      if(ior.lt.0) write(*,1001)
1001  format(
     + ' Input:E,v,ipla,fctm,fcm,g_f,g_f,gam1,gam2')
      nh = 8
      call dinput(d,md)
                  write(iow,1002) nh,(d(i),i=1,md)
      if(ior.lt.0)write(*  ,1002) nh,(d(i),i=1,md)
1002  format(5x,'concrete material data',/,
     +  5x,'length nh of h1,h2 ..................',i12,/,
     +  5x,'E     Youngs modulus ................',g12.5,/,
     +  5x,'v     Poissons ratio ................',g12.5,/,
     +  5x,'ipla  plasticity on=1 /off=0 ........',g12.5,/,
     +  5x,'f_ctm tensile strength ..............',g12.5,/,
     +  5x,'f_cm  compressive strength ..........',g12.5,/,
     +  5x,'G_f   energy release during cracking.',g12.5,/,
     +  5x,'G_c   energy release rate ...........',g12.5,/,
     +  5x,'gam_1 fitting parameter .............',g12.5,/,
     +  5x,'gam_2 fitting parameter .............',g12.5,/)
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine matm3d10(d,sig,cmat,e,ep,a1,a2,lgp,dvp)
c-----------------------------------------------------------------------
c     Linear elastic isotropy
c
c...  storage of epsilon and sigma
c     11  22  33  12  13  23
c
c...  storage of elasticity-matrix
c     1111  1122  1133  1112  1131  1123
c     2211  2222  2233  2212  2231  2223
c     3311  3322  3333  3312  3331  3323
c     1211  1222  1233  1212  1231  1223
c     3111  3122  3133  3112  3131  3123
c     2311  2322  2333  2312  2331  2323
c
c...  factor 2.0 only on E_12, E_13, E_23, not on sig, sd, cmat, IP
c     !! handle with care !!
c
c-----------------------------------------------------------------------
      USE eldata !!! necessary WW??
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*),sig(*),cmat(6,6),e(6),ep(6)
      dimension ee(6),sd(6)
c
c
c
c.....independent inputvalues
      xchl  = dvp**(1.d0/3.d0)
      xE    = d(1)
      ipla  = d(3)
      fctm  = d(4)
      fcm   = d(5)
      g_f   = d(6)/xchl
      g_c   = d(7)/xchl
      fit1  = d(8)
      fit2  = d(9)
c
      pi   = 4.d0*datan(1.d0)
c
c.....dependent inputvalues
      alf1 = dsqrt(2.d0/3.d0)*(fit1*fcm-fctm)/(fit1*fcm+fctm)
      alf2 = dsqrt(2.d0/3.d0)*(fit2-1.d0)/(2.d0*fit2-1.d0)
      alf3 = -1.d0/3.d0/alf1
      alf5 = -1.d0/3.d0/alf2
      bet1 = (2.d0*fit1*fcm)/(fit1*fcm+fctm)
      bet2 = fit2/(2.d0*fit2-1.d0)
c
c.... trial strains  E^el = E - E^pl
      do i=1,6
         ee(i)=0.d0
         ee(i)=e(i)-ep(i)
      end do
c
c.....elasticity matrix
      call mat10(cmat,d)
c
c.....trial stresses   sig(i)=cmat(i,j)*ee(j)
      call mvmul(cmat,ee,6,6,sig)
      if (ipla.eq.0) goto 100
c
c.....deviatoric stresses  sds = |(sd^tr)^2|
c     factor 2.0 on sd_ij * sd_ij, i ne j
      sd(1)= 2.0d0/3.0d0*sig(1)-1.0d0/3.0d0*sig(2)-1.0d0/3.0d0*sig(3)
      sd(2)=-1.0d0/3.0d0*sig(1)+2.0d0/3.0d0*sig(2)-1.0d0/3.0d0*sig(3)
      sd(3)=-1.0d0/3.0d0*sig(1)-1.0d0/3.0d0*sig(2)+2.0d0/3.0d0*sig(3)
      sd(4)= sig(4)
      sd(5)= sig(5)
      sd(6)= sig(6)
      sds  = dsqrt(        (sd(1)*sd(1)+sd(2)*sd(2)+sd(3)*sd(3))
     +              + 2.d0*(sd(4)*sd(4)+sd(5)*sd(5)+sd(6)*sd(6)))
      sdsd = sds*sds
c
c.....for num. stability in nearly hydrostatic case
      if (sds.lt.1.d-16) sds = 1.d-16
c
c.....trace of stresses  tr(sig)
      trsig = sig(1) + sig(2) + sig(3)
c
c.....f1 - tension failure (exponential softening)
      a_tu = g_f/fctm
cjs      q1 = dsqrt(2.d0/3.d0)*bet1*fctm*dexp(-a1/a_tu)
      if (-a1/a_tu.lt.-5.2983d0) then
        fact1 = 5.d-3
      else
        fact1 = dexp(-a1/a_tu)
      endif
      q1 = dsqrt(2.d0/3.d0)*bet1*fctm*fact1
      f1 = sds + alf1*trsig - q1
c
c.....f2 - compression failure (square hardening / exp. softening)
      a_ce = 0.0022d0-fcm/xE     !Winkler 01  S.68
c      a_ce = 4.d0*fcm/3.d0/xE
      a_cu = 2.d0/dsqrt(pi) * g_c/fcm
      if (a2.lt.a_ce) then
        y2 = bet2*fcm/3.d0*(1.d0 + 4.d0*a2/a_ce - 2.d0*a2*a2/a_ce/a_ce)
      else
cjs        y2 = bet2*fcm*dexp(-(a2-a_ce)*(a2-a_ce)/a_cu/a_cu)
        if (-(a2-a_ce)*(a2-a_ce)/a_cu/a_cu.lt.-5.2983d0) then
          fact2 = 5.d-3
c          if (lgp.eq.1) then
c            write ( * ,*) 'limit'
c            write (iow,*) 'limit'
c          end if
        else
          fact2 = dexp(-(a2-a_ce)*(a2-a_ce)/a_cu/a_cu)
        endif
        y2 = bet2*fcm*fact2
      endif
      q2 = dsqrt(2.d0/3.d0)*y2
      f2 = sds + alf2*trsig - q2
c
c.....f3 - inverted cone (num. stab. in nearly hydrost. tension case)
      q3 = -1.d0/(3.d0*alf1*alf1)*q1
      f3 = sds + alf3*trsig - q3
c
c.....f4 - yield sphere (square hardening / exp. softening)
      xL =      - (dsqrt(54.d0)*alf2 + 2.d0) * fit2 * y2/bet2
      xR = dsqrt(2.d0/3.d0 + 6.d0*alf2*alf2) * fit2 * y2/bet2
      f4 = dsqrt(sdsd + (trsig-xL)*(trsig-xL)/9.d0) -xR
c
c.....f5 - cone-sphere border
      q5 = alf5*xL
      f5 = sds + alf5*trsig - q5
c
c.....elastic - plastic switch
      if (a2.lt.a_ce) then       !hardening
        if (f1.le.1.0d-12.and.f2.le.1.0d-12.and.f3.gt.1.0d-12) then
          if (f5.le.1.0d-12) return
          if (f4.le.1.0d-10) return
          call sphere10(d,sig,sd,f4,sds,a2,ep,lgp,cmat,xchl)
          return
        else
          call cone10(d,sig,sd,f1,f2,f3,sds,trsig,a1,a2,ep,
     +                lgp,cmat,xchl)
          sd(1)= 2.d0/3.d0*sig(1)-1.d0/3.d0*sig(2)-1.d0/3.d0*sig(3)
          sd(2)=-1.d0/3.d0*sig(1)+2.d0/3.d0*sig(2)-1.d0/3.d0*sig(3)
          sd(3)=-1.d0/3.d0*sig(1)-1.d0/3.d0*sig(2)+2.d0/3.d0*sig(3)
          sd(4)= sig(4)
          sd(5)= sig(5)
          sd(6)= sig(6)
          sds  = dsqrt(        (sd(1)*sd(1)+sd(2)*sd(2)+sd(3)*sd(3))
     +                  + 2.d0*(sd(4)*sd(4)+sd(5)*sd(5)+sd(6)*sd(6)))
          y2 = bet2*fcm/3.d0*(1.d0+4.d0*a2/a_ce-2.d0*a2*a2/a_ce/a_ce)
          xL = -(dsqrt(54.d0)*alf2 + 2.d0) * fit2 * y2/bet2
          q5 = alf5*xL
          f5 = sds + alf5*trsig - q5
          if (f5.le.1.0d-12) return
          call sphere10(d,sig,sd,f4,sds,a2,ep,lgp,cmat,xchl)
        end if
      else       ! softening
        if (f5.le.1.0d-12) then
          if (f1.le.1.0d-12.and.f2.le.1.0d-12.and.f3.gt.1.0d-12) return
          call cone10(d,sig,sd,f1,f2,f3,sds,trsig,a1,a2,ep,
     +                lgp,cmat,xchl)
          sd(1)= 2.d0/3.d0*sig(1)-1.d0/3.d0*sig(2)-1.d0/3.d0*sig(3)
          sd(2)=-1.d0/3.d0*sig(1)+2.d0/3.d0*sig(2)-1.d0/3.d0*sig(3)
          sd(3)=-1.d0/3.d0*sig(1)-1.d0/3.d0*sig(2)+2.d0/3.d0*sig(3)
          sd(4)= sig(4)
          sd(5)= sig(5)
          sd(6)= sig(6)
          sds  = dsqrt(        (sd(1)*sd(1)+sd(2)*sd(2)+sd(3)*sd(3))
     +                  + 2.d0*(sd(4)*sd(4)+sd(5)*sd(5)+sd(6)*sd(6)))
cjs          y2 = bet2*fcm*dexp(-(a2-a_ce)*(a2-a_ce)/a_cu/a_cu)
          if (-(a2-a_ce)*(a2-a_ce)/a_cu/a_cu.lt.-5.2983d0) then
            fact2 = 5.d-3
          else
            fact2 = dexp(-(a2-a_ce)*(a2-a_ce)/a_cu/a_cu)
          end if
          y2 = bet2*fcm*fact2
          xL = -(dsqrt(54.d0)*alf2 + 2.d0) * fit2 * y2/bet2
          q5 = alf5*xL
          f5 = sds + alf5*trsig - q5
          if (f5.le.1.0d-12) return
          call sphere10(d,sig,sd,f4,sds,a2,ep,lgp,cmat,xchl)
        else
          if (f4.le.1.0d-10) return
          call sphere10(d,sig,sd,f4,sds,a2,ep,lgp,cmat,xchl)
          return
        end if
      end if
c
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine cone10(d,sig,sd,f1,f2,f3,sds,trsig,a1,a2,ep,
     +                  lgp,cmat,xchl)
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*),sig(*),sd(*),ep(*),cmat(6,6)
      dimension xn(6),xgna1(6),xna1(6),xgna2(6),xna2(6),xgna3(6),xna3(6)
      dimension xM(3,3)
      dimension cmat1(6,6),cmat2(6,6),cmat3(6,6),cmat4(6,6),cmat5(6,6)
      dimension pmat(6,6)
c.....independent inputvalues
      xE    = d(1)
      xk    = d(1)/(3.d0*(1.d0-2.d0*d(2)))
      xg    = d(1)/(2.d0*(1.d0+d(2)))
      fctm  = d(4)
      fcm   = d(5)
      g_f   = d(6)/xchl
      g_c   = d(7)/xchl
      fit1  = d(8)
      fit2  = d(9)
c
      pi   = 4.d0*datan(1.d0)
c
c      if (lgp.eq.1) then
c        write(*,*) 'onto cone'
c      endif
c
c.....dependent inputvalues
      a_tu = g_f/fctm
      a_ce = 0.0022d0-fcm/xE     !Winkler 01  S.68
c      a_ce = 4.d0*fcm/3.d0/xE
      a_cu = 2.d0/dsqrt(pi) * g_c/fcm
c
      alf1 = dsqrt(2.d0/3.d0)*(fit1*fcm-fctm)/(fit1*fcm+fctm)
      alf2 = dsqrt(2.d0/3.d0)*(fit2-1.d0)/(2.d0*fit2-1.d0)
      alf3 = -1.d0/3.d0/alf1
      bet1 = (2.d0*fit1*fcm)/(fit1*fcm+fctm)
      bet2 = fit2/(2.d0*fit2-1.d0)
c
c.....Startvalues
      c1 = 1.d0
      c2 = 1.d0
      c3 = 1.d0
      if (f1.le.1.0d-12)  c1 = 0.d0
      if (f2.le.1.0d-12)  c2 = 0.d0
      if (f3.gt.1.0d-12)  c3 = 0.d0      !opposite direction for f3
      gam1  = 0.d0
      gam2  = 0.d0
      gam3  = 0.d0
c
      a1h = a1
      a2h = a2
c
      newton = 0
c
c.....BEGIN LOCAL NEWTON ITERATION...................................
c
100   continue
      newton = newton + 1
      if (newton.gt.99) then
         write( * ,*) 'no converg. in local newton'
         write(iow,*) 'no converg. in local newton'
         stop
      end if
c
c.....frequently used factors
      sum1 = c1*gam1 + c2*gam2 + c3*gam3
      sum2 = 2.d0*xg*sum1
      sum3 = 9.d0*xk*(c1*gam1*alf1 + c2*gam2*alf2 + c3*gam3*alf3)
c
c.....yield conditions
cjs      q1 = dsqrt(2.d0/3.d0)*bet1*fctm*dexp(-a1/a_tu)
c      if (-a1/a_tu.lt.-29.9849d0) then
c        fact1 = 9.5d-14
      if (-a1/a_tu.lt.-5.2983d0) then
        fact1 = 5.d-3
      else if (-a1/a_tu.ge.2.d0) then
        fact1 = dexp(2.d0)
      else
        fact1 = dexp(-a1/a_tu)
      end if
      q1 = dsqrt(2.d0/3.d0)*bet1*fctm*fact1
      f1 = sds - sum2 + alf1 * (trsig - sum3) - q1
c
      if (a2.lt.a_ce) then
        y2 = bet2*fcm/3.d0*(1.d0 + 4.d0*a2/a_ce - 2.d0*a2*a2/a_ce/a_ce)
      else
cjs        y2 = bet2*fcm*dexp(-(a2-a_ce)*(a2-a_ce)/a_cu/a_cu)
        if (-(a2-a_ce)*(a2-a_ce)/a_cu/a_cu.lt.-5.2983d0) then
          fact2 = 5.d-3
        else
          fact2 = dexp(-(a2-a_ce)*(a2-a_ce)/a_cu/a_cu)
        end if
        y2 = bet2*fcm*fact2
      end if
      q2 = dsqrt(2.d0/3.d0)*y2
      f2 = sds - sum2 + alf2 * (trsig - sum3) - q2
c
      q3 = -1.d0/(3.d0*alf1*alf1)
cjs     +                   *dsqrt(2.d0/3.d0)*bet1*fctm*dexp(-a1/a_tu)
     +                   *dsqrt(2.d0/3.d0)*bet1*fctm*fact1
      f3 = sds - sum2 + alf3 * (trsig - sum3) - q3
c
c.....exit iteration?
      if ((dabs(c1*f1).lt.1.d-12).and.(dabs(c2*f2).lt.1.d-12)
     +                           .and.(dabs(c3*f3).lt.1.d-12)) goto 105
c
c.....derivations
cjs   dq1dg1 = -dsqrt(2.d0/3.d0)*bet1*fctm*dexp(-a1/a_tu)/a_tu
      dq1dg1 = -dsqrt(2.d0/3.d0)*bet1*fctm*fact1/a_tu
      df1dg1 = (-2.d0*xg - alf1*alf1*9.d0*xk)*c1 - dq1dg1
      df1dg2 = (-2.d0*xg - alf1*alf2*9.d0*xk)*c2
      df1dg3 = (-2.d0*xg - alf1*alf3*9.d0*xk)*c3
c
      if (a2.lt.a_ce) then
        dq2dg2 = dsqrt(2.d0/3.d0)*bet2*fcm/3.d0
     +                         *(4.d0/a_ce - 4.d0*a2/a_ce/a_ce)
      else
        dq2dg2 = -dsqrt(2.d0/3.d0)*bet2*fcm*2.d0*(a2-a_ce)/a_cu/a_cu
cjs     +                          *dexp(-(a2-a_ce)*(a2-a_ce)/a_cu/a_cu)
     +                          *fact2
      end if
      df2dg1 = (-2.d0*xg - alf2*alf1*9.d0*xk)*c1
      df2dg2 = (-2.d0*xg - alf2*alf2*9.d0*xk)*c2 - dq2dg2
      df2dg3 = (-2.d0*xg - alf2*alf3*9.d0*xk)*c3
c
      dq3dg1 = -1.d0/(3.d0*alf1*alf1)*dq1dg1
      df3dg1 = (-2.d0*xg - alf2*alf1*9.d0*xk)*c1 - dq3dg1
      df3dg2 = (-2.d0*xg - alf2*alf2*9.d0*xk)*c2
      df3dg3 = (-2.d0*xg - alf2*alf3*9.d0*xk)*c3
c
c.....M-Matrix and Inverse
      xM(1,1) = c1 * df1dg1 + 1.d0 - c1
      xM(1,2) = c1 * df1dg2
      xM(1,3) = c1 * df1dg3
      xM(2,1) = c2 * df2dg1
      xM(2,2) = c2 * df2dg2 + 1.d0 - c2
      xM(2,3) = c2 * df2dg3
      xM(3,1) = c3 * df3dg1
      xM(3,2) = c3 * df3dg2
      xM(3,3) = c3 * df3dg3 + 1.d0 - c3
      call pivot(xM,3,3,xM)
c
c.....next step values / exit values
      f1bar = c1*f1 + (1.d0-c1)*gam1
      f2bar = c2*f2 + (1.d0-c2)*gam2
      f3bar = c3*f3 + (1.d0-c3)*gam3
      dgam1 = -xM(1,1)*f1bar - xM(1,2)*f2bar - xM(1,3)*f3bar
      dgam2 = -xM(2,1)*f1bar - xM(2,2)*f2bar - xM(2,3)*f3bar
      dgam3 = -xM(3,1)*f1bar - xM(3,2)*f2bar - xM(3,3)*f3bar
c
      gam1 = gam1 + dgam1
      gam2 = gam2 + dgam2
      gam3 = gam3 + dgam3
c
      a1   = a1 + dgam1
      a2   = a2 + dgam2
c
      if (a1.lt.0.d0)  a1 = 0.d0
      if (a2.lt.0.d0)  a2 = 0.d0
c
c.....active/inactive yield-surface?
      if (gam1.le.0.d0) c1 = 0.d0
      if (gam1.gt.0.d0) c1 = 1.d0
      if (gam2.le.0.d0) c2 = 0.d0
      if (gam2.gt.0.d0) c2 = 1.d0
      if (gam3.lt.0.d0) c3 = 1.d0   !opposite direction for f3
      if (gam3.ge.0.d0) c3 = 0.d0
c
      goto 100
c
c.....END LOCAL NEWTON ITERATION.....................................
c
105   continue
c
c.....update equivalent plastic strains
      a1 = c1*a1 + (1.d0-c1)*a1h
      a2 = c2*a2 + (1.d0-c2)*a2h
c
c.....update stresses
c     factor 2.0 !
      do i=1,6
        xn(i)  = sd(i)/sds
        if (i.le.3) xgna1(i) = 2.d0*xg*xn(i) + 3.d0*xk*alf1
        if (i.le.3) xgna2(i) = 2.d0*xg*xn(i) + 3.d0*xk*alf2
        if (i.le.3) xgna3(i) = 2.d0*xg*xn(i) + 3.d0*xk*alf3
        if (i.gt.3) xgna1(i) = 2.d0*xg*xn(i)
        if (i.gt.3) xgna2(i) = 2.d0*xg*xn(i)
        if (i.gt.3) xgna3(i) = 2.d0*xg*xn(i)
c        if (i.gt.3) xgna1(i) = 2.d0*xg*xn(i) * 2.d0
c        if (i.gt.3) xgna2(i) = 2.d0*xg*xn(i) * 2.d0
      end do
      do i=1,6
        sig(i) = sig(i)-c1*gam1*xgna1(i)-c2*gam2*xgna2(i)
     +                 -c3*gam3*xgna3(i)
      end do
c
c.....update plastic strains
      do i=1,6
c       factor 2.0 on ep_ij, i ne j
        if (i.le.3) xna1(i)= xn(i) + alf1
        if (i.gt.3) xna1(i)= xn(i)
        if (i.le.3) xna2(i)= xn(i) + alf2
        if (i.gt.3) xna2(i)= xn(i)
        if (i.le.3) xna3(i)= xn(i) + alf3
        if (i.gt.3) xna3(i)= xn(i)
        if (i.le.3) ep(i)= ep(i) + c1*gam1*xna1(i) + c2*gam2*xna2(i)
     +                           + c3*gam3*xna3(i)
        if (i.gt.3) ep(i)= ep(i) +(c1*gam1*xna1(i) + c2*gam2*xna2(i)
     +                           + c3*gam3*xna3(i))*2.d0
      end do
c
c.....elasto-plastic tangent operator
c     values for f1
cjs      dq1dg1 = -dsqrt(2.d0/3.d0)*bet1*fctm*dexp(-a1/a_tu) /a_tu
      dq1dg1 = -dsqrt(2.d0/3.d0)*bet1*fctm*fact1 /a_tu
      term1  = c1*(  c1*2.d0*xg + 9.d0*xk*c1*alf1*alf1
     +             + c2*2.d0*xg + 9.d0*xk*c2*alf1*alf2
     +             + c3*2.d0*xg + 9.d0*xk*c3*alf1*alf3
     +             + dq1dg1)
c
c     values for f2
      if (a2.lt.a_ce) then
        dq2dg2 = dsqrt(2.d0/3.d0)*bet2*fcm/3.d0
     +                         *(4.d0/a_ce - 4.d0*a2/a_ce/a_ce)
      else
        dq2dg2 = -dsqrt(2.d0/3.d0)*bet2*fcm*2.d0*(a2-a_ce)/a_cu/a_cu
cjs     +                          *dexp(-(a2-a_ce)*(a2-a_ce)/a_cu/a_cu)
     +                          *fact2
      end if
      term2  = c2*(  c1*2.d0*xg + 9.d0*xk*c1*alf2*alf1
     +             + c2*2.d0*xg + 9.d0*xk*c2*alf2*alf2
     +             + c3*2.d0*xg + 9.d0*xk*c3*alf2*alf3
     +             + dq2dg2)
c
c     values for f3
      term3  = c3*(  c1*2.d0*xg + 9.d0*xk*c1*alf3*alf1
     +             + c2*2.d0*xg + 9.d0*xk*c2*alf3*alf2
     +             + c3*2.d0*xg + 9.d0*xk*c3*alf3*alf3
     +             - 1.d0/(3.d0*alf1*alf1)*dq1dg1)
c
      call pzero(pmat,6*6)
      pmat(1,1) =  2.d0/3.d0
      pmat(1,2) = -1.d0/3.d0
      pmat(1,3) = -1.d0/3.d0
      pmat(2,1) = -1.d0/3.d0
      pmat(2,2) =  2.d0/3.d0
      pmat(2,3) = -1.d0/3.d0
      pmat(3,1) = -1.d0/3.d0
      pmat(3,2) = -1.d0/3.d0
      pmat(3,3) =  2.d0/3.d0
c      pmat(4,4) =  1.d0
c      pmat(5,5) =  1.d0
c      pmat(6,6) =  1.d0
      pmat(4,4) =  0.5d0
      pmat(5,5) =  0.5d0
      pmat(6,6) =  0.5d0
c
      fac1 = 0.d0
      fac2 = 0.d0
      fac3 = 0.d0
      fac4 = 0.d0
      fac5 = 0.d0
      if (c1.eq.1.d0) then
        fac1 = fac1 + 9.d0*xk*xk*alf1*alf1/term1
        fac2 = fac2 + 6.d0*xg*xk*alf1/term1
        fac3 = fac2
        fac4 = fac4 + 4.d0*xg*xg/term1
        fac5 = fac5 + 4.d0*gam1*xg*xg/sds
      end if
      if (c2.eq.1.d0) then
        fac1 = fac1 + 9.d0*xk*xk*alf2*alf2/term2
        fac2 = fac2 + 6.d0*xg*xk*alf2/term2
        fac3 = fac2
        fac4 = fac4 + 4.d0*xg*xg/term2
        fac5 = fac5 + 4.d0*gam2*xg*xg/sds
      end if
      if (c3.eq.1.d0) then
        fac1 = fac1 + 9.d0*xk*xk*alf3*alf3/term3
        fac2 = fac2 + 6.d0*xg*xk*alf3/term3
        fac3 = fac2
        fac4 = fac4 + 4.d0*xg*xg/term3
        fac5 = fac5 + 4.d0*gam3*xg*xg/sds
      end if
c
      call mat10(cmat,d)
      call pzero(cmat1,6*6)
      do i=1,3
        do j=1,3
          cmat1(i,j) = fac1
        end do
      end do
c
      call pzero(cmat2,6*6)
      call pzero(cmat3,6*6)
      do i=1,6
        cmat2(i,1) = xn(i)*fac2
        cmat2(i,2) = cmat2(i,1)
        cmat2(i,3) = cmat2(i,1)
c
        cmat3(1,i) = xn(i)*fac3
        cmat3(2,i) = cmat3(1,i)
        cmat3(3,i) = cmat3(1,i)
        do j=1,6
          cmat4(i,j) = xn(i)*xn(j)*fac4
          cmat5(i,j) = (pmat(i,j) - xn(i)*xn(j))*fac5
        end do
      end do
c
      do i=1,6
        do j=1,6
          cmat(i,j) = cmat(i,j) - cmat1(i,j) - cmat2(i,j) - cmat3(i,j)
     +                          - cmat4(i,j) - cmat5(i,j)
        end do
      end do
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sphere10(d,sig,sd,f,sds,a,ep,lgp,cmat,xchl)
c-----------------------------------------------------------------------
      USE eldata !!! necessary WW??
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*),sig(*),sd(*),ep(*),cmat(6,6)
      dimension xM(2,2),ten1(6),ten1T(1,6),dfds(6),dfdsT(1,6)
      dimension dfdss(6,6),dfdsy(6),dfdsyT(1,6),theta(6,6)
      dimension xN1(6),xN1T(1,6),xN2(6),xN2T(1,6)
      dimension xN11(6,6),xN12(6,6),xN21(6,6),xN22(6,6),Zmat(2,2)
c
c
c      if (lgp.eq.1) then
c        write(*,*) 'onto sphere'
c      endif
c
c.....independent Inputvalues
      Ecm   = d(1)
      xk    = d(1)/(3.d0*(1.d0-2.d0*d(2)))
      xg    = d(1)/(2.d0*(1.d0+d(2)))
      fcm   = d(5)
      Gc    = d(7)/xchl
      fit2  = d(9)
c
      pi    = 4.d0*datan(1.d0)
c
c.....dependent Inputvalues
      alfa  = dsqrt(2.d0/3.d0)*(fit2-1.d0)/(2.d0*fit2-1.d0)
      trsig = sig(1) + sig(2) + sig(3)
      sdsd  = sds*sds
c
c.....square hardening / softening
      a_ce = 0.0022d0-fcm/Ecm     !Winkler 01  S.68
c      a_ce = 4.d0*fcm/(3.d0*Ecm)
      a_cu = 2.d0/dsqrt(pi) * gc/fcm
      if (a.lt.a_ce) then
         y   = fcm/3.d0 * (1.d0 +4.d0*a/a_ce -2.d0*a*a/(a_ce*a_ce))
      else
cjs         y   = fcm*dexp(-(a-a_ce)*(a-a_ce)/a_cu/a_cu)
        if (-(a-a_ce)*(a-a_ce)/a_cu/a_cu.lt.-5.2983d0) then
          fact2 = 5.d-3
        else
          fact2 = dexp(-(a-a_ce)*(a-a_ce)/a_cu/a_cu)
        end if
        y = fcm*fact2
      end if
c
c.....origin, radius and drivates
      xL   = -(dsqrt(54.d0)*alfa + 2.d0) * fit2 * y
      dLdy = -(dsqrt(54.d0)*alfa + 2.d0) * fit2
      xR   = dsqrt(2.d0/3.d0 + 6.d0*alfa*alfa) * fit2 * y
      dRdy = dsqrt(2.d0/3.d0 + 6.d0*alfa*alfa) * fit2
c
c.....startvalues
      gam    = 0.d0
      newton = 0
c
c.....frequently used factors
      trSL  = trsig-xL
      fac1  = xR+gam*2.d0*xg
      fac2  = xR+gam*xk
      fac3  = xR*dRdy*2.d0*xg
      fac4  = xR*dRdy*xk
      fac5  = (xR/fac1)**2 *sdsd +(xR/fac2)**2 /9.d0 *trSL*trSL
      fac6  = fac3*gam/fac1**3 *sdsd + fac4*gam/fac2**3 *trSL*trSL /9.d0
     +        -(xR/fac2)**2 *trSL*dLdy /9.d0
c
c
c.....BEGIN LOCAL NEWTON ITERATION...................................
c
100   continue
      newton = newton + 1
      if (newton.gt.99) then
         write( * ,*) 'no converg. in local newton [sphere02]'
         write(iow,*) 'no converg. in local newton [sphere02]'
         stop
      end if
c
c.....yield function and derivates
      f     = dsqrt(fac5) -xR
c
      dfdy  = fac6/dsqrt(fac5) - dRdy
      dfdg  = -xR*xR*(2.d0*xg/fac1**3 *sdsd
     +             + xk/fac2**3 * trSL*trSL/9.d0)/dsqrt(fac5)
c
      dfdyy = ( dRdy*dRdy*gam*2.d0*xg*(fac1-3.d0*xR)/fac1**4 *sdsd
     +         +dRdy*dRdy*gam*xk*(fac2-3.d0*xR)
     +                                  /fac2**4 *trSL*trSL/9.d0
     +         -fac4*gam/fac2**3 *4.d0*trSL*dLdy /9.d0
     +         +(xR*dLdy/fac2)**2 /9.d0 )  /dsqrt(fac5)
     +        -fac6*fac6 /(fac5**1.5)
      dfdyg = ( fac3*(fac1-6.d0*gam*xg)/fac1**4 *sdsd
     +         +fac4*(fac2-3.d0*gam*xk)/fac2**4 *trSL*trSL/9.d0
     +         +2.d0*xR*xR*xk/fac2**3 * trSL*dLdy/9.d0 ) /dsqrt(fac5)
     +        +fac6*( xR*xR*2.d0*gam/fac1**3 *sdsd
     +               +xR*xR*xk/fac2**3 *trSL*trSL /9.d0) /fac5**1.5
c
      if ((a - gam*dfdy).lt.a_ce) then      !hardening
        fac  = (a - gam*dfdy)/a_ce
        g    = y - fcm/3.d0 * (1.d0 +4.d0*fac -2.d0*fac*fac)
        dgdy = 1.d0 + 4.d0/3.d0 *fcm*gam*dfdyy*(1.d0 - fac) / a_ce
        dgdg = 4.d0/3.d0 *fcm*(dfdy + gam*dfdyg)*(1.d0 - fac) / a_ce
      else                                  !softening
        fac  = (a - gam*dfdy - a_ce)/a_cu
        if (dabs(fac).ge.3.39307d0) fac = 3.39307d0
        g    = y - fcm*dexp(-fac*fac)
        dgdy = 1.d0 - 2.d0*fac*fcm*gam*dfdyy*dexp(-fac*fac) / a_cu
        dgdg = 2.d0*fac*fcm*(-dfdy-gam*dfdyg)*dexp(-fac*fac) / a_cu
      end if
c
ctest
c    if (lgp.eq.1) then
c      write( * ,*) 'f4=   ',f,  ' g4=   ',g
c    endif
ctest
c.....exit iteration?
      if ((dabs(f).lt.1.d-10).and.(dabs(g).lt.1.d-10)) goto 105
c      if (dabs(f).lt.1.d-12) f = 1.d-13
c
c.....M-Matrix and Inverse
      xM(1,1) = dfdg
      xM(1,2) = dfdy
      xM(2,1) = dgdg
      xM(2,2) = dgdy
      call invert(xM,2,2)
c
c.....next step values / exit values
      dgam = -xM(1,1)*f -xM(1,2)*g
      dy   = -xM(2,1)*f -xM(2,2)*g
c
      gam = gam + dgam
      y   = y   + dy
      if (dabs(y).lt.1.d-12) y = 1.d-12
c
c.....update origin, radius and drivates
      xL   = -(dsqrt(54.d0)*alfa + 2.d0) * fit2 * y
      dLdy = -(dsqrt(54.d0)*alfa + 2.d0) * fit2
      xR   = dsqrt(2.d0/3.d0 + 6.d0*alfa*alfa) * fit2 * y
      dRdy = dsqrt(2.d0/3.d0 + 6.d0*alfa*alfa) * fit2
c
c.....update factors
      trSL  = trsig-xL
      fac1  = xR+gam*2.d0*xg
      fac2  = xR+gam*xk
      fac3  = xR*dRdy*2.d0*xg
      fac4  = xR*dRdy*xk
      fac5  = (xR/fac1)**2 *sdsd +(xR/fac2)**2 /9.d0 *trSL*trSL
      fac6  = fac3*gam/fac1**3 *sdsd + fac4*gam/fac2**3 *trSL*trSL /9.d0
     +        -(xR/fac2)**2 *trSL*dLdy /9.d0
c
      goto 100
c
c.....END LOCAL NEWTON ITERATION.....................................
c
105   continue
c
c.....update equivalent plastic strains
      a     = a - dfdy*gam
c
c.....update stresses
      trsig = xR/fac2 * (trsig + gam*xk*xL/xR)
      do i=1,6
        sd(i) = sd(i)*xR/fac1
        if (i.le.3) sig(i) = sd(i) + trsig/3.d0
        if (i.gt.3) sig(i) = sd(i)
      end do
      sds  = dsqrt(        (sd(1)*sd(1)+sd(2)*sd(2)+sd(3)*sd(3))
     +              + 2.d0*(sd(4)*sd(4)+sd(5)*sd(5)+sd(6)*sd(6)))
      sdsd = sds*sds
c
c.....frequently used factors and tensors
      trSL = trsig-xL
      faca = sdsd + trSL*trSL/9.d0
      facb = dsqrt(faca)
      facc = trSL*dLdy /9.d0 /faca**1.5
      do i=1,6
        if (i.le.3) ten1(i) = sd(i) + trSL/9.d0
        if (i.gt.3) ten1(i) = sd(i)
        ten1T(1,i) = ten1(i)
      end do
c
c.....gradients for the elasto-plastic tangent operator
      do i=1,6
        dfds(i)    = ten1(i) / facb
        dfdsT(1,i) = dfds(i)
      end do
c
      call matmulf(ten1,ten1T,6,1,6,dfdss)
      do i=1,6
        do j=1,6
          dfdss(i,j) = - dfdss(i,j) / faca**1.5
        end do
      end do
      dfdss(1,1) = (  2.d0/3.d0 + 1.d0/9.d0 ) / facb + dfdss(1,1)
      dfdss(1,2) = ( -1.d0/3.d0 + 1.d0/9.d0 ) / facb + dfdss(1,2)
      dfdss(1,3) = ( -1.d0/3.d0 + 1.d0/9.d0 ) / facb + dfdss(1,3)
      dfdss(2,1) = ( -1.d0/3.d0 + 1.d0/9.d0 ) / facb + dfdss(2,1)
      dfdss(2,2) = (  2.d0/3.d0 + 1.d0/9.d0 ) / facb + dfdss(2,2)
      dfdss(2,3) = ( -1.d0/3.d0 + 1.d0/9.d0 ) / facb + dfdss(2,3)
      dfdss(3,1) = ( -1.d0/3.d0 + 1.d0/9.d0 ) / facb + dfdss(3,1)
      dfdss(3,2) = ( -1.d0/3.d0 + 1.d0/9.d0 ) / facb + dfdss(3,2)
      dfdss(3,3) = (  2.d0/3.d0 + 1.d0/9.d0 ) / facb + dfdss(3,3)
      dfdss(4,4) = (  1.d0                  ) / facb + dfdss(4,4)
      dfdss(5,5) = (  1.d0                  ) / facb + dfdss(5,5)
      dfdss(6,6) = (  1.d0                  ) / facb + dfdss(6,6)
c
      dfdy  = -trSL*dLdy /9.d0 /facb - dRdy
c
      dfdyy = dLdy*dLdy /9.d0 /facb - (trSL*(-dLdy)/9.d0)**2 / faca**1.5
c
      do i=1,6
        if (i.le.3) dfdsy(i) = ten1(i)*facc -dLdy/9.d0/facb
        if (i.gt.3) dfdsy(i) = ten1(i)*facc
        dfdsyT(1,i) = dfdsy(i)
      end do
c
c.....hardening tangent modul H
      if (a.lt.a_ce) then
        dyda  = fcm/3.d0 * (4.d0/a_ce -4.d0*a/a_ce/a_ce )
      else
cjs        dyda = -2.d0*fcm*(a-a_ce)/a_cu/a_cu
cjs     +              *dexp(-(a-a_ce)*(a-a_ce)/a_cu/a_cu)
        if (-(a-a_ce)*(a-a_ce)/a_cu/a_cu.lt.-5.2983d0) then
          fact2 = 5.d-3
        else
          fact2 = dexp(-(a-a_ce)*(a-a_ce)/a_cu/a_cu)
        end if
        dyda = -2.d0*fcm*(a-a_ce)/a_cu/a_cu * fact2
      end if
c      xH = -dyda
      xH = dyda
c
c.....further elements for the elasto-plastic tangent operator
      call mat10(cmat,d)
      call invert(cmat,6,6)     ! NB: cmat -> cmat^-1
      do i=1,6
        do j=1,6
          theta(i,j) = cmat(i,j) + gam * dfdss(i,j)
        end do
      end do
      call invert(theta,6,6)
c
      call mvmul(theta, dfds,6,6,xN1)
      call mvmul(theta,dfdsy,6,6,xN2)
      do i=1,6
         xN1T(1,i) = xN1(i)
         xN2T(1,i) = xN2(i)
      end do
      call matmulf(xN1,xN1T,6,1,6,xN11)
      call matmulf(xN1,xN2T,6,1,6,xN12)
      call matmulf(xN2,xN1T,6,1,6,xN21)
      call matmulf(xN2,xN2T,6,1,6,xN22)
c
      call mvmul(dfdsT ,xN1,1,6,z11)
      call mvmul(dfdsT ,xN2,1,6,z12)
      call mvmul(dfdsyT,xN2,1,6,z22)
      Zmat(1,1) = z11
      Zmat(1,2) = z12 -  dfdy/gam
      Zmat(2,1) = z12 -  dfdy/gam
      Zmat(2,2) = z22 - dfdyy/gam - 1.d0 /xH /gam/gam
      call invert(Zmat,2,2)
c
c.....elasto-plastic tangent operator
      call pzero(cmat,6*6)
      do i=1,6
        do j=1,6
          cmat(i,j) = theta(i,j) - Zmat(1,1) * xN11(i,j)
     +                           - Zmat(1,2) * xN12(i,j)
     +                           - Zmat(2,1) * xN21(i,j)
     +                           - Zmat(2,2) * xN22(i,j)
        end do
      end do
c
c.....plastic strains
      do i=1,6
c        ep(i)= ep(i) + gam * dfds(i)
        if (i.le.3) ep(i)= ep(i) + gam * dfds(i)
        if (i.gt.3) ep(i)= ep(i) + gam * dfds(i) *2.d0
      end do
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine mat10(cmat,d)
c-----------------------------------------------------------------------
c.....elasticity matrix plain stress element
c.....d(1)=E  d(2)=xnu
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension cmat(6,6),d(*)
      e1  = d(1)
      xnu = d(2)
c
      call pzero(cmat,6*6)
c
c.... local iostropy elasticity matrix
      cc = e1/((1.0d0+xnu)*(1.0d0-2.0d0*xnu))
      tt = (1.0d0-2.0d0*xnu)/2.0d0
c
c.... material parameters
      cmat(1,1) = cc*(1.0d0-xnu)
      cmat(1,2) = cc*xnu
      cmat(1,3) = cc*xnu
      cmat(2,2) = cc*(1.0d0-xnu)
      cmat(2,3) = cc*xnu
      cmat(3,3) = cc*(1.0d0-xnu)
      cmat(4,4) = cc*tt
      cmat(5,5) = cc*tt
      cmat(6,6) = cc*tt
c
c.... symmetry
      cmat(2,1) = cmat(1,2)
      cmat(3,1) = cmat(1,3)
      cmat(3,2) = cmat(2,3)
c
      return
      end
c