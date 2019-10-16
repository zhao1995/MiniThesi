      subroutine mate3d04(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c     Purpose: calculate S and C for 
c              a small strain elasto-(visco-)plastic isotropic mat. law 
c              E = Ee + Ep     (E = Ee + Evp)
c
c     Inputs:
c         h1(nh)      - history array h1
c         d(md)       - local d-array
c         Eps         - strains  
c         isw         - solution option from element 
c
c     Input material parameter:                
c         E           - Young's modulus
c         v           - Poisson's ratio
c         y_o         - initial yield stress
c         y_i         - yield stress at t=infty
c         xk          - linear hardening modul
c         xd          - exponential hardening
c         eta         - viscosity parameter
c
c     Outputs:
c         md = 7      - number of used data for control of d-array 
c         nh = 7      - number of history parameter at Gauss-Point
c         h2(nh)      - history array h2
c         Sig         - 2nd Piola-Kirchhoff Stress
c         Cmat        - algorithimic consistent tangent modulus
c         plout(10)   - plot data    
c
c     Allocation of d-array:
c         d(1) = E    - Young's modulus
c         d(2) = v    - Poisson's ratio
c         d(3) = y_o  - initial yield stress
c         d(4) = y_i  - yield stress at t=infty
c         d(5) = xk   - linear hardening modul
c         d(6) = xd   - exponential hardening
c         d(7) = eta  - viscosity parameter
c         yield condition: f = g - y_o+ xk*a + (y_i-y_o)*(1-exp(-xd*a))
c                          g = sqrt(3/2 Sd:Sd) - sd deviatoric stress
c                          a = equivalent plastic strain
c     
c     History-terms 7
c         h(1) E_px
c         h(2) E_py
c         h(3) E_pz
c         h(4) E_pxy
c         h(5) E_pxz
c         h(6) E_pyz
c         h(7) a
c
c     References: 
c     [1] Mises, von R.; Mechanik der plastischen Formänderung von 
c         Metallen; ZAMM 8(3)(1928)161-185
c     [2] Schuett, J.;Theorie und Finite-Element-Implementierung eines
c         elasto-plastischen 3D Stoffgesetzes, Institut für Baustatik
c         Karlsruhe
c     [3] Perzyna, P., Fundamental Problems in Viscoplasticity,
c         Adv.Appl.Mech. 243-373 (1966)
c
c     (c) s.klinkel                                             May,2001
c-----------------------------------------------------------------------
      USE tdata
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),Eps(*),Sig(*),Cmat(6,6),E_p(6),
     +          plout(10)
c
      if(isw.eq.1) then
c....  input data
       call mati3d04(d,md,nh)
c
      else
       xm   = d(1)/(2.d0*(1.d0+d(2))) 
       xd   = d(6)
       eta  = d(7)
       if(dabs(dt).lt.1.d-20.and.eta.gt.1.d-20)stop 'set dt !!!'
       if(dabs(dt).gt.1.d-20)then
        fact = dabs(eta/dt/xm) 
        if(dabs(xd).gt.1.d-10.and.fact.gt.1.d-10)stop 
     +               'stop: visco-plasticity only with linear hardening'
       end if

c....  Get plastic variables from history variables at time n 
       do i = 1,6
          E_p(i) = h1(i)
       end do
       a = h1(7)
c
c      compute elastic-plastic tangent modulus and stresses
       call matm3d04(d,Eps,E_p,a,sig,cmat)
c
c....  store history variables at time n+1 in h2-array
       do i = 1,6
          h2(i) = E_p(i)
       end do
       h2(7) = a
c       
c....  setup plot data       
       plout(1)=a  
       plout(2)=dsqrt(eps(1)*eps(1)+eps(2)*eps(2)+eps(3)*eps(3)  
     +              +(eps(4)*eps(4)+eps(5)*eps(5)+eps(6)*eps(6))*0.5d0)   
      end if
c
      return
      end
c
      subroutine mati3d04(d,md,nh)
c-----------------------------------------------------------------------
c
c     mate3d04 small elastoplastic strains isotropic E = Ee + Ep
c     input material parameters 
c     
c-----------------------------------------------------------------------
c
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(7)

      md=7
      if(ior.lt.0) write(*,1001)
1001  format(
     + ' Input: E,nu,Y_0,Y_inf,K,d,eta')
      nh=7
      call dinput(d,md)
      if(d(2).le.-1.d0.or.d(2).ge.0.5d0)stop ' -1 < nue < 0.5  ! '
      
                  write(iow,1002) nh,(d(i),i=1,md)
      if(ior.lt.0)write(*  ,1002) nh,(d(i),i=1,md)
1002  format(5x,'Small strain elastoplastic material data',/,
     +  5x,'length nh of h1,h2 ...........',i12,/,
     +  5x,'Youngs modulus E .............',g12.5,/,
     +  5x,'Poissons ratio v .............',g12.4,/,
     +  5x,'Initial yield stress    y_0 ..',g12.5,/,
     +  5x,'Yield stress at t=infty y_i ..',g12.5,/,
     +  5x,'Linear hardening modul xk ....',g12.5,/,
     +  5x,'Exponential hardening  xd ....',g12.5,/,
     +  5x,'viscosity parameter eta ......',g12.5,/)
c
      return
      end 
c
      subroutine matm3d04(d,ec,ep,a,sig,cmat)
c-----------------------------------------------------------------------
c
c     mate3d04 small elastoplastic strains isotropic E = Ee + Ep
c     calculate C, Sigma_V 
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*),ep(6),sig(6),cmat(6,6),ee(6),sd(6),ec(6)
c    
      xE   = d(1)
      xn   = d(2)
      y_o  = d(3)
      y_i  = d(4)
      xk   = d(5)
      xd   = d(6)
c
c.....elasticity matrix
      call matm3d01(Cmat,xE,xn)
c      
c.... trial strains  E^el = E - E^pl
      ee = ec - ep
c
c.....trial stresses   sig = cmat*ee
      sig = matmul(Cmat,ee)
      if (dabs(y_o/xE).lt.1.e-15) return
c
c.....yield stress 
      y = y_o + xk*a + (y_i-y_o)*(1-dexp(-xd*a))
c
c.... compute g^trial
      trsig = (sig(1)+sig(2)+sig(3))/3.d0
      sd(1) = sig(1) - trsig 
      sd(2) = sig(2) - trsig
      sd(3) = sig(3) - trsig
      sd(4) = sig(4)
      sd(5) = sig(5)
      sd(6) = sig(6)
      g_tr = dsqrt( 3.d0/2.d0*(sd(1)*sd(1)+sd(2)*sd(2)+sd(3)*sd(3))
     +             +     3.d0*(sd(4)*sd(4)+sd(5)*sd(5)+sd(6)*sd(6)) ) 
c      
c.....Yield condition   f = g^trial - y
      f = g_tr - y
      if(f.gt.1.d-10*y)call plas3d04(d,sig,sd,trsig,f,g_tr,y,a,ep,cmat)
c
      return
c
      end
c
      subroutine plas3d04(d,sig,sd,trsig,f,g_tr,y,a,ep,cmat)
c-----------------------------------------------------------------------
      USE iofile
      USE tdata
      implicit double precision (a-h,o-z)
      dimension d(*),sig(*),sd(*),ep(*),cmat(6,6)
c
      xk = d(1)/(3.d0*(1.d0-2.d0*d(2)))       ! bulk modulus 
      xm = d(1)/(2.d0*(1.d0+d(2)))            ! shear modulus
      y_o  = d(3)
      y_i  = d(4)
      xh   = d(5)
      xd   = d(6)
      eta  = d(7)
      etaddt = 0.d0
      if(dabs(dt).gt.1.d-20)etaddt = eta/dt 
c.....start values
      gam = 0.d0
      an  = a
c
c.....begin local Newton iteration...................................
      nmax = 20
      if(dabs(xd).lt.1.d-12.or.dabs(eta).gt.1.d-12)nmax = 1 
      do newton = 1,nmax 
        if (newton.eq.20) then
           write( * ,*) 'more than 20 iterations in plas3d04'
           write(iow,*) 'more than 20 iterations in plas3d04'
           stop
        end if
        dgam = f /(3.d0*xm + xh + (y_i-y_o)*xd*exp(-xd*a) + etaddt )
        gam = gam + dgam
        a = an + gam
        y = y_o + xh*a + (y_i-y_o)*(1.d0 - exp(-xd*a))
        f = g_tr - 3.d0*xm*gam - y
c           write( * ,*) newton,f
c           write(iow,*) newton,f
        if (dabs(f).lt.1.d-10*y) go to 100 
      enddo
c.....end local Newton iteration.....................................
c
c.....elasto-plastic tangent operator
100   continue
      dyda = xh + (y_i-y_o)*xd*exp(-xd*a) + etaddt
      xnt =              sd(1)*sd(1)+sd(2)*sd(2)+sd(3)*sd(3)+
     +             2.d0*(sd(4)*sd(4)+sd(5)*sd(5)+sd(6)*sd(6))
      beta  = 1.d0 - 3.d0*xm*gam/g_tr
      beta1 = 2.d0*xm*beta
      beta2 = 2.d0*xm*(1.d0/(1.d0+dyda/(3.d0*xm)) + beta - 1.d0)/xnt
            
      cmat(1,1) = xk + beta1*2.d0/3.d0
      cmat(1,2) = xk - beta1*1.d0/3.d0
      cmat(1,3) = cmat(1,2)
      cmat(2,1) = cmat(1,2)
      cmat(2,2) = cmat(1,1)
      cmat(2,3) = cmat(1,2)
      cmat(3,1) = cmat(1,3)
      cmat(3,2) = cmat(2,3)
      cmat(3,3) = cmat(1,1)
      cmat(4,4) = beta1*0.5d0
      cmat(5,5) = cmat(4,4)
      cmat(6,6) = cmat(4,4)
      do i=1,6
         do j=1,i
            cmat(i,j) = cmat(i,j) - beta2*sd(i)*sd(j)
            cmat(j,i) = cmat(i,j)
         end do
      end do
c
c.....update stresses and plastic strains
      fac = 3.d0/2.d0*gam/g_tr
      do i=1,3
         sig(i) = beta * sd(i) + trsig
         ep(i) = ep(i) + fac*sd(i)
      end do
      fac = fac*2.d0
      do i=4,6
         sig(i) = beta * sd(i) 
         ep(i) = ep(i) + fac*sd(i)
      end do
c
      return
      end
