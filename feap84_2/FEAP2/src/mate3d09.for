      subroutine mate3d09(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c     Purpose: calculate S and C for
c              a small strain isotropic damage material law
c
c     Inputs:
c         h1(nh)      - history array h1
c         h1(1)       - r_n - Damage Threshold or tau_eps_n in case of Impl-Ex scheme
c         h1(2)       - d_n - Damage variable
c         h1(3)       - dr in case of Impl-Ex scheme or tau_eps_n in case of viscodamage
c         d(md)       - local d-array
c         Eps         - strains
c         isw         - solution option from element
c
c     Input material parameters:
c         E           - Young's modulus
c         nu          - Poisson's ratio
c         y_0         - ultimate stress
c         G_f         - Fracture energy
c         implex      - 0 = exact integration, != 0 implex (only without viscosity)
c        (eta         - Viscosity parameter         only for matm3d09_2) 
c        (alpha       - 1 = Backward Euler Method   only for matm3d09_2)
c
c     Outputs:
c         md = 8      - number of used data for control of d-array
c         nh = 4      - number of history parameter at Gauss-Point
c         h2(nh)      - history array h2
c         sig         - stresses
c         Cmat        - tangential stiffness Matrix
c         plout(10)   - plot data
c
c     Allocation of d-array:
c       for  matm3d09_1 & matm3d09_2
c         d(1): e1                 : Young's modulus
c         d(2): nu                 : Poisson's ratio
c         d(3): r_0 = y_0/sqrt(E)  : initial damage threshold parameter
c         d(4): G_f                : Fracture energy
c         d(5): y_0                : ultimate stress
c         d(6): 1 = implex, 0 = implicit                  
c       for  matm3d09_3
c         d(7): eta  : Viscosity Parameter
c         d(8): alpha: 1 = Backward Euler Method
c
c     Impl-Ex scheme based on:
c     An implicit/explicit integration scheme to increase cimputability of
c     non-linear material and contact/friction problems, J. Oliver, A.E. Huespe, J.C. Cante
c     Computer methods in applied mechanics and engineering, 197 (2008) p. 1865-1889
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),Eps(*),Sig(*),Cmat(6,6),plout(10)
c
      if(isw.eq.1) then
c....   input data
        call mati3d09(d,md,nh)
      else if (d(7).eq.0.and.d(6).eq.0) then
c....   Isotropic Damage Model, exact integration
        call matm3d09_1(h1,h2,cmat,eps,sig,plout,d,dvp)
      elseif(d(6).ne.0) then
c....   Isotropic Damage Model Impl_Ex scheme
        call matm3d09_2(h1,h2,cmat,eps,sig,plout,d,dvp,isw)   
      else
c....   Isotropic Viscodamage Model
        call matm3d09_3(h1,h2,cmat,eps,sig,plout,d,dvp)
      end if

      return
      end 
c
      subroutine matm3d09_1(h1,h2,cmat,eps,sig,plout,d,dvp)
c-----------------------------------------------------------------------
c
c                 Isotropic Damage Model
c
c-----------------------------------------------------------------------

      USE eldata
      USE iofile
      implicit double precision(a-h,o-z)
      double precision :: lc,lambda
      dimension d(*),cmat(6,6),eps(6),sig_eff(6),sig(6),plout(*),
     +          sig_tensor(6,6), cmd(6,6), h1(*),h2(*)

      e1   = d(1)      ! Young's modulus
      xnu  = d(2)      ! Poisson's ratio
      r_0  = d(3)      ! initial damage threshold
      G_f  = d(4)      ! Fracture energy
      y_0  = d(5)      ! ultimate stress

c....  copy time n data from h1-array
      r_n  = h1(1)     ! Damage threshold
      d_n  = h1(2)     ! Damage value

c---------------------------------------------------------------------
c      Initial conditions
c---------------------------------------------------------------------
      call matm3d01(cmat,e1,xnu)

c....  effective stresses
      sig_eff = matmul(cmat,eps)

c....  effective stress Tensor (dyadic product of sig_eff)
       do i=1,6
        do j=1,6
         sig_tensor(i,j) = sig_eff(i)*sig_eff(j)
        enddo
       enddo

c....  free HELMHOLTZ-Energy psi_0
       psi_0 = 0.5d0*dot_product(sig_eff,eps)

c....  characteristic element lenght and Parameter H_e; dvp=V^e/n_int
       lc = dvp**(1.d0/3.d0)
       lambda = G_f*d(1)/(y_0*y_0)
       if (lc.ge.lambda) stop 'lc > lambda: mesh refinement necessary!'
       H_e = 1.d0/( (G_f*d(1)/(lc*y_0*y_0)) -0.5d0 )

c---------------------------------------------------------------------
c      Damage criterion
c---------------------------------------------------------------------
      r_n = max(r_0,r_n)               ! Damage threshold

      tau_eps_np1 = dsqrt(2.d0*psi_0)  ! Norm in strain space
c                                       (equivalent strain)
      F_eps = (tau_eps_np1 - r_n)      ! yield function in strain space
      TOL   = 1.d-8*dsqrt(e1)
      if (F_eps.le.TOL.and. F_eps.gt.-TOL) F_eps = TOL ! shift

      if(F_eps .le. 0.d0) then !....................elastic or unloading case
       d_np1 = d_n
       r_np1 = r_n
       sig = (1.d0-d_np1)*sig_eff
       cmd = (1.d0-d_np1)*cmat

      else                      !.............................damage case
       r_np1 = tau_eps_np1
       d_np1 = 1.d0 - (r_0/r_np1) * dexp(H_e * (1.d0 - (r_np1/r_0)))
       sig   = (1.d0-d_np1)*sig_eff
       dd_np1= ((r_0 + H_e * r_np1) / (r_np1*r_np1)
     +          *dexp(H_e*(1.d0-r_np1/r_0)))
       cmd = (1.d0-d_np1)*cmat - dd_np1*(1.d0/tau_eps_np1)*sig_tensor
      endif

c....  copy cmd-matrix into cmat-matrix
      cmat = cmd

c....  save time n+1 data from local array to h2-array
      h2(1) = r_np1
      h2(2) = d_np1
c....  store inelastic data for plotting
      plout(1) = r_np1
      plout(2) = d_np1
      plout(3) = dsqrt(eps(1)*eps(1)+eps(2)*eps(2)+eps(3)*eps(3)
     +               +(eps(4)*eps(4)+eps(5)*eps(5)+eps(6)*eps(6))*0.5d0)
      return
      end 
c
      subroutine matm3d09_3(h1,h2,cmat,eps,sig,plout,d,dvp)
c-----------------------------------------------------------------------
c
c         Isotropic Viscodamage Model
c
c-----------------------------------------------------------------------

      USE tdata
      implicit double precision(a-h,o-z)
      double precision :: lc,dtt,lambda
      dimension d(*),cmat(6,6),eps(6),sig_eff(6),sig(6),plout(*),
     +          sig_tensor(6,6), cmd(6,6), h1(*),h2(*)

      e1   = d(1)      ! Young's modulus
      xnu  = d(2)      ! Poisson's ratio
      r_0  = d(3)      ! initial damage threshold
      G_f  = d(4)      ! Fracture energy
      y_0  = d(5)      ! ultimate stress
      eta  = d(7)      ! Viscosity Parameter
      alpha= d(8)      ! 1 = Backward Euler

c....  copy time n data from h1-array 
      r_n  = h1(1)    ! Damage threshold
      d_n  = h1(2)    ! Damage value
      tau_eps_n = h1(3)

c---------------------------------------------------------------------
c     Initial conditions
c---------------------------------------------------------------------
      call matm3d01(cmat,e1,xnu)

c....  effective stresses
       sig_eff = matmul(cmat,eps)

c....  effective stress Tensor (dyadic product of sig_eff)
       do i=1,6
        do j=1,6
         sig_tensor(i,j) = sig_eff(i)*sig_eff(j)
        enddo
       enddo

c....  free HELMHOLTZ-Energy psi_0
       psi_0 = 0.5d0*dot_product(sig_eff,eps)

c....  characteristic element lenght and Parameter H_e; dvp=V^e/n_int
       lc = dvp**(1.d0/3.d0)
       lambda = G_f*d(1)/(y_0*y_0)
       if (lc.ge.lambda) stop 'lc > lambda: mesh refinement necessary!'
       H_e = 1.d0/( (G_f*d(1)/(lc*y_0*y_0)) -0.5d0 )

c---------------------------------------------------------------------
c     Damage criterion
c---------------------------------------------------------------------
      r_n = max(r_0,r_n)               ! Damage threshold

      tau_eps_np1 = dsqrt(2.d0*psi_0)  ! Norm in strain space
c                                       (equivalent strain)
      dtt=abs(dt)
      tau_eps_alpha = (1.d0-alpha)*tau_eps_n + alpha * tau_eps_np1
      r_np1_trial = (((1.d0-(dtt/eta) * (1.d0-alpha)))*r_n + (dtt/eta)*
     + tau_eps_alpha) / (1.d0 + alpha*(dtt/eta))
      r_np1_trial = max(r_0,r_np1_trial)
      r_np1 = (1.d0-alpha) * r_n + alpha * r_np1_trial !

      if(r_np1_trial .le. r_n) then !....................elastic
       d_np1 = d_n
       r_np1 = r_n
       sig = (1.d0-d_np1)*sig_eff
       cmd = (1.d0-d_np1)*cmat

      else
       r_np1 = r_np1_trial
       d_np1 = 1.d0 - (r_0/r_np1) * dexp(H_e * (1.d0 - (r_np1/r_0)))
       sig   = (1.d0-d_np1)*sig_eff
       dd_np1= ((r_0 + H_e * tau_eps_np1) /
     +          (tau_eps_np1*tau_eps_np1)
     +          *dexp(H_e*(1.d0-r_np1/r_0)))
       cmd = (1.d0-d_np1)*cmat-((alpha*dtt/eta)/(1.d0 + alpha*dtt/eta))*
     +       dd_np1*(1.d0/tau_eps_np1)*sig_tensor
      endif

c....  copy cm^d-matrix into cmat-matrix
      cmat = cmd

c....  save time n+1 data from local array to h2-array
      h2(1) = r_np1
      h2(2) = d_np1
      h2(3) = tau_eps_np1
c....  plot data
      plout(1) = r_np1
      plout(2) = d_np1
      plout(3) = dsqrt(eps(1)*eps(1)+eps(2)*eps(2)+eps(3)*eps(3)
     +               +(eps(4)*eps(4)+eps(5)*eps(5)+eps(6)*eps(6))*0.5d0)
      return
c
      end 
c-----------------------------------------------------------------------
c
c     Isotropic Damage Model Impl-Ex integration scheme
c
c-----------------------------------------------------------------------
      subroutine matm3d09_2(h1,h2,cmat,eps,sig,plout,d,dvp,isw)
      USE implstep
      USE tdata
      implicit double precision(a-h,o-z) 
      dimension d(*),cmat(6,6),eps(6),sig_eff(6),sig(6),plout(*)
      dimension h1(*),h2(*)
      
      
c...  input parameters      
      e1  = d(1)      ! Young's modulus
      xnu = d(2)      ! Poisson's ratio
      r_0 = d(3)      ! initial damage threshold
      G_f = d(4)      ! Fracture energy
      y_0 = d(5)      ! ultimate stress
      
c...  history variables
      r_n         = h1(1)     ! Damage threshold
      d_n         = h1(2)     ! Damage value
      dr          = h1(3)     ! Rate of damage threshold at n
      tau_eps_n   = h1(4) 
      dr  = dr*dt             ! explicit rate integration, increment for internal variable at current step
      
c---------------------------------------------------------------------
c Initial conditions
c---------------------------------------------------------------------
      call matm3d01(cmat,e1,xnu)
      
      sig_eff = matmul(cmat,eps)

c---------------------------------------------------------------------
c Damage criterion
c---------------------------------------------------------------------
      tau_eps_n   = max(r_0,tau_eps_n)
      tau_eps_np1 = dsqrt(dot_product(sig_eff,eps))     ! equivalent strain equals exact solution
                                                        ! for internal variable at current step
      
      xlc = dvp**(1.d0/3.d0)  ! characteristic length
      xla = G_f*e1/y_0/y_0
      if (xlc.ge.xla) stop 'lc > lambda: mesh refinement necessary!'
      
      if(dr.gt.0.d0) then          ! damage only in case of positiv increment
          H_e = 1.d0/( (G_f*e1/(xlc*y_0*y_0)) -0.5d0 )          
          r_np1 = tau_eps_n + dr
          d_np1 = 1.d0 - (r_0/r_np1) * dexp(H_e * (1.d0 - (r_np1/r_0)))
      else
          r_np1 = r_n
          d_np1 = d_n
      endif

      dr    = (tau_eps_np1 - tau_eps_n)/dt  ! Rate assumption based on exact solution for next step
      
      tau_eps_np1 = max(tau_eps_np1,tau_eps_n)
      dr    = (tau_eps_np1 - tau_eps_n)/dt 
      
c...  Update stresses and algrithmic       
      sig   = (1.d0-d_np1)*sig_eff
      cmat  = (1.d0-d_np1)*cmat
      
c...  store history variables
      h2(1) = r_np1
      h2(2) = d_np1
      h2(3) = dr
      h2(4) = tau_eps_np1
      
c...  Plot data      
      plout(1) = r_np1
      plout(2) = d_np1
      plout(3) = dsqrt(eps(1)*eps(1)+eps(2)*eps(2)+eps(3)*eps(3)
     +               +(eps(4)*eps(4)+eps(5)*eps(5)+eps(6)*eps(6))*0.5d0)
      
      
c--------------------------------------------------------------------------------
c     Adaptive time stepping
c--------------------------------------------------------------------------------
      
      if(isw.eq.19) then
          fac2 = dabs(dr - h1(3))
          fac2 = fac2/r_0
          
!$OMP ATOMIC
          tm = max(tm,fac2)
      endif
      
      return
      end
      
      
      subroutine mati3d09(d,md,nh)
c-----------------------------------------------------------------------
c     input material parameters for small strain isotropic damage model
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*),dd(7)
c
      md=8
      nh=4
      if(ior.lt.0) write(*,1001)
1001  format(' Input: E,v,y_0,G_f,implex(,eta,alpha) ')
      call dinput(dd,7)

      d(1) = dd(1)              ! Young's modulus
      d(2) = dd(2)              ! Poisson's ratio
      d(3) = dd(3)/dsqrt(dd(1)) ! initial damage threshold
      d(4) = dd(4)              ! Fracture energy
      d(5) = dd(3)              ! ultimate stress
      d(6) = dd(5)              ! 1 = Impl_Ex integration
      d(7) = dd(6)              ! Viscosity Parameter
      d(8) = dd(7)              ! 1 = Backward Euler
                  write(iow,1002) nh,(d(i),i=1,md)
      if(ior.lt.0)write(*  ,1002) nh,(d(i),i=1,md)
      if(d(7).ne.0.and.d(6).ne.0) then ! check parameters for integration model
          write(iow,1003)
          write(*,1003)
          stop
      endif
      
1002  format(5x,'Small strain isotropic damage model',/,
     + 5x,'length nh of h1,h2............',i12,/,
     + 5x,'Youngs modulus E..............',g12.4,/,
     + 5x,'Poissons ratio v..............',g12.4,/,
     + 5x,'Initial damage threshold r_0..',f12.4,/,
     + 5x,'Fracture Energy G_f...........',f12.4,/,
     + 5x,'Ultimate stress y_0 ..........',f12.4,/,
     + 5x,'1=Impl_ex, 0=std .............',f12.4,/,
     + 5x,'Viscosity Parameter eta.......',f12.4,/,
     + 5x,'Parameter alpha ..............',f12.4)
1003  format('Impl-Ex scheme with Viscosity not possible,
     + check input: d(5) and d(6)')
c
      return
      end

