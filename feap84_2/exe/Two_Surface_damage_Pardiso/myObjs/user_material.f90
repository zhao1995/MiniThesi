!       module sharedConstants01
!         implicit none
!         save
! 
! !        parameter for comparison with Lukas Element
!         double precision, parameter:: cP1=0.8d4            !0.8d4
!         double precision, parameter:: linternal=1.2d-3      !1.2d-3
!         double precision, parameter:: eta=0.0d0
! 
!         double precision, parameter:: lambda=5826.923077d0      !lambda (Lame constant)
!         double precision, parameter:: mu=3884.615385d0         !mu (Lame constant)
!         integer, parameter::          ndam=1                 !ndam=1: strain equivalence, ndam=2: energy equivalence
!         double precision, parameter:: eta_cc=1.0d0                  !crack closure parameter: eta_cc=0: no crack closure, eta_cc=1: full crack closure
!         double precision, parameter:: eps_01=0.d0               !eps_0,xx (irreversible strain)
!         double precision, parameter:: eps_02=0.d0              !eps_0,yy (irreversible strain)
!         double precision, parameter:: eps_03=0.d0               !eps_0,zz (irreversible strain)
!         double precision, parameter:: eps_04=0.d0               !sqrt(2)*eps_0,xy (irreversible strain)
!         double precision, parameter:: eps_05=0.d0               !sqrt(2)*eps_0,yz (irreversible strain)
!         double precision, parameter:: eps_06=0.d0               !sqrt(2)*eps_0,xz (irreversible strain)
!         integer, parameter::          h_switch=4               !switch for hardening: 1=linear, 2=quadratic, 3= exponential, 4= -> Infty
!         double precision, parameter:: YC0=0.0354d0                !initial damage surface   
!         double precision, parameter:: k1=0.06d0                 !hardening variable (linear term)
!         double precision, parameter:: k2=10.d0                 !hardening variable (quadratic term)
!         double precision, parameter:: k0=10.d0                 !hardening variable (exponential)
!         double precision, parameter:: alpha=10.d0              !hardening variable (exponential)
!         double precision, parameter:: Hpen=1.d12                   !penalty parameter
!         double precision, parameter:: D0=0.995d0               !critical damage value
!         double precision, parameter:: kappa=0.d0                !parameter for gradient extension
!         integer, parameter::          loc_tang=1                !switch for local tangent, 1= analytical, 2= numerical
!         integer, parameter::          glo_tang=1                !switch for consistent tangent, 1= analytical
!         
! 
! 
!       end module sharedConstants01





subroutine SmallStrainUmat(d,eps,sig,CC)!SW: Here, you were supposed to implement a small strain umat based on the infinitesimal strain, i.e., with sigma and epsilon, not tau and FF.
!SW: arguments should be: Umat(eps, sigma, CC), where CC is the stiffness tensor
!SW: pseudo-code of this subroutine: 1. set up CC for infinitesimal theory for small strains. 2. compute sigma=matmul(CC,eps). That's it. NOTHING ELSE should happen in this Umat.
implicit none 
double precision, intent(in) :: d(*), eps(6)
double precision, intent(out) :: sig(6), CC(6,6)

double precision :: nu, lambda, konsE, mu
integer :: i,j,k
konsE=210.d3!d(1)!SW2: Okay, I never work with the 'd'-vector. If you want to do that, you certainly have to modify something in the element subroutine. But I don't know what. At the moment we should simply set E and nu to some values. We can later - when everything works - add the vector 'd' if you want.
nu=.3d0!d(2)
lambda= (nu*konsE)/((1.d0+nu)*(1.d0-(2.d0*nu)))
mu=konsE/(2.d0*(1.d0+nu))
 
    CC = 0.d0
    do j=1,6!SW2: changed 3 to 6: in normalized Voigt notation, there is 2*mu on the whole diagonal.
      CC(j,j) = 2.0d0*mu
    end do
!     CC(1:3,1:3) = CC(1:3,1:3) + lambda!SW2: can you do that? I would have thought one needs a do loop here. Let's first try with a do-loop and then later try your version, when we are sure that everything works.
    do i=1,3; do j=1,3; CC(i,j)=CC(i,j)+lambda; end do; end do
!     do k=4,6
!       CC(k,k) = mu
!     end do
   
sig= matmul(CC,eps)

end subroutine SmallStrainUmat


    
      
subroutine umatIsoGradDam(eps,DChi,h1,hn,sig,dSigDeps,&
    &dSigDDChi,pi,dPiDDChi,gradDChi,Xi,dXiDGradDChi,nHistVars,d)

implicit none

include  'tdata.h'
include  'swvars.h' 
! include 'elauto.h'

integer,intent(in)::nHistVars
double precision,intent(in):: eps(6), DChi, hn(nHistVars),gradDChi(3), d(*)
double precision,intent(out)::h1(nHistVars)
double precision,intent(inout)::dSigDeps(6,6), dSigDDChi(6),dPiDDChi,Xi(3)
double precision,intent(inout)::dXiDGradDChi(3,3), sig(6),pi

double precision::dDDChi, dDDeps(6)


double precision::E			!E=8000.d0
double precision::eta			!eta=8.00d0
double precision::nu			!nu=0.30d0
double precision::HChi			!HChi=3.0d3
double precision::yc			!yc=0.03d0
double precision::teta			!teta=0.12d0
double precision,parameter::tolD=0.0005d0
double precision::l			!l=1.3d-3


double precision::lambda, mu, SpurEps, epsquad, PP, PPr, sig0(6)
double precision::DD,Dn, piH, piHPr, epsn(6),phitr,psi0

integer::i,j,inelStep
!       write(*,*) 'sig=' hn(1:6)

! call schreibeb(eps,6,1,'eps ')

E	= d(11)
nu	= d(12)
yc	= d(15)
teta	= d(17)
HChi	= d(19)
eta	= d(21)
l	= d(23)

!write(*,*) E, nu, yc, teta, eta
!write(*,*) HChi

! call schreibeb(eps,6,1,'Eps1')

Dn=hn(10) 

!!!Xi
do i=1,3
  dXiDGradDChi(i,i)=dXiDGradDChi(i,i)+E*(l**2)
end do
Xi=Xi+E*(l**2)*gradDChi

do i=1,6; epsn(i)=hn(10+i); end do

lambda=E*nu/((1.d0-2.d0*nu)*(1.d0+nu))
mu=E/(2.d0*(1.d0+nu))
!write(*,*)'lambda',lambda
!write(*,*)'mu',mu


SpurEps=0.d0; epsquad=0.d0
do i=1,3; SpurEps=SpurEps+eps(i); end do
do i=1,6; epsquad=epsquad+eps(i)*eps(i); end do

psi0=lambda/2.d0*(SpurEps**2)+mu*epsquad

call getPAndPPr(Dn,DChi,HChi,PP,PPr) ! HChi*(D-DChi); HChi

phitr=psi0-PP-yc
sig0=2.d0*mu*eps
do i=1,3; sig0(i)=sig0(i)+lambda*SpurEps; end do

    
inelStep=1
if(phitr.le.0.d0) inelStep=0
if(duVnshSW.eq.1) then
  inelStep=int(hn(17)+1.d-6)!old step inelastic or not
!                 write(*,*), 'inelStep=', inelStep
endif

if(inelStep.eq.0)then
!if(phitr<=0.d0)then
!write(*,*)'1.Fall phi<0'
  DD=Dn
  dDDChi=0.d0
  do i=1,6; dDDeps(i)=0.d0; end do
  h1(17)=0.d0!elasticstep
else
  h1(17)=1.d0!inelastic step
! write(*,*)'2.Fall phi>0'
  DD=DChi+(psi0-yc)/HChi
  dDDChi=1.d0
  do i=1,6; dDDeps(i)=sig0(i)/HChi; end do
  if(DD>(1.d0-tolD))then
    DD=1.d0-tolD
    dDDChi=0.d0
    do i=1,6; dDDeps(i)=0.d0; end do;
  end if
end if




call getPiHAndPiHPr(DChi,teta,piH,piHPr)

call getPAndPPr(DD,DChi,HChi,PP,PPr)
pi=pi+piH-PP

sig=sig+(1.d0-DD)*sig0+eta/dt*(eps-epsn)
!write(*,*)'sig1',sig(1)
!write(*,*)'sig2',sig(2)
!write(*,*)'sig3',sig(3)

do i=1,3
  do j=1,3
    dSigDeps(i,j)=dSigDeps(i,j)+lambda*(1.d0-DD)
  end do
end do
do i=1,6
  dSigDeps(i,i)=dSigDeps(i,i)+2.d0*mu*(1.d0-DD)
end do
do i=1,6
  do j=1,6
    dSigDeps(i,j)=dSigDeps(i,j)-sig0(i)*dDDeps(j)
  end do
  dSigDeps(i,i)=dSigDeps(i,i)+eta/dt
end do

dSigDDChi=dSigDDChi-dDDChi*sig0

dPiDDChi=dPiDDChi+piHPr+PPr*(-dDDChi+1.d0)

if(piH/=piH) write(*,*) 'ERROR NaN'

h1(2:7)=sig(:)
h1(8)=piH !
h1(9)=PP
h1(10)=DD
h1(11:16)=eps(:)



!       write(*,*) 'dt=', dt 
!       call schreibeb(sig,6,1,'sig ')
!       write(*,*) 'D=', DD
!       call schreibeb(dSigDeps,6,6,'DSDE')
!       call schreibeb(dSigDDChi,6,1,'dSDX')
!       call schreibeb(dPiDDChi,1,1,'dPDX')
!       call schreibeb(-sig0,6,1,'dSD ')  

end subroutine umatIsoGradDam


!     #############################################################
!     Material SUBROUTINE
!     #############################################################
      subroutine umatIsoGradDam_CC(eps,Dchi,h1,hn,sig,CC,dSigDDChi,&
      Pi,dPiDDchi,gradDChi,Xi,dXiDGradDChi,nHistVars,d)                                     

!       use sharedConstants01

      
      implicit none

      include  'swvars.h'
      include  'elauto.h'
      include  'counts.h'
      include  'tdata.h'

      integer nHistVars
      double precision  sig(6),eps(6),dEps(6),DChi,Pi,gradDChi(3),Xi(3),dXiDGradDChi(3,3)
      double precision  dSigDeps(6,6),dSigDDChi(6),dPiDeps(6),dPiDDchi,d(*),CC(6,6)
      double precision,dimension(nHistVars)::hn,h1


      double precision  res(2),rightSide(2),solVec(2)
      double precision  normRes,resTol
      double precision  loctan(2,2),loctan_inv(2,2)
      double precision  macauley,heavyside

      double precision eins33(3,3),eins6(6)
      double precision DGP,lambdabar,D_acc,Yp
      double precision DGP_n,Yp_n,lambdabar_n,D_acc_n
      double precision eps_star(6),tr_eps_star
      double precision eps_star_33(3,3),eigvec(3,3),eigval(3)
      double precision eps_star_pos(3),eps_star_neg(3)
      double precision eps_star_pos_square
      double precision eps_star_neg_square,eps_star_pos_33(3,3)
      double precision trafo(3,3),eps_star_pos_6(6),eps_star_neg_33(3,3)
      double precision eps_star_neg_6(6),trafo_t(3,3),YY,qq,qqpr,phi,Hp
      double precision eins6_dyad_eins6(6,6),eins4s(3,3,3,3)
      double precision eins4sdyadNNNNi_pos(3,3,3,3)
      double precision eins4sdyadNNNNi_neg(3,3,3,3)
      double precision eins4sdyadNNNN_pos(3,3,3,3)
      double precision eins4sdyadNNNN_neg(3,3,3,3),ni(3),NNi(3,3)
      double precision NNNNi(3,3,3,3),eins4sdyadNNNNi(3,3,3,3)
      double precision nj(3),NNij(3,3),NNji(3,3),NNNNijij(3,3,3,3)
      double precision NNNNijji(3,3,3,3),NNNNjiji(3,3,3,3)
      double precision NNNNjiij(3,3,3,3),eins4sdyadNNNNijji(3,3,3,3)
      double precision eins4sdyadNNNNijij(3,3,3,3)
      double precision eins4sdyadNNNNjiji(3,3,3,3)
      double precision eins4sdyadNNNNjiij(3,3,3,3)
      double precision eins4sdyadNNNNi_pos66(6,6)
      double precision eins4sdyadNNNNi_neg66(6,6)
      double precision eins4sdyadNNNN_pos66(6,6)
      double precision eins4sdyadNNNN_neg66(6,6)
      double precision dresdeps(2,6),dresDDchi(2),dSigdD(6),krel1(2,6)
      double precision krel2stern
      double precision krel2(2),eps_0_6(6)
      double precision dSigdDdDdeps(6,6),eins4s_66(6,6)
      double precision eins4sdyadNNNNi_neg_66(6,6)
      double precision eins4sdyadNNNNi_pos_66(6,6),test1,test2
      double precision checkheavy1, checkheavy2,test
      double precision xres(2),xres_n(2),ntan(2,2)
      double precision krel1stern(6),dSigdD33(3,3),krel1stern33(3,3)
      double precision dSigdDdDdeps4(3,3,3,3),dSigdDdDdeps4_66(6,6)
      double precision difftan(2,2),eins4sdyadNNNN(6,6)

      double precision eps_star_tang(6),tr_eps_star_tang
      double precision eps_star_33_tang(3,3),eigvec_tang(3,3)
      double precision eps_star_pos_tang(3),eps_star_neg_tang(3)
      double precision eps_star_pos_square_tang,eigval_tang(3)
      double precision eps_star_neg_square_tang
      double precision trafo_tang(3,3),eps_star_pos_6_tang(6)
      double precision eps_star_neg_33_tang(3,3),eps_star_neg_6_tang(6)
      double precision trafo_t_tang(3,3),YY_tang,eps_tang(6)
      double precision eps_star_pos_33_tang(3,3),delta,dYdeps(6)
      double precision dresdeps_test(6),check_reset_1,check_reset_2
      double precision check_reset_3,tol_eq,disturb

      integer damage_indicator,critical_damage,damage_indicator_n,rot
      integer critical_damage_n,kiter,switch_ev,inelStep
      integer i,j,k,l
!----------------------Materialparameter Lukas--------------------------------------
      double precision eta_vis,IPIV(2),dblDmy,fibAng, beta, Nf(3), m33(3,3), TT(3,3)
      double precision alpha_1, alpha_2, beta_1
      integer INFO
!----------------------Materialparameter Marek--------------------------------------

      double precision:: cP1                      !0.8d4
      double precision:: linternal                !1.2d-3
      double precision:: eta                      !0.0d0

      double precision:: lambda                   !lambda (Lame constant) 5826.923077d0
      double precision:: mu                       !mu (Lame constant) 3884.615385d0
      integer::          ndam                     !ndam=1: strain equivalence, ndam=2: energy equivalence
      double precision:: eta_cc                   !crack closure parameter: eta_cc=0: no crack closure, eta_cc=1: full crack closure
      double precision:: eps_01=0.d0              !eps_0,xx (irreversible strain)
      double precision:: eps_02=0.d0              !eps_0,yy (irreversible strain)
      double precision:: eps_03=0.d0              !eps_0,zz (irreversible strain)
      double precision:: eps_04=0.d0              !sqrt(2)*eps_0,xy (irreversible strain)
      double precision:: eps_05=0.d0              !sqrt(2)*eps_0,yz (irreversible strain)
      double precision:: eps_06=0.d0              !sqrt(2)*eps_0,xz (irreversible strain)
      integer::          h_switch                 !switch for hardening: 1=linear, 2=quadratic, 3= exponential, 4= -> Infty
      double precision:: YC0                      !initial damage surface   0.0354d0
      double precision:: k1                       !hardening variable (linear term) 0.06d0
      double precision:: k2=10.d0                 !hardening variable (quadratic term)
      double precision:: k0=10.d0                 !hardening variable (exponential)
      double precision:: alpha=10.d0              !hardening variable (exponential)
      double precision:: Hpen=1.d8                !penalty parameter
      double precision:: D0                       !critical damage value
      double precision:: kappa=0.d0               !parameter for gradient extension
      integer::          loc_tang=1               !switch for local tangent, 1= analytical, 2= numerical
!       integer, parameter::          glo_tang=1               !switch for consistent tangent, 1= analytical
      
      lambda	=d(11)*d(12)/((1.d0-2.d0*d(12))*(1.d0+d(12)))	! stiffness constant Lambda
      mu	=d(11)/(2.d0*(1.d0+d(12)))			! stiffness constant mu
      fibAng	=d(14)						! angle of CFRP-Fibers
      YC0	=d(15)						! damage theshold
      k1	=d(17)						! slope of hardening q_h=k1*D
      cP1	=d(19)						! H_chi for cuppling of D and D_chi
      eta	=d(21)						! pseudoviscosity
      linternal	=d(23)						! internal length for Gradient
      eta_cc	=d(31)						! crack closure parameter: eta_cc=0: no crack closure, eta_cc=1: full crack closure
      D0	=d(32)						! maximal damage
      h_switch	=d(33)						! switch for hardening: 1=linear, 2=quadratic, 3= exponential, 4= -> Infty
      ndam	=d(34)						! ndam=1: strain equivalence, ndam=2: energy equivalence
      alpha_1	=d(35)						! parameter for fiber dependence of damage gradient alpha_1=1: complete isotropy
      alpha_2	=d(36)						! parameter for fiber dependence of damage gradient alpha_2=1: complete fiber dependence
      beta_1	=0.0d0						! parameter for activation of Crac Closure on shear terms (0=NO,1=YES)

      dblDmy=fibAng/180.d0*3.14159265359d0
      Nf=(/cos(dblDmy), sin(dblDmy) , 0.d0/)
      call product_T1dyadT1(Nf,Nf, m33)
      
      do i=1,6; dEps(i)=eps(i)-hn(i+10); end do
      eta_vis = d(21)
      
      if (h_switch .eq. 4) then
       Hp=0.d0
      else
       Hp=Hpen
      endif

!     tolerance for residuum
      resTol=1.d-8

!     tolerance for equal eigenvalues
      tol_eq=1.d-10
      disturb=1.d-8


!     unit tensor
      eins33(1:3,1:3)=0.d0
      do i=1,3
        eins33(i,i)=1.d0
      enddo

!     unit tensor in voigt notation
      eins6(1:6)=0.d0
      do i=1,3
       eins6(i)=1.d0
      enddo
      
!       Calculate Tensor T for Direction dependence of Xi
      TT = alpha_1 * m33 + alpha_2 * eins33

!     irreversible strain eps_0 in normalized voigt notation
      eps_0_6(1)=eps_01
      eps_0_6(2)=eps_02
      eps_0_6(3)=eps_03
      eps_0_6(4)=eps_04
      eps_0_6(5)=eps_05
      eps_0_6(6)=eps_06
      

!     Extraction of the internal variables
      DGP_n=hn(10)
      lambdabar_n=0.d0              !lambdadot_n=hn(15) kann auch benutzt werden, dann allerdings in manchen Fällen schlechtere Konvergenz
      damage_indicator_n=hn(17)
      D_acc_n=hn(10)
      
!     Initializing the internal variables
      DGP=DGP_n
      lambdabar=lambdabar_n
      damage_indicator=damage_indicator_n
      D_acc=D_acc_n

!     computation of the difference of eps and eps_0
      eps_star(1)=eps(1)-eps_01
      eps_star(2)=eps(2)-eps_02
      eps_star(3)=eps(3)-eps_03
      eps_star(4)=eps(4)-eps_04
      eps_star(5)=eps(5)-eps_05
      eps_star(6)=eps(6)-eps_06

!     computation of the trace of eps_star
      tr_eps_star=eps_star(1)+eps_star(2)+eps_star(3)
      call transfer_6to33_nv(eps_star, eps_star_33)
!     computation of the eigenvalues and eigenvectors of eps_star
      eigvec=eps_star_33
      call eig3(eigvec,eigval,rot)

!     computation of the positive part of eps_star
      eps_star_pos=0.d0
      eps_star_pos(1)=macauley(eigval(1))
      eps_star_pos(2)=macauley(eigval(2))
      eps_star_pos(3)=macauley(eigval(3))
      eps_star_pos_square=dot_product(eps_star_pos,eps_star_pos)

      eps_star_pos_33(1:3,1:3)=0.d0
      do i=1,3
        eps_star_pos_33(i,i)=eps_star_pos(i)
      enddo
!     transformation back to the original basis
      trafo=eigvec
      eps_star_pos_33=matmul(trafo,eps_star_pos_33)
      do i=1,3
       do j=1,3
        trafo_t(i,j)=trafo(j,i)
       enddo
      enddo
      eps_star_pos_33=matmul(eps_star_pos_33,trafo_t)
      call transfer_33to6_nv(eps_star_pos_33, eps_star_pos_6)
      

!     computation of the negative part of eps_star
      eps_star_neg=0.d0
      eps_star_neg(1)=-macauley(-eigval(1))
      eps_star_neg(2)=-macauley(-eigval(2))
      eps_star_neg(3)=-macauley(-eigval(3))
      eps_star_neg_square=dot_product(eps_star_neg,eps_star_neg)

      eps_star_neg_33(1:3,1:3)=0.d0
      do i=1,3
        eps_star_neg_33(i,i)=eps_star_neg(i)
      enddo
!     transformation back to the original basis
      eps_star_neg_33=matmul(trafo,eps_star_neg_33)
      eps_star_neg_33=matmul(eps_star_neg_33,trafo_t)
      call transfer_33to6_nv(eps_star_neg_33, eps_star_neg_6)



!    begin trial step
!    computation of the damage driving force Y
      YY=0.5d0*ndam*(1.d0-DGP)**(ndam-1.d0)*lambda*(macauley(tr_eps_star))**2.d0&
      + 0.5d0*ndam*(1.d0-(1.d0-eta_cc)*DGP)**(ndam-1.d0)*lambda*(1.d0-eta_cc)*&
      (macauley(-tr_eps_star))**2.d0 +mu*ndam*(1.d0-DGP)**(ndam-1.d0)*eps_star_pos_square&
      +mu*ndam*(1.d0-(1.d0-beta_1*eta_cc)*DGP)**(ndam-1.d0)*(1.d0-beta_1*eta_cc)*eps_star_neg_square&
      - Hp*macauley(DGP-D0)+ kappa + cP1*(DChi-DGP)

!     computation of the hardening variable qq
      if (h_switch .eq. 1) then
       qq=k1*D_acc
      elseif (h_switch .eq. 2) then
       qq=k1*D_acc+k2*D_acc**2.d0
      elseif (h_switch .eq. 3) then
       qq=k0*(1.d0-dexp(-alpha*D_acc))
      elseif (h_switch .eq. 4) then
       call getQhAndQhPr(D_acc,qq,qqpr,k1,D0)
       qq=-qq
       qqpr=-qqpr
      end if

!     computation of the damage function
      phi=YY-(YC0+qq)
      
      
      inelStep=1
!     damage function is evaluated
      if (phi.le. 1.d-8) then
        inelStep=0
      endif

      if(duVnshSW.eq.1) then
        inelStep=int(hn(17))!old step inelastic or not
      endif



      if (inelStep.eq.0) then
      !elastic step
!      write(*,*) 'elastic step'
      damage_indicator=0

!     computation of the consistent tangent
      call product_T1dyadT1_nv(eins6,eins6,eins6_dyad_eins6)
      call sopertion_T2sT2(eins33,eins33,eins4s)
      call transfer_4_to_66(eins4s,eins4s_66,5)

      eins4sdyadNNNNi_pos(1:3,1:3,1:3,1:3)=0.d0
      eins4sdyadNNNNi_neg(1:3,1:3,1:3,1:3)=0.d0
      eins4sdyadNNNN_pos(1:3,1:3,1:3,1:3)=0.d0
      eins4sdyadNNNN_neg(1:3,1:3,1:3,1:3)=0.d0

!     check for equal eigenvalues (if equal: perturbation)
      switch_ev=0
      check_reset_1=0.d0
      check_reset_2=0.d0
      check_reset_3=0.d0


      if (dabs(eigval(1)-eigval(2)).lt.(tol_eq)) then
       switch_ev=1
       check_reset_1=1.d0 
       eigval(1)=eigval(1)+disturb
       if (dabs(eigval(2)-eigval(3)).lt.(tol_eq)) then
        check_reset_2=1.d0 
        eigval(2)=eigval(2)-disturb
       endif
       if (dabs(eigval(1)-eigval(3)).lt.(tol_eq)) then
        check_reset_3=1.d0 
        eigval(1)=eigval(1)+disturb
       endif
      elseif (dabs(eigval(1)-eigval(3)).lt.(tol_eq)) then
       switch_ev=2
       check_reset_1=1.d0 
       eigval(1)=eigval(1)+disturb
       if (dabs(eigval(3)-eigval(2)).lt.(tol_eq)) then
        check_reset_2=1.d0 
        eigval(3)=eigval(3)-disturb
       endif
       if (dabs(eigval(1)-eigval(2)).lt.(tol_eq)) then
        check_reset_3=1.d0 
        eigval(1)=eigval(1)+disturb
       endif
      elseif (dabs(eigval(2)-eigval(3)).lt.(tol_eq)) then
       switch_ev=3
       check_reset_1=1.d0 
       eigval(2)=eigval(2)+disturb
       if (dabs(eigval(3)-eigval(1)).lt.(tol_eq)) then
        check_reset_2=1.d0 
        eigval(3)=eigval(3)-disturb
       endif
       if (dabs(eigval(2)-eigval(1)).lt.(tol_eq)) then
        check_reset_3=1.d0 
        eigval(2)=eigval(2)+disturb
       endif
      endif

!     error messages for not succesful perturbations
      if (dabs(eigval(1)-eigval(2)).lt.(tol_eq)) then
        write(*,*) 'Perturbation of the eigenvalues does not work!!!' 
      elseif (dabs(eigval(1)-eigval(3)).lt.(tol_eq)) then
        write(*,*) 'Perturbation of the eigenvalues does not work!!!' 
      elseif (dabs(eigval(2)-eigval(3)).lt.(tol_eq)) then
        write(*,*) 'Perturbation of the eigenvalues does not work!!!' 
      endif



      do i=1,3
       ni(1:3)=eigvec(1:3,i) !small n^star (eigenvector of eps_star)
       call product_T1dyadT1(ni,ni, NNi) !NNi: capital N^star (dyadic product of the small n^stars)
       call product_T2dyadT2(NNi,NNi, NNNNi) !NNNNi: dyadic product of the capital N^stars)
       call product_T4xxT4(eins4s,NNNNi, eins4sdyadNNNNi)

       if (heavyside(eigval(i)).gt.0.5d0) then                            !ensures that for eigval(i)=0 stiffness terms are not accounted in eins4sdyadNNNNi_pos and eins4sdyadNNNNi_neg at the same time
         checkheavy1=1.d0
         checkheavy2=0.d0
       else
         checkheavy1=0.d0
         checkheavy2=1.d0  
       endif

       eins4sdyadNNNNi_pos(1:3,1:3,1:3,1:3)=eins4sdyadNNNNi_pos(1:3,1:3,1:3,1:3)&
       + eins4sdyadNNNNi(1:3,1:3,1:3,1:3)*checkheavy1
    

      call transfer_4_to_66(eins4sdyadNNNNi_pos,eins4sdyadNNNNi_pos_66,5)     




       eins4sdyadNNNNi_neg(1:3,1:3,1:3,1:3)=eins4sdyadNNNNi_neg(1:3,1:3,1:3,1:3)&
       + eins4sdyadNNNNi(1:3,1:3,1:3,1:3)*checkheavy2

      call transfer_4_to_66(eins4sdyadNNNNi_neg,eins4sdyadNNNNi_neg_66,5)

       do j=1,3
        if (j.eq.i) cycle
        nj(1:3)=eigvec(1:3,j) !small n^star (eigenvector of eps_star)
        call product_T1dyadT1(ni,nj, NNij) !NNij: capital N^star (dyadic product of the small n^stars)
        call product_T1dyadT1(nj,ni, NNji) !NNji: capital N^star (dyadic product of the small n^stars)  

        call product_T2dyadT2(NNij,NNij, NNNNijij) !NNNN: dyadic product of the capital N^stars)
        call product_T2dyadT2(NNij,NNji, NNNNijji) !NNNN: dyadic product of the capital N^stars)
        call product_T2dyadT2(NNji,NNij, NNNNjiij) !NNNN: dyadic product of the capital N^stars)
        call product_T2dyadT2(NNji,NNji, NNNNjiji) !NNNN: dyadic product of the capital N^stars)
        call product_T4xxT4(eins4s,NNNNijij, eins4sdyadNNNNijij)
        call product_T4xxT4(eins4s,NNNNijji, eins4sdyadNNNNijji)
        call product_T4xxT4(eins4s,NNNNjiij, eins4sdyadNNNNjiij)
        call product_T4xxT4(eins4s,NNNNjiji, eins4sdyadNNNNjiji)

        eins4sdyadNNNN_pos(1:3,1:3,1:3,1:3)=eins4sdyadNNNN_pos(1:3,1:3,1:3,1:3)&
        + 0.5d0*checkheavy1*eigval(i)/(eigval(i)-eigval(j))*&
        (eins4sdyadNNNNijij(1:3,1:3,1:3,1:3)+eins4sdyadNNNNijji(1:3,1:3,1:3,1:3)&
        +eins4sdyadNNNNjiij(1:3,1:3,1:3,1:3)+eins4sdyadNNNNjiji(1:3,1:3,1:3,1:3))

        eins4sdyadNNNN_neg(1:3,1:3,1:3,1:3)=eins4sdyadNNNN_neg(1:3,1:3,1:3,1:3)&
        + 0.5d0*checkheavy2*eigval(i)/(eigval(i)-eigval(j))*&
        (eins4sdyadNNNNijij(1:3,1:3,1:3,1:3)+eins4sdyadNNNNijji(1:3,1:3,1:3,1:3)&
        +eins4sdyadNNNNjiij(1:3,1:3,1:3,1:3)+eins4sdyadNNNNjiji(1:3,1:3,1:3,1:3))

      call transfer_4_to_66(eins4sdyadNNNN_pos,eins4sdyadNNNN_pos66,5)
      call transfer_4_to_66(eins4sdyadNNNN_neg,eins4sdyadNNNN_neg66,5)


       enddo
      enddo

! correction of the eigenvalues, which were perturbated before because of equal eigenvalues
      if (switch_ev.eq.1) then
       eigval(1)=eigval(1)-tol_eq*check_reset_1-tol_eq*check_reset_3
       eigval(2)=eigval(2)+tol_eq*check_reset_2
      elseif (switch_ev.eq.2) then
       eigval(1)=eigval(1)-tol_eq*check_reset_1-tol_eq*check_reset_3
       eigval(3)=eigval(3)+tol_eq*check_reset_2
      elseif (switch_ev.eq.3) then
       eigval(2)=eigval(2)-tol_eq*check_reset_1-tol_eq*check_reset_3
       eigval(3)=eigval(3)+tol_eq*check_reset_2
      endif



      call transfer_4_to_66(eins4sdyadNNNNi_pos,eins4sdyadNNNNi_pos66,5)
      call transfer_4_to_66(eins4sdyadNNNNi_neg,eins4sdyadNNNNi_neg66,5)
      call transfer_4_to_66(eins4sdyadNNNN_pos,eins4sdyadNNNN_pos66,5)
      call transfer_4_to_66(eins4sdyadNNNN_neg,eins4sdyadNNNN_neg66,5)

      eins4sdyadNNNN=eins4sdyadNNNN_pos66+eins4sdyadNNNN_neg66
     


!    final result of the consistent tangents for an elastic step
      dSigDeps(1:6,1:6)=(1.d0-DGP)**ndam*lambda*heavyside(tr_eps_star)*&
       eins6_dyad_eins6(1:6,1:6) +(1.d0-(1.d0-eta_cc)*DGP)**ndam*lambda*&
       (1.d0-heavyside(tr_eps_star))*eins6_dyad_eins6(1:6,1:6) +2.d0*mu*&
       (1.d0-DGP)**ndam*eins4sdyadNNNNi_pos66(1:6,1:6) + 2.d0*mu*(1.d0-&
       (1.d0-beta_1*eta_cc)*DGP)**ndam*eins4sdyadNNNNi_neg66(1:6,1:6) +2.d0*mu*&
       (1.d0-DGP)**ndam*eins4sdyadNNNN_pos66(1:6,1:6) + 2.d0*mu*(1.d0-&
       (1.d0-beta_1*eta_cc)*DGP)**ndam*eins4sdyadNNNN_neg66(1:6,1:6)


      dSigDDchi(1:6)=0.d0
      dPiDEps(1:6)=0.d0
      dPiDDchi=cP1 !+ d(17) !MF: added hardening for micromorphic variable
!     end computation of the consistent tangent


!**************************************
!     damage step      
!**************************************     
      else !if:damage
!       write(*,*) 'damage step'
       damage_indicator=1

       do kiter=1,30 !kiter-loop

!    computation of the damage driving force Y
      YY=0.5d0*ndam*(1.d0-DGP)**(ndam-1.d0)*lambda*(macauley(tr_eps_star))&
      **2.d0+0.5d0*ndam*(1.d0-(1.d0-eta_cc)*DGP)**(ndam-1.d0)*lambda*&
      (1.d0-eta_cc)*(macauley(-tr_eps_star))**2.d0 +mu*ndam*(1.d0-DGP)**&
      (ndam-1.d0)*eps_star_pos_square +mu*ndam*(1.d0-(1.d0-beta_1*eta_cc)*DGP)**&
      (ndam-1.d0)*(1.d0-beta_1*eta_cc)*eps_star_neg_square - Hp*macauley(DGP-D0)+ kappa + cP1*(DChi-DGP)

!     computation of the hardening variable qq
       if (h_switch .eq. 1) then
        qq=k1*D_acc
       elseif (h_switch .eq. 2) then
        qq=k1*D_acc+k2*D_acc**2.d0
       elseif (h_switch .eq. 3) then
        qq=k0*(1.d0-dexp(-alpha*D_acc))
       elseif (h_switch .eq. 4) then
        call getQhAndQhPr(D_acc,qq,qqpr,k1,D0)
        qq=-qq
        qqpr=-qqpr
       end if

       call compute_residuum01(DGP,DGP_n,lambdabar_n,lambdabar,YY,qq,YC0,res)

!      compute error 
       call normVec(res,2, normRes)

!      criterion for convergence
       if (normRes.lt.1.d-8) exit

!      too many iterations?
       if(kiter.gt.10) then; rmeas=1.5d0; exit; endif

       do i = 1,2
        rightSide(i) = -1.d0*res(i)		    
       end do

!-----------------------------------------------------------------------
!      local tangent (analytically)
       if (loc_tang .eq. 1) then
       loctan(1,1)=1.d0
       loctan(1,2)=-1.d0
       loctan(2,1)=-0.5d0*ndam*(ndam-1.d0)*(1.d0-DGP)**(ndam-2.d0)*&
       lambda*(macauley(tr_eps_star))**2.d0 -0.5d0*ndam*(ndam-1.d0)&
       *(1.d0-(1.d0-eta_cc)*DGP)**(ndam-2.d0)*lambda*(1.d0-eta_cc)&
       **(2.d0)*(macauley(-tr_eps_star))**2.d0-mu*ndam*(ndam-1.d0)&
       *(1.d0-DGP)**(ndam-2.d0)*eps_star_pos_square-mu*ndam*&
       (ndam-1.d0)*(1.d0-(1.d0-beta_1*eta_cc)*DGP)**(ndam-2.d0)*(1.d0-beta_1*eta_cc)&
       **(2.d0)*eps_star_neg_square-Hp*heavyside(DGP-D0) - cP1
       
       if (h_switch .eq. 1) then
        loctan(2,2)=-k1
       elseif (h_switch .eq. 2) then
        loctan(2,2)=-(k1+2.d0*k2*D_acc)
       elseif (h_switch .eq. 3) then
        loctan(2,2)=-k0*alpha*dexp(-alpha*D_acc)
       elseif (h_switch .eq. 4) then
        call getQhAndQhPr(D_acc,dblDmy,loctan(2,2),k1,D0)
       end if
       end if
!      end local tangent (analytically)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!      local tangent (numerically)
       if (loc_tang .eq. 2) then
!         xres(1)=DGP
!         xres(2)=lambdabar
!         xres_n(1)=DGP_n
!         xres_n(2)=lambdabar_n
! 
!         call numTangent01(xres,xres_n,res,tr_eps_star,eps_star_pos_square,eps_star_neg_square,Dchi,ntan,D0)
! !         do i=1,2
! !          do j=1,2
! !           difftan(i,j)=loctan(i,j)-ntan(i,j)
! !          enddo
! !         enddo
!         loctan=ntan
        write(*,*)'!!!FOR NUMERICAL TANGENT, THE CODE HAS TO BE ADJUSTED!!!'
        end if
!      end local tangent (numerically)

!-----------------------------------------------------------------------

!      solution of the linear equation system: LP-neu: Lapac Routine verwendet
       solVec=rightSide

       call DGESV( 2, 1, loctan, 2, IPIV, solVec, 2, INFO)


       if (INFO.NE.0) then
         write(*,*)"WARNING: INFO from DGESV=", INFO

         
       elseif (INFO.eq.0) then
!      update of the internal variables
         DGP=DGP+solVec(1)
         lambdabar=lambdabar+solVec(2)
         D_acc=D_acc_n+lambdabar
       endif

       enddo !kiter-Loop   

!     computation of the consistent tangent
      call product_T1dyadT1_nv(eins6,eins6,eins6_dyad_eins6)
      call sopertion_T2sT2(eins33,eins33,eins4s)

      eins4sdyadNNNNi_pos(1:3,1:3,1:3,1:3)=0.d0
      eins4sdyadNNNNi_neg(1:3,1:3,1:3,1:3)=0.d0
      eins4sdyadNNNN_pos(1:3,1:3,1:3,1:3)=0.d0
      eins4sdyadNNNN_neg(1:3,1:3,1:3,1:3)=0.d0

!     check for equal eigenvalues (if equal: perturbation)
      switch_ev=0
      check_reset_1=0.d0
      check_reset_2=0.d0
      check_reset_3=0.d0


      if (dabs(eigval(1)-eigval(2)).lt.(tol_eq)) then
       switch_ev=1
       check_reset_1=1.d0 
       eigval(1)=eigval(1)+disturb
       if (dabs(eigval(2)-eigval(3)).lt.(tol_eq)) then
        check_reset_2=1.d0 
        eigval(2)=eigval(2)-disturb
       endif
       if (dabs(eigval(1)-eigval(3)).lt.(tol_eq)) then
        check_reset_3=1.d0 
        eigval(1)=eigval(1)+disturb
       endif
      elseif (dabs(eigval(1)-eigval(3)).lt.(tol_eq)) then
       switch_ev=2
       check_reset_1=1.d0 
       eigval(1)=eigval(1)+disturb
       if (dabs(eigval(3)-eigval(2)).lt.(tol_eq)) then
        check_reset_2=1.d0 
        eigval(3)=eigval(3)-disturb
       endif
       if (dabs(eigval(1)-eigval(2)).lt.(tol_eq)) then
        check_reset_3=1.d0 
        eigval(1)=eigval(1)+disturb
       endif
      elseif (dabs(eigval(2)-eigval(3)).lt.(tol_eq)) then
       switch_ev=3
       check_reset_1=1.d0 
       eigval(2)=eigval(2)+disturb
       if (dabs(eigval(3)-eigval(1)).lt.(tol_eq)) then
        check_reset_2=1.d0 
        eigval(3)=eigval(3)-disturb
       endif
       if (dabs(eigval(2)-eigval(1)).lt.(tol_eq)) then
        check_reset_3=1.d0 
        eigval(2)=eigval(2)+disturb
       endif
      endif

!     error messages for not succesful perturbations
      if (dabs(eigval(1)-eigval(2)).lt.(tol_eq)) then
        write(*,*) 'Perturbation of the eigenvalues does not work!!!' 
      elseif (dabs(eigval(1)-eigval(3)).lt.(tol_eq)) then
        write(*,*) 'Perturbation of the eigenvalues does not work!!!' 
      elseif (dabs(eigval(2)-eigval(3)).lt.(tol_eq)) then
        write(*,*) 'Perturbation of the eigenvalues does not work!!!' 
      endif

      do i=1,3
       ni(1:3)=eigvec(1:3,i) !small n^star (eigenvector of eps_star)
       call product_T1dyadT1(ni,ni, NNi) !NNi: capital N^star (dyadic product of the small n^stars) 
       call product_T2dyadT2(NNi,NNi, NNNNi) !NNNNi: dyadic product of the capital N^stars)
       call product_T4xxT4(eins4s,NNNNi, eins4sdyadNNNNi)

       if (heavyside(eigval(i)).gt.0.5d0) then                            !ensures that for eigval(i)=0 stiffness terms are not accounted in eins4sdyadNNNNi_pos and eins4sdyadNNNNi_neg at the same time
         checkheavy1=1.d0
         checkheavy2=0.d0
       else
         checkheavy1=0.d0
         checkheavy2=1.d0  
       endif

       eins4sdyadNNNNi_pos(1:3,1:3,1:3,1:3)=eins4sdyadNNNNi_pos(1:3,1:3,1:3,1:3)&
       + eins4sdyadNNNNi(1:3,1:3,1:3,1:3)*checkheavy1

       eins4sdyadNNNNi_neg(1:3,1:3,1:3,1:3)=eins4sdyadNNNNi_neg(1:3,1:3,1:3,1:3)&
       + eins4sdyadNNNNi(1:3,1:3,1:3,1:3)*checkheavy2

       do j=1,3
        if (j.eq.i) cycle
        nj(1:3)=eigvec(1:3,j) !small n^star (eigenvector of eps_star)
        call product_T1dyadT1(ni,nj, NNij) !NNij: capital N^star (dyadic product of the small n^stars)
        call product_T1dyadT1(nj,ni, NNji) !NNji: capital N^star (dyadic product of the small n^stars)  

        call product_T2dyadT2(NNij,NNij, NNNNijij) !NNNN: dyadic product of the capital N^stars)
        call product_T2dyadT2(NNij,NNji, NNNNijji) !NNNN: dyadic product of the capital N^stars)
        call product_T2dyadT2(NNji,NNij, NNNNjiij) !NNNN: dyadic product of the capital N^stars)
        call product_T2dyadT2(NNji,NNji, NNNNjiji) !NNNN: dyadic product of the capital N^stars)
        call product_T4xxT4(eins4s,NNNNijij, eins4sdyadNNNNijij)
        call product_T4xxT4(eins4s,NNNNijji, eins4sdyadNNNNijji)
        call product_T4xxT4(eins4s,NNNNjiij, eins4sdyadNNNNjiij)
        call product_T4xxT4(eins4s,NNNNjiji, eins4sdyadNNNNjiji)

        eins4sdyadNNNN_pos(1:3,1:3,1:3,1:3)=eins4sdyadNNNN_pos(1:3,1:3,1:3,1:3)&
           + 0.5d0*checkheavy1*eigval(i)/(eigval(i)-eigval(j))*(eins4sdyadNNNNijij(1:3,1:3,1:3,1:3)&
           +eins4sdyadNNNNijji(1:3,1:3,1:3,1:3)+eins4sdyadNNNNjiij(1:3,1:3,1:3,1:3)+eins4sdyadNNNNjiji(1:3,1:3,1:3,1:3))

        eins4sdyadNNNN_neg(1:3,1:3,1:3,1:3)=eins4sdyadNNNN_neg(1:3,1:3,1:3,1:3)&
         + 0.5d0*checkheavy2*eigval(i)/(eigval(i)-eigval(j))*(eins4sdyadNNNNijij(1:3,1:3,1:3,1:3)+&
         eins4sdyadNNNNijji(1:3,1:3,1:3,1:3)+eins4sdyadNNNNjiij(1:3,1:3,1:3,1:3)+ eins4sdyadNNNNjiji(1:3,1:3,1:3,1:3))

       enddo
      enddo

! correction of the eigenvalues, which were perturbated before because of equal eigenvalues
      if (switch_ev.eq.1) then
       eigval(1)=eigval(1)-tol_eq*check_reset_1-tol_eq*check_reset_3
       eigval(2)=eigval(2)+tol_eq*check_reset_2
      elseif (switch_ev.eq.2) then
       eigval(1)=eigval(1)-tol_eq*check_reset_1-tol_eq*check_reset_3
       eigval(3)=eigval(3)+tol_eq*check_reset_2
      elseif (switch_ev.eq.3) then
       eigval(2)=eigval(2)-tol_eq*check_reset_1-tol_eq*check_reset_3
       eigval(3)=eigval(3)+tol_eq*check_reset_2
      endif


      call transfer_4_to_66(eins4sdyadNNNNi_pos,eins4sdyadNNNNi_pos66,5)
      call transfer_4_to_66(eins4sdyadNNNNi_neg,eins4sdyadNNNNi_neg66,5)
      call transfer_4_to_66(eins4sdyadNNNN_pos,eins4sdyadNNNN_pos66,5)
      call transfer_4_to_66(eins4sdyadNNNN_neg,eins4sdyadNNNN_neg66,5)

!    
      dSigDeps(1:6,1:6)=(1.d0-DGP)**ndam*lambda*heavyside(tr_eps_star)*eins6_dyad_eins6(1:6,1:6)&
      +(1.d0-(1.d0-eta_cc)*DGP)**ndam*lambda*(1.d0-heavyside(tr_eps_star))*eins6_dyad_eins6(1:6,1:6)&
      +2.d0*mu*(1.d0-DGP)**ndam*eins4sdyadNNNNi_pos66(1:6,1:6) + 2.d0*mu*(1.d0-(1.d0-beta_1*eta_cc)*DGP)**&
      ndam*eins4sdyadNNNNi_neg66(1:6,1:6)+ 2.d0*mu*(1.d0-DGP)**ndam*eins4sdyadNNNN_pos66(1:6,1:6)&
      + 2.d0*mu*(1.d0-(1.d0-beta_1*eta_cc)*DGP)**ndam*eins4sdyadNNNN_neg66(1:6,1:6)
!       call schreibeb(dSigDeps,6,6,'sti1')

!    computation of dres_deps
      dresdeps(1,1:6)=0.d0                                                  !dres1deps
      dresdeps(2,1:6)=ndam*(1.d0-DGP)**(ndam-1.d0)*lambda  *macauley(tr_eps_star)*&    !dres2deps=dYdeps
         eins6(1:6)-ndam*(1.d0-(1.d0-eta_cc)*DGP)**(ndam-1.d0)*(1.d0-eta_cc)*&
         lambda*macauley(-tr_eps_star)*eins6(1:6) + 2.d0*mu*ndam*(1.d0-DGP)**&
         (ndam-1.d0)*eps_star_pos_6(1:6) + 2.d0*mu*ndam*(1.d0-(1.d0-beta_1*eta_cc)*DGP)**&
         (ndam-1.d0)*(1.d0-eta_cc)*eps_star_neg_6(1:6)

!    computation of dres_dDChi
      dresdDChi(1)=0.d0
      dresdDChi(2)=cP1

!     computation of dsigmadD
       dSigdD(1:6)=-ndam*(1.d0-DGP)**(ndam-1.d0)*lambda*macauley(tr_eps_star)*eins6(1:6)&
       + ndam*(1.d0-(1.d0-eta_cc)*DGP)**(ndam-1.d0)*(1.d0-eta_cc)*lambda*macauley(-tr_eps_star)*&
       eins6(1:6) - 2.d0*mu*ndam*(1.d0-DGP)**(ndam-1.d0)*eps_star_pos_6(1:6) - 2.d0*mu*ndam*(1.d0&
       -(1.d0-beta_1*eta_cc)*DGP)**(ndam-1.d0)*(1.d0-beta_1*eta_cc)*eps_star_neg_6(1:6)

!      new computation of kloc for converged values
!-----------------------------------------------------------------------
!      local tangent (analytically)
       if (loc_tang .eq. 1) then
       loctan(1,1)=1.d0
       loctan(1,2)=-1.d0
       loctan(2,1)=-0.5d0*ndam*(ndam-1.d0)*(1.d0-DGP)**(ndam-2.d0)*lambda*(macauley(tr_eps_star))**2.d0&
          -0.5d0*ndam*(ndam-1.d0)*(1.d0-(1.d0-eta_cc)*DGP)**(ndam-2.d0)*lambda*(1.d0-eta_cc)**(2.d0)*&
          (macauley(-tr_eps_star))**2.d0 -mu*ndam*(ndam-1.d0)*(1.d0-DGP)**(ndam-2.d0)*eps_star_pos_square&
          -mu*ndam*(ndam-1.d0)*(1.d0-(1.d0-beta_1*eta_cc)*DGP)**(ndam-2.d0)*(1.d0-beta_1*eta_cc)**(2.d0)*eps_star_neg_square&
          -Hp*heavyside(DGP-D0) - cP1
          
       if (h_switch .eq. 1) then
        loctan(2,2)=-k1
       elseif (h_switch .eq. 2) then
        loctan(2,2)=-(k1+2.d0*k2*D_acc)
       elseif (h_switch .eq. 3) then
        loctan(2,2)=-k0*alpha*dexp(-alpha*D_acc)
       elseif (h_switch .eq. 4) then
        call getQhAndQhPr(D_acc,dblDmy,loctan(2,2),k1,D0)
       end if
       end if
!      end local tangent (analytically)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!      local tangent (numerically)
       if (loc_tang .eq. 2) then
        xres(1)=DGP
        xres(2)=lambdabar
        xres_n(1)=DGP_n
        xres_n(2)=lambdabar_n

        call numTangent01(xres,xres_n,res,tr_eps_star,eps_star_pos_square,eps_star_neg_square,Dchi, ntan)
!         do i=1,2
!          do j=1,2
!           difftan(i,j)=loctan(i,j)-ntan(i,j)
!          enddo
!         enddo
        loctan=ntan
        end if
!       call schreibeb(difftan,2,2,'dtan')
!      end local tangent (numerically)
!-----------------------------------------------------------------------


       call invert(loctan,2,2)
       loctan_inv(1:2,1:2)=-loctan(1:2,1:2)
       krel1=matmul(loctan_inv,dresdeps)
       krel1stern(1:6)=krel1(1,1:6)
       do i=1,6
        do j=1,6
          dSigdDdDdeps(i,j)=dSigdD(i)*krel1stern(j)
        enddo
       enddo

 
       krel2=matmul(loctan_inv,dresdDchi)
       krel2stern=krel2(1)
!     final result of the consistent tangents for a damage step
       dSigDeps(1:6,1:6)=dSigDeps(1:6,1:6) + dSigdDdDdeps(1:6,1:6)

       do i=1,6
        dSigDDchi(i)=dSigdD(i)*krel2stern
        dPiDeps(i)=-cP1*krel1stern(i) 
       enddo

       dPiDDchi=cP1*(1.d0-krel2stern) !+ d(17) !MF: added hardening for micromorphic variable

       
!     end computation of the consistent tangent

      end if !if:damage

!     stress computation
      sig(1:6)=sig(1:6)+(1.d0-DGP)**ndam*lambda*(macauley(tr_eps_star))&
      *eins6(1:6)+2.d0*mu*(1.d0-DGP)**ndam*eps_star_pos_6(1:6) - (1.d0-&
      (1.d0-eta_cc)*DGP)**ndam*lambda*(macauley(-tr_eps_star))*eins6(1:6)&
      +2.d0*mu*(1.d0-(1.d0-beta_1*eta_cc)*DGP)**ndam*eps_star_neg_6(1:6)+lambda&
      *(eps_01+eps_02+eps_03)*eins6(1:6)+2.d0*mu*eps_0_6(1:6)
      

       
     
!----------LP-Viscosity-----------------------
     sig(1:6)=sig(1:6)+eta_vis/dt*dEps(1:6)
     
     do i=1,6
      dSigDeps(i,i)=dSigDeps(i,i) + eta_vis/dt
     enddo
!----------LP-Viscosity-END----------------------
      CC(1:6,1:6) = CC(1:6,1:6) + dSigDeps(1:6,1:6)

!     computation of Pi
      Pi=cP1*(DChi-DGP) !+ d(17)*DChi   !MF: added hardening for micromorphic variable
      Xi(1:3)=Xi(1:3)+d(11)*(linternal**2)*gradDChi(1:3)
      
      
      dXiDGradDChi=dXiDGradDChi+d(11)*(linternal**2)*TT
      
      Xi=matmul(TT,Xi)

      h1(2:7)  = sig(1:6)
      h1(8)    = qq
      h1(9)    = Pi
      h1(10)   = DGP
      h1(11:16)= eps(1:6)
      h1(17)   = damage_indicator

      RETURN
      END      








!-------------------------------------------------------------
function heavyside(a) result(b)

      implicit none 
      double precision, intent(in):: a
      double precision b 

      b = sign(0.5d0,a) + 0.5
 
      return
end

SUBROUTINE matinv2(A,B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    double precision, intent(in) :: A(2,2)   !! Matrix
    double precision, intent(out):: B(2,2)   !! Inverse matrix
    double precision             :: detinv
    
    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
!     write(*,*) 'Inverse_det', detinv

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
!     write(*,*) 'Inverse'
end
!-------------------------------------------------------------
!-------------------------------------------------------------
SUBROUTINE GAUSS(a,b,n,x)
!
!  GAUßsches Eliminationsverfahren mit Spaltenpivotisierung
!        [BRAESS, Dietrich]
! 
!      n       : Größe des GS
!      a(n,n)  : "Steifigkeits"-Matrix, wird überschrieben
!      b(n)    : rechte Seite, wird überschrieben
!      x(n)    : Lösungsvektor
! 
      IMPLICIT NONE
      INTEGER n,I,J,K,M,ind(n)
      DOUBLE PRECISION a(n,n),a_(n,n),b(n),b_(n),d(n),d_(n),x(n)
      DOUBLE PRECISION Zeilensumme
!  Diagonalmatrix d, die a zeilenweise äquilibriert

      DO I=1,n
          x(I)=0.D0
          Zeilensumme=0.D0
          DO K=1,n
              Zeilensumme=Zeilensumme+Abs(a(I,K))
          END DO

          if (Zeilensumme .lt. 1.d-10) then !ADDED BY VIVIAN TINI
             d(I)=1.d15 !ADDED BY VIVIAN TINI
          else if (Zeilensumme .ge. 1.d-10) then !ADDED BY VIVIAN TINI
             d(I)=1.D0/Zeilensumme ! Original line
          end if

      END DO

      DO J=1,n
          M=J
          DO I=J+1,n
              IF (Abs(d(I)*a(I,J)).GT.Abs(d(M)*a(M,J))) THEN
                  M=I
              END IF
          END DO
          ind(J)=m
          IF (a(M,J).EQ.0.D0) THEN
              STOP 'Abbruch GAUSS'
              write(*,*)'Abbruch GAUSS'
!               CALL XIT  ! ADDED BY VIVIAN TINI (Calling user subroutine XIT to ensure proper exit)
          END IF
          DO K=J,n
              a_(J,K)=a(J,K)
          END DO
          DO K=J,n
              a(J,K)=a(M,K)
          END DO
          DO K=J,n
              a(M,K)=a_(J,K)
          END DO
          d_(J)=d(J)
          d(J)=d(M)
          d(M)=d_(J)
          DO I=J+1,n
!               write(*,*) 'before div a JJ'
              a(I,J)=a(I,J)/a(J,J)
!               write(*,*) 'after div aJJ'
              DO K=j+1,n
                  a(I,K)=a(I,K)-a(I,J)*a(J,K)
              END DO
          END DO
      END DO
      RETURN
END
!-------------------------------------------------------------
!-------------------------------------------------------------
subroutine normVec(v,n, norm) 

      implicit none
      integer  i,n
      double precision v(n),norm

      norm = 0.d0
      do i=1,n
         norm = norm + v(i)*v(i)
      end do	
      norm = sqrt(norm)
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine compute_residuum01(DGP,DGP_n,lambdabar_n,lambdabar,YY,qq,YC0,res)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc]
!       use sharedConstants01
      implicit none

      double precision DGP,DGP_n,lambdabar_n,lambdabar,YY,qq,YC0,res(2)

!c...  computes the residuum vector

!     compute residuum 1
      res(1)=DGP-DGP_n-lambdabar
!       if (sign(1.d0,YY).lt.(-1d-8)) then
!        write(*,*) 'Warnung: sign(YY) ist negativ!!!!!!!!'
!       endif

!     compute residuum 2
      res(2)=YY-(YC0+qq)

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc]
!-------------------------------------------------------------
subroutine product_T4xxT4(aa4,bb4, res4)
      implicit none

      integer i,j,k,l,m,n
      real*8  aa4(3,3,3,3),bb4(3,3,3,3)!input
      real*8  res4(3,3,3,3)!output

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              res4(i,j,k,l)=0.d0
              do m=1,3
                do n=1,3
                  res4(i,j,k,l)=res4(i,j,k,l)+aa4(i,j,m,n)*bb4(m,n,k,l)
                end do
              end do
            end do
          end do
        end do
      end do

      return
      end 
!-------------------------------------------------------------
!-------------------------------------------------------------
subroutine product_T2dyadT2(aa33,bb33, res4)
      implicit none

      integer i,j,k,l
      real*8  aa33(3,3),bb33(3,3)!input
      real*8  res4(3,3,3,3)!output


      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              res4(i,j,k,l)=aa33(i,j)*bb33(k,l)
            end do
          end do
        end do
      end do

      return
      end
!-------------------------------------------------------------
!-------------------------------------------------------------
subroutine transfer_4_to_66(A4,A66,switch)

      implicit none

      integer switch
  
      double precision :: A4(3,3,3,3),A66(6,6)

! Local variables
       
      integer i,j,k,l,m,n
 
      double precision :: factor1,factor2,a,b

! Define factors according to preferred representation

      if (switch .EQ. 1) then

! Contravariant : A^(ab)

      factor1 = 1.0D+00
      factor2 = 1.0D+00

      elseif (switch .EQ. 2) then

      ! Mixed type : A^a_b

      factor1 = 1.0D+00
      factor2 = 2.0D+00

      elseif (switch .EQ. 3) then

      ! Mixed type : A_a^b

      factor1 = 2.0D+00
      factor2 = 1.0D+00

      elseif (switch .EQ. 4) then

! Covariant : A_(ab)

      factor1 = 2.0D+00
      factor2 = 2.0D+00

      elseif (switch .EQ. 5) then

! Normalized : A^(ab) = A^a_b = A_a^b = A^(ab)

      factor1 = sqrt(2.0D+00)
      factor2 = sqrt(2.0D+00)

      else
  
      write(*,*) 'Wrong switch argument in transfer_4_to_66!'
      stop

      endif  

      do m = 1,6
      if (m .EQ. 1) then
       i = 1
       j = 1
       a = 1.0D+00
      elseif (m .EQ. 2) then
       i = 2
       j = 2
       a = 1.0D+00
      elseif (m .EQ. 3) then
       i = 3
       j = 3
       a = 1.0D+00
      elseif (m .EQ. 4) then
       i = 1
       j = 2
       a = factor1
      elseif (m .EQ. 5) then
       i = 2
       j = 3
       a = factor1
      elseif (m .EQ. 6) then
        i = 3
        j = 1
        a = factor1
      endif     
     
      do n = 1,6
       if (n .EQ. 1) then
	 k = 1
         l = 1
         b = 1.0D+00
       elseif (n .EQ. 2) then
         k = 2
         l = 2
         b = 1.0D+00
       elseif (n .EQ. 3) then
         k = 3
         l = 3
         b = 1.0D+00
       elseif (n .EQ. 4) then
         k = 1
         l = 2
         b = factor2
       elseif (n .EQ. 5) then
         k = 2
         l = 3
         b = factor2
       elseif (n .EQ. 6) then
         k = 3
         l = 1
         b = factor2
       endif
        
       A66(m,n) = a*b*A4(i,j,k,l)
        
      end do

      end do

      end subroutine transfer_4_to_66
!-------------------------------------------------------------
!-------------------------------------------------------------
SUBROUTINE sopertion_T2sT2(aa33, bb33, res4)
      IMPLICIT NONE
      INTEGER i, j, k, l
      DOUBLE PRECISION, INTENT(IN)	 :: aa33(3,3), bb33(3,3)
      DOUBLE PRECISION, INTENT(OUT)	 :: res4(3,3,3,3)
  
      DO i=1,3
        DO j=1,3
          DO k= 1,3
            DO l=1,3
	   res4(i,j,k,l)=0.5d0*(aa33(i,k)*bb33(j,l)+aa33(i,l)*bb33(j,k))
	    END DO
	  END DO
        END DO
      END DO
   
      RETURN
      END SUBROUTINE
!-------------------------------------------------------------
!-------------------------------------------------------------
      subroutine product_T1dyadT1_nv(a6,b6, res66)
      implicit none

      integer i,j,k,l
      real*8  a6(6),b6(6)!input
      real*8  res66(6,6)!output


      do i=1,6
        do j=1,6
           res66(i,j)=a6(i)*b6(j)
        end do
      end do

      return
      end
!-------------------------------------------------------------
!-------------------------------------------------------------
      subroutine product_T1dyadT1(a3,b3, res33)
      implicit none
! c
      integer i,j,k,l
      real*8  a3(3),b3(3)!input
      real*8  res33(3,3)!output
! c
! c
      do i=1,3
        do j=1,3
           res33(i,j)=a3(i)*b3(j)
        end do
      end do
! c
      return
      end
!-------------------------------------------------------------
!-------------------------------------------------------------
      function macauley(a) result(b)
      implicit none 
      double precision, intent(in):: a
      double precision  b

      b = 0.5d0*(a + abs(a))
      return
      end




! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc]
! c
      subroutine transfer_6to33_nv(a6, b33)
      implicit none
! c
! c.... transfer 6-dim. Voigt notation to 3x3 matrix
! c
      integer i
      real*8  b33(3,3) !output
      real*8  a6(6) !input
      double precision, parameter::w2=sqrt(2.d0)
! c
      do i=1,3
        b33(i,i)=a6(i)
      end do
      b33(1,2)=a6(4)/w2
      b33(2,1)=a6(4)/w2
      b33(2,3)=a6(5)/w2
      b33(3,2)=a6(5)/w2
      b33(3,1)=a6(6)/w2
      b33(1,3)=a6(6)/w2
! c
      return
      end

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc]
! c                                                                        
      subroutine transfer_33to6_nv(b33, a6)
      implicit none
! c
! c.... transfer 3x3 matrix to 6-dim. Voigt notation
! c
      integer i
      double precision ::  a6(6) !output
      double precision ::  b33(3,3) !input
      double precision, parameter::w2=sqrt(2.d0)
! c
      do i=1,3
        a6(i)=b33(i,i)
      end do
      a6(4)=b33(1,2)*w2
      a6(5)=b33(2,3)*w2
      a6(6)=b33(3,1)*w2
! c
      return
      end
      


subroutine getPAndPPr(DD,DChi,HChi,PP,PPr)
implicit none
double precision,intent(in)::DD
double precision,intent(in)::DChi
double precision,intent(in)::HChi
double precision,intent(out)::PP
double precision,intent(out)::PPr

PP=HChi*(DD-DChi)
PPr=HChi

end subroutine getPAndPPr


subroutine getPiHAndPiHPr(DChi,teta,piH,piHPr)
implicit none
double precision,intent(in)::DChi
double precision,intent(in)::teta
double precision,intent(out)::piH
double precision,intent(out)::piHPr

piH=teta*DChi
piHPr=teta

end subroutine getPiHAndPiHPr






      
subroutine umatIsoGradDamFibre(eps,DChi,h1,hn,sig,dSigDeps,&
    &dSigDDChi,pi,dPiDDChi,gradDChi,Xi,dXiDGradDChi,nHistVars,d)

implicit none

include  'tdata.h'
include  'swvars.h' 
! include 'elauto.h'

integer,intent(in)::nHistVars
double precision,intent(in):: eps(6), DChi, hn(nHistVars),gradDChi(3), d(*)
double precision,intent(out)::h1(nHistVars)
double precision,intent(inout)::dSigDeps(6,6), dSigDDChi(6),dPiDDChi,Xi(3),dXiDGradDChi(3,3), sig(6),pi

double precision::dDDChi, dDDEpsf, Nf(3), dblDmy, epsf, Mf(6), DEff

double precision::EModf			!EModf=130.0d3 !130000.0d0
double precision::eta			!eta=10.0d0
double precision::HChi			!HChi=3.0d4
double precision::Yc0			!Yc0=20.d0
double precision::l			!l=0.25d0
double precision, parameter::tol=1.d-9
double precision::fibAng		!fibreAngle=00.d0!in degree
!double precision::HH

double precision::qh, qhPr, Y0, ff
double precision::DD,Dn
!---Neu---
double precision::D_crit,K1

integer::i,j,inelStep
integer, parameter::outp=0

double precision, parameter :: sq2=dsqrt(2.d0)

! Nf=(/0.5d0*sq2, 0.5d0*sq2, 0.d0/)
! Nf=(/0.d0, 1.d0, 0.d0/)
!-- Materialinput --

EModf	=d(13)
fibAng	=d(14)
Yc0	=d(16)
K1	=d(18)
HChi	=d(20)
eta	=d(22)
l	=d(24)
D_crit	=d(32)

!-- End Material input --

dblDmy=fibAng/360.d0*2.d0*3.14159265359d0
Nf=(/cos(dblDmy), sin(dblDmy) , 0.d0/)
Mf(1)=Nf(1)*Nf(1); Mf(2)=Nf(2)*Nf(2); Mf(3)=Nf(3)*Nf(3) 
Mf(4)=sq2*Nf(2)*Nf(3); Mf(5)=sq2*Nf(1)*Nf(3); Mf(6)=sq2*Nf(1)*Nf(2)


!!!Xi
dblDmy=0.d0
do i=1,3
  dblDmy=dblDmy+Nf(i)*gradDChi(i)
  do j=1,3
    dXiDGradDChi(i,j)=dXiDGradDChi(i,j)+EModf*(l**2)*Nf(i)*Nf(j)
  end do
end do
Xi=Xi+EModf*(l**2)*dblDmy*Nf

epsf=0.d0; do i=1,6; epsf=epsf+eps(i)*Mf(i); end do

Y0=.5d0*EModf*epsf*epsf

Dn=hn(18); DD=Dn

call getQhAndQhPr(DD,qh,qhPr,K1,D_crit)
ff=Y0+qh+HChi*(DChi-DD)-Yc0

inelStep=1
if(ff.le.0.d0) inelStep=0
if(duVnshSW.eq.1) inelStep=int(hn(19)+1.d-6)
if(inelStep.eq.0.or.epsf.lt.0.d0) then
  dDDEpsf=0.d0; dDDChi=0.d0; dSigDDChi=dSigDDChi
  DEff=0.d0
  if(epsf.gt.0.d0) DEff=Dn
  h1(19)=0.d0
else
  if(outp.eq.1) write(*,*) 'inel. fibre step, DD=', DD, 'Y0=',Y0, 'epsf=', epsf, 'eps=', eps
  h1(19)=1.d0
  i=1
  do
    i=i+1; if(i.gt.100) then; if(outp.eq.1) write(*,*) 'too many fib. it.'; exit; endif
    call getQhAndQhPr(DD,qh,qhPr,K1,D_crit)
    ff=Y0+qh+HChi*(DChi-DD)-Yc0
    if(dabs(ff).lt.tol) exit
    DD=DD+ff/(HChi-qhPr)
  end do
  dblDmy=(HChi-qhPr)
  dDDEpsf=EModf*epsf/dblDmy; dDDChi=HChi/dblDmy; dSigDDChi=dSigDDChi-dDDChi*EModf*epsf*Mf; DEff=DD
endif

dblDmy=0.d0; do i=1,6; dblDmy=dblDmy+(eps(i)-hn(10+i))*Mf(i); end do
sig=sig + ( (1.d0-DEff)*EModf*epsf + eta/dt*dblDmy ) * Mf
dblDmy=EModf*(1.d0-DEff-epsf*dDDEpsf)
do i=1,6
!   dSigDeps(i,i)=dSigDeps(i,i)+eta/dt
  do j=1,6
    dSigDeps(i,j)=dSigDeps(i,j)+(dblDmy+eta/dt)*Mf(i)*Mf(j)
  end do
end do
! dSigDDChi=dSigDDChi-HChi*dDDEpsf*Mf
pi=pi+HChi*(DChi-DD)
dPiDDChi=dPiDDChi+HChi*(1.d0-dDDChi)
h1(18)=DD
end subroutine umatIsoGradDamFibre



subroutine getQhAndQhPr(DD,qh,qhPr,H2,DHat)
implicit none
double precision,intent(in)::DD, H2, DHat
double precision,intent(out)::qh, qhPr

double precision :: epsi, D0, dblDmy,m0
double precision, parameter:: kappa=1.d0
!double precision, parameter:: DHat=0.9995d0
!double precision, parameter:: HH=0.12d0
double precision, parameter:: m0Aim=1.d8

epsi=(H2*DHat*kappa/m0Aim)**(1.d0/(1.d0+kappa))
D0=DHat-epsi
dblDmy=(DHat-DD)**kappa
if(DD.lt.D0) then !1.eq.1.or.
  qh=-H2*DD/dblDmy; qhPr=-H2*(DHat-DD+kappa*DD)/dblDmy/(DHat-DD)
else
  m0=H2*(DHat-D0+kappa*D0)/(DHat-D0)**(kappa+1.d0)
  qh=-H2*D0/( (DHat-D0)**kappa ) - m0*(DD-D0); qhPr=-m0         ! Aenderung -m0 zu +m0
endif

end subroutine getQhAndQhPr
! !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc]
!       subroutine numTangent01(xres,xres_n,R0,tr_eps_star,eps_star_pos_square,eps_star_neg_square,Dchi,ntan,D0)
! !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc]
!       use sharedConstants01
!       implicit none
! 
!       integer i,j,k,nres
! 
!       double precision  xres(2),xres_n(2),R0(2),RV(2)
!       double precision  ntan(2,2),DGP,DGP_n,lambdabar,lambdabar_n
!       double precision  xv(2),delta,D_acc,YY,macauley,qq,qqpr,Hp
!       double precision  tr_eps_star,eps_star_pos_square
!       double precision  eps_star_neg_square,Dchi
! 
! !     Initialize
!       do i=1,2
!        do j=1,2
!         ntan(i,j)=0.d0
!        enddo
!       enddo
!       
! 
!       delta=1.d-8
!       do i = 1,2
!          xv(i) = xres(i)
!       end do
!       
!       if (h_switch .eq. 4) then
!        Hp=0
!       else
!        Hp=Hpen
!       endif
! 
!       DGP_n=xres_n(1)
!       lambdabar_n=xres_n(2)     
! 
!       do i = 1,2
!         xv(i) = xres(i) + delta
!         DGP=xv(1)
!         lambdabar=xv(2)
!         D_acc=DGP_n + lambdabar
! 
! !    computation of the damage driving force Y
!       YY=0.5d0*ndam*(1.d0-DGP)**(ndam-1.d0)*lambda*(macauley(tr_eps_star))**2.d0 &
!        +0.5d0*ndam*(1.d0-(1.d0-eta_cc)*DGP)**(ndam-1.d0)*lambda*(1.d0-eta_cc)*(macauley(-tr_eps_star))**2.d0&
!        +mu*ndam*(1.d0-DGP)**(ndam-1.d0)*eps_star_pos_square +mu*ndam*(1.d0-(1.d0-eta_cc)*DGP)**&
!        (ndam-1.d0)*(1.d0-eta_cc)*eps_star_neg_square - Hp*macauley(DGP-D0)+ kappa + cP1*(DChi-DGP)
! 
! !     computation of the hardening variable qq
!       if (h_switch .eq. 1) then
!        qq=k1*D_acc
!       elseif (h_switch .eq. 2) then
!        qq=k1*D_acc+k2*D_acc**2.d0
!       elseif (h_switch .eq. 3) then
!        qq=k0*(1.d0-dexp(-alpha*D_acc))
!       elseif (h_switch .eq. 4) then
!        call getQhAndQhPr(D_acc,qq,qqpr,k1,D0)
!        qq=-qq
!        qqpr=-qqpr
!       end if
! 
! 
!         call compute_residuum01(DGP,DGP_n,lambdabar_n,lambdabar,YY,qq,RV)
!  
!          do j = 1,2
!             ntan(j,i) = (RV(j) - R0(j))/delta
!          end do
!          xv(i) = xv(i) - delta        
!       end do
!       return
!       end        
!       
!       
      
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc]      