! module testMod
! implicit none
! integer testIntgr
! end module testMod

subroutine elmt01(      &
                   d,   &    ! element data parameters
                   ul,  &    ! element nodal solution parameters
                   xl,  &    ! element nodal reference coordinates
                   ix,  &    ! element global node numbers
                   tl,  &    ! element nodal temperature values
                   s,   &    ! element matrix (e.g. stiffness, mass)
                   p,   &    ! element vector (e.g. residual, mass)
                   ndf, &    ! number of unknowns (max per node)
                   ndm, &    ! space dimension of mesh 
                   nst, &    ! size of element arrays s and p
                   isw  &    ! task parameter to control computation
                 )
! use testMod
implicit none
                                
include 'cdata.h'
include 'tdata.h'
include 'pointer.h'
include 'hdata.h'
include 'comblk.h'
include 'swvars.h'
include 'eldata.h'
 
!---------------- subroutine I/O parameter ----------------!
integer, intent(in) :: ndf, ndm, nst, isw
integer :: ix(*)
double precision, intent(in) :: ul(ndf,nen,*), xl(ndm,nen), tl(*)
double precision, intent(out) :: s(nst,nst), p(nst) !TN:p->r for residual!SW2:Changed it back. If there are names by feap for quantities we should use them. This prevents confusions
double precision :: kUU(ndm*nen,ndm*nen), fU(ndm*nen)

!---------------- Uebergabeparameter SHP----------------! 
! double precision :: SHP(ndm+1,nen),XJAC

!---------------- Zwischenwerte ----------------!
integer :: nTens, frstStp, calcX, Nintp, pos1, pos2
integer :: i,j,k,l,m,o,ip            !zum Durchzählen 
double precision :: xlc(ndm,nen), tempFF(ndm,ndm), FF(ndm,ndm), tau_dxd(ndm,ndm)
double precision :: gg(nen,nen), Btau(ndm*nen), KMat(ndm*nen,ndm*nen)!SW2: recalled 'stiffTens' to KMat: it's the material tangent matrix (not a tensor)
double precision, allocatable :: SW(:,:), tau_Voigt(:), ca(:,:) 

double precision, allocatable::SHP(:,:,:), XJAC(:), omega(:)

!--------for P0-extension:----------------------
double precision :: theta, pp, omegaEl, U2Pr,dblDmy


!---------------- Konstanten ----------------!
double precision, parameter :: w2=dsqrt(2.d0)
integer, parameter::outpt=0! if 1, print comments on screen
integer, parameter::elType=1!usually 1
integer, parameter::p0On=0
integer, parameter::nHistVars=1+16+2


!++++++++++++++++++++++++++++++++++gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!if one comments everything in between 'gradient-extension' 'end gradient-extension' out, one gets the standard element
integer, parameter::gradExtOn=1!0 or 1. 1 if gradient extended, otherwise 0->standard element without gradient extension
integer, parameter::nPhi=2!dimension of phi
integer :: ndfPhi
double precision, allocatable :: B0(:,:), pi(:), xi(:), dTauDPhi(:,:), cDPi(:,:), dPiDPhi(:,:), dXiDGPhi(:,:)
double precision, allocatable :: kPhiU(:,:), kUPhi(:,:), kPhiPhi(:,:), fPhi(:)
double precision, allocatable :: dblMat(:,:), dblMat2(:,:), phi(:),gradPhiVoigt(:)!, phiHat(:)
!++++++++++++++++++++++++++++++end gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  

!++++++++++++++ User Input ++++++++++++++++++
logical :: pcomp
double precision:: d(*), param(19)
character :: text(2)*15
logical :: inputs
logical :: errck
logical :: tinput
!++++++++++++ end User Input ++++++++++++++++

! write(*,*) '----------------------------------'
! write(*,*) 'ISW==', isw
! write(*,*) '----------------------------------'
nTens=3
if(ndm.eq.3) nTens=6
! testIntgr=testIntgr+1
call getNp(elType,ndm,nen,Nintp)
!write(*,*) 'ndf=', ndf, 'ndm=', ndm, 'nst=', nst, 'nen=', nen, 'nHistVars=', nHistVars, 'Nintp=', Nintp

if (isw.eq.1) then
!   testIntgr=-10
  nh1=nHistVars*Nintp 
  
  inputs = .TRUE.

  do while (inputs)
  errck = .TRUE.
    do while (errck)
      errck = tinput(text,2,param,15)
    end do
  
  if ( pcomp(text(1),'dama',4) ) then

      if ( pcomp(text(2),'mate',4) ) then ! material parameters
      
	d(11) = param(1)   ! Youngs modulus
        d(12) = param(2)   ! Poissons ratio
        d(13) = param(3)   ! Fiber stiffness
        d(14) = param(4)   ! Fiber angle
        d(15) = param(5)   ! Damage Surface matrix
        d(16) = param(6)   ! Damage Surface fiber
        d(17) = param(7)   ! Hardening parameter/ initial slope H_matrix
        d(18) = param(8)   ! Hardening parameter/ initial slope H_fiber
        d(19) = param(9)   ! energy penalizing parameter HChi matrix
        d(20) = param(10)  ! energy penalizing parameter HChi fiber
        d(21) = param(11)  ! eta matrix
        d(22) = param(12)  ! eta fiber
        d(23) = param(13)  ! char. length matrix
        d(24) = param(14)  ! char. length fiber
        !d(27) = param(15)  ! energy penalizing parameter MEnPen2*(DChi_2-D_2)
        !d(28) = param(16)  ! Slope after D_crit m0
        !d(29) = param(17)  ! pseudoviscosity eta
        
        
      elseif ( pcomp(text(2),'comp',4) ) then ! computation parameters
	
	d(30) = param(1)	! CC? 1=yes 0=no
	d(31) = param(2)	! CC parameter 1=full 0=none
	d(32) = param(3)	! D_crit
	d(33) = param(4)	! hardening switch (1=lin, 2=quad, 3=exp, 4=qq->inf)
	d(34) = param(5)	! switch (1=strain equivalence, 2=energy equivalence) 
	d(35) = param(6)	! factor alpha_1 for fiber dependence of xi
	d(36) = param(7)	! factor alpha_2 for fiber dependence of xi
	
	!write(*,*) d(30), d(31)
      
      endif
  
  elseif ( pcomp(text(1),'    ',4) ) then    
      inputs = .FALSE.
    else
      write(*,*) 'ERROR: Unknown identifier found! STOP!'
      STOP
    endif      
  end do
end if

if(isw.eq.3 .or. isw.eq.6) then

  allocate(SHP(Nintp,ndm+1,nen)); allocate(XJAC(Nintp)); allocate(omega(Nintp))

  allocate(tau_Voigt(nTens))
  allocate(ca(nTens,nTens))
!++++++++++++++++++++++++++++++++++gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  if(gradExtOn.eq.1) then
    allocate(B0(ndm,nen)); allocate(pi(nPhi)); allocate(xi(ndm*nPhi))
    allocate(dTauDPhi(nTens,nPhi)); allocate(cDPi(nPhi,nTens))
    allocate(dPiDPhi(nPhi,nPhi)); allocate(dXiDGPhi(ndm*nPhi,ndm*nPhi)); allocate(kPhiU(nPhi*nen,ndm*nen))
    allocate(kUPhi(ndm*nen,nPhi*nen)); allocate(kPhiPhi(nPhi*nen,nPhi*nen)); allocate(fPhi(nPhi*nen))
    i=nPhi
    if(nPhi.lt.ndm) i=ndm
    allocate(dblMat(nen*i,nen*i)); allocate(dblMat2(nen*i,nen*i)); 
    allocate(gradPhiVoigt(ndm*nPhi)); allocate(phi(nPhi))!; allocate(phiHat(nPhi*nen))
  endif
!++++++++++++++++++++++++++++++end gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  

  frstStp = 0
  calcX = 0

!-----------NEU-LP----------------------
  duVnshSW=1
  do I=1,nen; do j=1,2
    if(dabs(ul(j,I,2)).gt.1.d-10) duVnshSW=0
  end do; end do
!-----------NEU-ENDE-LP----------------------


  if (dabs(hr(nh1)).lt.1.d-16) then; frstStp=1; if(outpt.eq.1) write(*,*) 'frstStp=1'; endif

  if (frstStp==1.or.p0On==1) calcX=1

  xlc = xl + ul(1:ndm,:,1)   !xlc in Momentankonfig. !TN:(1:ndm,1:nen,1) not necessary
!   write(*,*) 'ul='
!   do i=1,ndm; write(*,*) ul(i,:,1); end do
!   write(*,*) 'xlc='
!   do i=1,ndm; write(*,*) xlc(i,:); end do

  allocate(SW(ndm+1,Nintp))
  call getIntPtData(nen,Nintp,ndm,elType,SW)

  kUU = 0.d0; fU = 0.d0
!++++++++++++++++++++++++++++++++++gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  if(gradExtOn.eq.1) then
!     do I=1,nen; phiHat(nPhi*(I-1)+1:nPhi*I)=ul(ndm+1:ndf,I,1); end do
    kUPhi=0.d0; kPhiU=0.d0; kPhiPhi=0.d0; fPhi=0.d0
  endif
!++++++++++++++++++++++++++++++end gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  

  do ip=1,Nintp
    call getSHP(xlc,ndm,nen,elType,calcX,XJAC(ip),SHP(ip,:,:),SW(1:ndm,ip))
    
    pos1=nh1+(ip-1)*nHistVars 
    pos2=nh2+(ip-1)*nHistVars

    if (frstStp==1) hr(pos1)=XJAC(ip)*SW(ndm+1,ip)!SW: Remember: NEVER STORE ANYTHING IN THE HISTORY VARIABLES AT THE BEGINNING OF THE TIME STEP (='pos1') EXCEPT FROM THIS ONE SPECIAL CASE!
    omega(ip)=hr(pos1); hr(pos2)=omega(ip)
!     write(*,*) 'omega=', omega(ip), 'ip=', ip
!     write(*,*) ' '
!     write(*,*) 'check if omega is properly stored in histVars. After that delete this write command., omega=', omega
    if (frstStp.eq.1.and.dabs(sum(ul(1,:,1))).gt.1.d-12) then
      write(*,*) 'This should not happen. If it does happen, we have to think again about the code.'
      stop
    endif
  end do

!----------------P0-Extension---------------------------------------------------------------------------------  
  if(p0On.eq.1) then
    omegaEl=0.d0; theta=0.d0
    dblMat(1,:)=0.d0
    do ip=1,Nintp
      dblDmy=XJAC(ip)*SW(ndm+1,ip)
      theta=theta+dblDmy
      omegaEl=omegaEl+omega(ip)
      do I=1,nen
        k=ndm*(I-1)+1; l=ndm*I
        dblMat(1,k:l)=dblMat(1,k:l)+SHP(ip,:,I)*dblDmy
      end do
    end do
    theta=theta/omegaEl
    call getUPrAndU2Pr(theta, pp, U2Pr)
    dblDmy=U2Pr/omegaEl
    do i=1,ndm*nen; do j=1,ndm*nen
      kUU(i,j)=dblDmy*dblMat(1,i)*dblMat(1,j)
    end do; end do
  endif
!-------------------------------------------------------------------------------------------------------------
    
  do ip=1,Nintp    
    pos1=nh1+(ip-1)*nHistVars 
    pos2=nh2+(ip-1)*nHistVars
    tempFF=matmul(xl,transpose(SHP(ip,1:ndm,:)))
!     FF=tempFF
!     write(*,*) 'FF='; do i=1,ndm; write(*,*) FF(i,:); end do
!     write(*,*) 'SHP='; do i=1,ndm+1; write(*,*) SHP(ip,i,:); end do
!     write(*,*) 'invert durch Lapack routine ersetzen'
!     stop
!     call invert(FF,ndm,ndm)
!     call schreibeb(tempFF,3,3,'IVGr')
    call CalcInverse(ndm,tempFF,FF)
!     call schreibeb(FF,3,3,'Grad')
!++++++++++++++++++++++++++++++++++gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    if(gradExtOn.eq.1) then
      B0=matmul(transpose(FF),SHP(ip,1:ndm,:))
!       subroutine getBPhiTMat(oo,BB,vv,nen,nPhi,ndm,m)
!       call getBPhiTMat(gradPhiVoigt,B0,phiHat,nen,nPhi,ndm,1)
      phi=0.d0      
      do i=1,nPhi
        do j=1,nen
          phi(i)=phi(i)+SHP(ip,ndm+1,j)*ul(ndm+i,j,1)
        end do
        gradPhiVoigt(ndm*(i-1)+1:ndm*i)=matmul(B0,ul(ndm+i,:,1))
      end do
    endif
!++++++++++++++++++++++++++++++end gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    if(outpt.eq.1) then; write(*,*) 'FF='; do j=1,ndm; write(*,*) FF(j,:); end do; endif
!     call CalcInverse(ndm,tempFF,FF)
    
    if(gradExtOn.eq.0) then 
      call UmatInterface(d,FF,ndm,nTens,tau_Voigt,ca,hr(pos1:pos1+nHistVars-1),hr(pos2:pos2+nHistVars-1),nHistVars)
    else
!++++++++++++++++++++++++++++++++++gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
      call UmatInterfaceGradExt(d,FF,phi,gradPhiVoigt,ndm,nTens,nPhi,tau_Voigt,pi,xi,ca,cDPi,dTauDPhi,dPiDPhi,&
          &dXiDGPhi,hr(pos1:pos1+nHistVars-1),hr(pos2:pos2+nHistVars-1),nHistVars, isw)
!++++++++++++++++++++++++++++++end gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    endif
    
!----------------P0-Extension---------------------------------------------------------------------------------  
    if(p0On.eq.1) then
      dblDmy=pp*XJAC(ip)*SW(ndm+1,ip)/omega(ip)!p*J
      do i=1,ndm
        tau_Voigt(i)=tau_Voigt(i)+dblDmy
        do j=1,ndm; ca(i,j)=ca(i,j)+dblDmy; end do
      end do
      dblDmy=2.d0*dblDmy
      j=3; if(ndm.eq.3) j=6
      do i=1,j; ca(i,i)=ca(i,i)-dblDmy; end do
    endif
!-------------------------------------------------------------------------------------------------------------
    
    call trnsfrmVec2Mat(tau_dxd,tau_Voigt,ndm,nTens)
    
    gg = matmul(matmul(transpose(SHP(ip,1:ndm,:)),tau_dxd),SHP(ip,1:ndm,:))
    
    do o = 0,nen-1
      do j = 0,nen-1
	do k = 1,ndm
	l = ndm*o + k
	m = ndm*j + k
	kUU(l,m) = kUU(l,m) + gg(o+1,j+1)*omega(ip)!SW3: omega was missing
	end do !k
      end do !j
    end do !n
    
    call generate_K(KMat,SHP(ip,1:ndm,1:nen),ca,nen,nTens,ndm) !then we don't have to allocate and deallocate more variables
    kUU = kUU + KMat*omega(ip)
!     write(*,*) 'tau_Voigt=', tau_Voigt
    call getBEpsTMat(Btau,SHP(ip,1:ndm,1:nen),tau_Voigt,nen,nTens,ndm,1)
    fU = fU - Btau*omega(ip)
    
!++++++++++++++++++++++++++++++++++gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    if(gradExtOn.eq.1) then
!       subroutine getNPhiTMat(oo,NN,vv,nen,nPhi,m)
      call getNPhiTMat(dblMat(1:nPhi*nen,1:nTens),SHP(ip,ndm+1,1:nen),transpose(dTauDPhi),nen,nPhi,nTens)
      call getBEpsTMat(dblMat2(1:ndm*nen,1:nPhi*nen),SHP(ip,1:ndm,1:nen),transpose(dblMat(1:nPhi*nen,1:nTens)),&
         &nen,nTens,ndm,nPhi*nen)
      kUPhi=kUPhi+omega(ip)*dblMat2(1:ndm*nen,1:nPhi*nen)

      call getNPhiTMat(dblMat(1:nPhi*nen,1:nTens),SHP(ip,ndm+1,1:nen),cDPi,nen,nPhi,nTens)
      call getBEpsTMat(dblMat2(1:ndm*nen,1:nPhi*nen),SHP(ip,1:ndm,1:nen),transpose(dblMat(1:nPhi*nen,1:nTens)),&
            &nen,nTens,ndm,nPhi*nen)
      kPhiU=kPhiU+omega(ip)*transpose(dblMat2(1:ndm*nen,1:nPhi*nen))
      
      call getNPhiTMat(dblMat(1:nPhi*nen,1:nPhi),SHP(ip,ndm+1,1:nen),transpose(dPiDPhi),nen,nPhi,nPhi)
      call getNPhiTMat(dblMat2(1:nPhi*nen,1:nPhi*nen),SHP(ip,ndm+1,1:nen),transpose(dblMat(1:nPhi*nen,1:nPhi)),nen,nPhi,nPhi*nen)
      kPhiPhi=kPhiPhi+omega(ip)*dblMat2(1:nPhi*nen,1:nPhi*nen)
      
!       subroutine getBPhiTMat(oo,BB,vv,nen,nPhi,ndm,m)
      call getBPhiTMat(dblMat(1:nPhi*nen,1:ndm*nPhi),B0,transpose(dXiDGPhi),nen,nPhi,ndm,ndm*nPhi)
      call getBPhiTMat(dblMat2(1:nPhi*nen,1:nPhi*nen),B0,transpose(dblMat(1:nPhi*nen,1:ndm*nPhi)),nen,nPhi,ndm,nPhi*nen)
      kPhiPhi=kPhiPhi+omega(ip)*dblMat2(1:nPhi*nen,1:nPhi*nen)
      
      call getBPhiTMat(dblMat(1:nPhi*nen,1),B0,xi,nen,nPhi,ndm,1)
      call getNPhiTMat(dblMat2(1:nPhi*nen,1),SHP(ip,ndm+1,1:nen),pi,nen,nPhi,1)
      fPhi=fPhi-omega(ip) * ( dblMat(1:nPhi*nen,1) + dblMat2(1:nPhi*nen,1) )
!       kPhiU=0.d0; kPhiPhi=0.d0
    endif
!++++++++++++++++++++++++++++++end gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  end do !end Gausspunkt
  
  if(gradExtOn.eq.0) then
    s=kUU; p=fU
  else
!++++++++++++++++++++++++++++++++++gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    ndfPhi=ndf-ndm
    p = 0.d0                               
    do i=1,nen
      do j=1,ndm
        p(ndf*(i-1)+j) = fU(ndm*(i-1)+j)
      end do
      do k=1,ndfPhi
        p(ndf*(i-1)+ndm+k) = fPhi(ndfPhi*(i-1)+k)
      end do
    end do
  
    do i=0,nen-1                                                                                    
      do j=0,nen-1
        do k=1,ndm; do l=1,ndm;
          s(ndf*i+k, ndf*j+l)=kUU(ndm*i+k,ndm*j+l)
        end do; end do
        do k=1,ndm; do l=1,ndfPhi
          s(ndf*i+k, ndf*j+ndm+l)=kUPhi(ndm*i+k,ndfPhi*j+l)
        end do; end do
        do k=1,ndfPhi; do l=1,ndm
          s(ndf*i+ndm+k, ndf*j+l)=kPhiU(ndfPhi*i+k,ndm*j+l)
        end do; end do
        do k=1,ndfPhi; do l=1,ndfPhi
          s(ndf*i+ndm+k, ndf*j+ndm+l)=kPhiPhi(ndfPhi*i+k,ndfPhi*j+l)
        end do; end do
      end do    
    end do
!++++++++++++++++++++++++++++++end gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++      
  endif
!   write(*,*) 'fU=', fU

  deallocate(SHP); deallocate(XJAC); deallocate(omega)
  
  deallocate(SW)
  deallocate(ca)
  deallocate(tau_Voigt)

!++++++++++++++++++++++++++++++++++gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  if(gradExtOn.eq.1) then
    deallocate(B0); deallocate(pi); deallocate(xi); deallocate(dTauDPhi); deallocate(cDPi)
    deallocate(dPiDPhi); deallocate(dXiDGPhi); deallocate(kPhiU)
    deallocate(kUPhi); deallocate(kPhiPhi); deallocate(fPhi)
    deallocate(dblMat); deallocate(dblMat2); deallocate(gradPhiVoigt); deallocate(phi)!; deallocate(phiHat)
  endif
!++++++++++++++++++++++++++++++end gradient-extension+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  
end if !end isw

end subroutine elmt01
      
      

      
      
      
      
subroutine getUPrAndU2Pr(theta, pp, U2Pr)
implicit none

double precision, intent(in) :: theta
double precision, intent(out) :: pp, U2Pr
double precision, parameter :: kappa0 = 1.d4

pp=kappa0*(theta-1.d0); U2Pr=kappa0

end subroutine getUPrAndU2Pr
      
      
      
      
      
subroutine getNp(elType, ndm, nen,Np)
implicit none 

integer, intent(in):: elType, nen, ndm
integer, intent(out) :: Np

if (elType==1) then 

  if (ndm==2) then
  
    if (nen==3) then 
      Np=1     !lin. triangle
    else if (nen==4) then 
      Np=4      !bilin. quad.
    else 
      write(*,*) 'ERROR.UnknownElementType2D.Np'
      stop
    end if  
    
  else if (ndm==3) then
    
    if (nen==4) then 
      Np=1     !lineare Tetr.
    else if (nen==10) then 
      Np=5 !quadratische Tetraeder
    else if (nen==8) then 
      Np=8      !trilin. hexa.
    else 
      write(*,*) 'ERROR.ERROR.UnknownElementType3D.Np'
      stop
    end if
    
  else 
    write(*,*) 'ERROR.ERROR.UnknownElementTypeWeder2Dnoch3D.Np'
    stop
  end if  
  
else 
  write(*,*) 'ERROR.ERROR.UnknownElement.Np'  
  stop
end if
end subroutine getNp





subroutine getIntPtData(nen,Np,ndm,elType,SW)
implicit none

integer, intent(in) :: nen,ndm, elType, Np 
integer :: i
double precision, intent(out) :: SW(ndm+1,Np)
double precision :: aa,bb,cc,ee,ff,gg

if (elType==1) then 

  if (ndm==2) then   
  
    if (nen==3) then  !Np=1
      SW(1:3,1)=(/ 1.d0/3.d0, 1.d0/3.d0, .5d0 /) 
    else if (nen==4) then  !Np=4 
!       write(*,*) 'bilin quad.'
      aa=1.d0/dsqrt(3.d0)
      SW(1,1:4)= (/ -aa, aa, aa, -aa /)    
      SW(2,1:4)=(/ -aa, -aa, aa, aa /)    
      do i=1,nen
        SW(3,i)=1.d0
      end do
    else  
      write(*,*) 'ERROR.keine.Gausspkt.fuer2DElemente'
      stop
    end if
    
  else if (ndm==3) then 
  
    if (nen==4) then !Np=1
      SW(1:4,1)=(/ .25d0, .25d0, .25d0, 1.d0/6.d0 /)
    else if (nen==8) then !Np=8
      aa=1.0d0/dsqrt(3.d0)
      SW(1,1:8)=(/ -aa, aa, aa, -aa, -aa, aa, aa, -aa /)
      SW(2,1:8)= (/ -aa, -aa, aa, aa, -aa, -aa, aa, aa /)
      SW(3,1:8)=(/ -aa, -aa, -aa, -aa, aa, aa, aa, aa/)
      do i=1,nen
        SW(4,i)=1.d0
      end do
    
    else if (nen==10) then !Np=5
      bb=.5d0
      cc=.25d0
      gg=1.d0/6.d0
      ee=3.d0/40.d0
      ff=2.d0/15.d0
      
      SW(1,1:5)=(/ cc, gg, gg, gg, bb/)
      SW(2,1:5)= (/ cc, gg, gg, bb, gg /)
      SW(3,1:5)=(/ cc, gg, bb, gg, gg /)
      SW(4,1:5)= (/ -ff, ee, ee, ee, ee /)
  
!       !TN3: Np=4 Alternative!SW: didn't check that alternative!!! NIntp must also be set to 4 in 'getNp', if this alternative is used!!!!!!!
!       aa=.25d0-.05d0*dsqrt(5.d0)
!       bb=.25d0+.15d0*dsqrt(5.d0)
!       cc=1.d0/24.d0
!       
!       SW(1,1:4)=(/ bb, aa, aa, aa /)
!       SW(2,1:4)= (/ aa, bb, aa, aa /)
!       SW(3,1:4)=(/ aa, aa, bb, aa /)
!       SW(4,1:4)= (/ cc, cc, cc, cc /)    
    else  
      write(*,*) 'ERROR.keine.Gausspkt.fuer3DElemente' 
      stop
    end if
   
  
  else  
    write(*,*) 'ERROR.keine.Ahnung.was.fuer.Dimension'  
    stop
  end if
  
else 
  write(*,*) 'ERROR.keine.Ahnung.was.fuer.ein.elType'   
  stop
end if

end subroutine getIntPtData


       
subroutine trnsfrmVec2Mat(AMat,AVec,ndm,nTens)

  implicit none
  double precision :: AMat(ndm,ndm), AVec(nTens), sqrt22 !TN:changed it into Fortran95
  integer :: i,j,k,ndm,nTens !TN: "::" was missing
  sqrt22=.5d0*dsqrt(2.d0)
  if(ndm.eq.2) then
    do i=1,2; AMat(i,i)=AVec(i); end do
    AMat(1,2)=sqrt22*AVec(3); AMat(2,1)=AMat(1,2)
  else
    do i=1,3; AMat(i,i)=AVec(i); end do
    AMat(2,3)=sqrt22*AVec(4); AMat(1,3)=sqrt22*AVec(5)
    AMat(1,2)=sqrt22*AVec(6);
    AMat(3,2)=AMat(2,3); AMat(3,1)=AMat(1,3); AMat(2,1)=AMat(1,2)
  endif

end subroutine



subroutine trnsfrmMat2Vec(AMat, AVec,ndm,nTens)

  implicit none
  double precision ::  AMat(ndm,*), AVec(nTens), sqrt2
  integer ::  i,j,k,ndm,nTens !TN:same here
  sqrt2=dsqrt(2.d0)
  if(ndm.eq.2) then
    do i=1,2; AVec(i)=AMat(i,i); end do
    AVec(3)=sqrt2*AMat(1,2)
  else
    do i=1,3; AVec(i)=AMat(i,i); end do
    AVec(4)=sqrt2*AMat(2,3); AVec(5)=sqrt2*AMat(1,3)
    AVec(6)=sqrt2*AMat(1,2)
  endif

end subroutine


subroutine getSHP(xlc,ndm,nen,elType,calcX,XJAC,SHP,xi) 
implicit none

integer, intent(in) :: ndm,nen,elType,calcX
double precision, intent(in) :: xlc(ndm,nen),xi(ndm)
double precision, intent(out) :: XJAC,SHP(ndm+1, nen)
double precision :: invj(ndm,ndm), SHPXi(ndm+1,nen),jj(ndm,ndm)                    !lokales Zeug
integer::i

call getSHPXi(ndm,nen, elType, xi, SHPXi)
! write(*,*) 'SHPXi='
! do i=1,ndm+1; write(*,*) SHPXi(i,:); end do
invj=matmul(xlc,transpose(SHPXi(1:ndm,:)))!SW: In fact, this is still j here,  but is inverted below, therefore it is already called invj here.
! write(*,*) 'j='
! do i=1,ndm; write(*,*) invj(i,:); end do

if (calcX==1) then
  if(ndm==2) then 
    XJAC=invj(1,1)*invj(2,2)-invj(1,2)*invj(2,1)
  else 
    XJAC=invj(1,1)*(invj(2,2)*invj(3,3)-invj(3,2)*invj(2,3)) &
    +invj(1,2)*(invj(2,3)*invj(3,1)-invj(3,3)*invj(2,1)) &
    +invj(1,3)*(invj(2,1)*invj(3,2)-invj(3,1)*invj(2,2))
  end if
end if

! jj=invj
! call invert(invj,ndm,ndm)
call CalcInverse(ndm,invj,jj)
! jj=matmul(jj,invj)
! write(*,*) 'ndm=', ndm, 'test='; do i=1,ndm; write(*,*) jj(i,:); end do

SHP(1:ndm,:)=matmul(transpose(jj),SHPXi(1:ndm,:))
SHP(ndm+1,:)=SHPXi(ndm+1,:)
!     write(*,*) 'XJAC=', XJAC, 'SHP1='; do i=1,ndm+1; write(*,*) SHP(i,:); end do

end subroutine getSHP






subroutine getSHPXi(ndm,nen, elType, xi, SHPXi)
implicit none

integer, intent(in) :: ndm, nen, elType
double precision, intent(in) :: xi(ndm)
double precision, intent(out) :: SHPXi(ndm+1,nen)
double precision :: aa,bb,cc,dd,ee,ff,gg, lambda,c1_8 !lokales Zeug 




if (elType==1) then 
  
  if (ndm==2) then !2D
  
    if (nen==3) then !Dreieckselement mit linearen Ansätzen
      SHPXi(1,:)=(/ -1.d0, 1.d0, 0.d0 /)
      SHPXi(2,:)=(/ -1.d0, 0.d0, 1.d0 /)
      SHPXi(3,:)=(/ 1.d0-xi(1)-xi(2),xi(1),xi(2) /)
      
    else if (nen==4) then !Viereckselement mit linearen Ansätzen
      aa=0.25d0*(1.d0-xi(1))
      bb=0.25d0*(1.d0+xi(1))
      cc=1.d0-xi(2)
      ee=1.d0+xi(2)
      ff=0.25d0*cc
      gg=0.25d0*ee
      
      SHPXi(1,:)= (/-ff, ff, gg, -gg /)
      SHPXi(2,:)= (/-aa, -bb, bb, aa /)
      SHPXi(3,:)=(/ aa*cc, bb*cc, bb*ee, aa*ee /)!SW2: replaced 'ndm+1' by '3'. Also everywhere below by appropriate number. Faster. 
      
    else 
      write(*,*) 'ERROR.keine.Ahnung.was.fuer.ein.2DElement.SHPXi' 
      stop    
    end if 
  
  else if (ndm==3) then   !3D 
  
    if (nen==4) then !Tetraeder mit linearen Ansätzen
      SHPXi(1,:)=(/ -1.d0, 1.d0, 0.d0, 0.d0 /)
      SHPXi(2,:)=(/ -1.d0, 0.d0, 1.d0, 0.d0 /)
      SHPXi(3,:)=(/ -1.d0, 0.d0, 0.d0, 1.d0 /)
      SHPXi(4,:)= (/ (1.d0-xi(1)-xi(2)-xi(3)), xi(1), xi(2), xi(3)/)
    
    else if (nen==8) then !Hexaeder mit trilinearen Ansätzen
      c1_8=0.125d0 
      aa=c1_8*(1.d0-xi(1)); bb=c1_8*(1.d0+xi(1))
      cc=1.d0-xi(2); dd=1.d0+xi(2); ee=1.d0-xi(3); ff=1.d0+xi(3) 
      
      SHPXi(1,:)=(/ -c1_8*cc*ee, c1_8*cc*ee, c1_8*dd*ee, -c1_8*dd*ee, -c1_8*cc*ff, c1_8*cc*ff, c1_8*dd*ff, -c1_8*dd*ff /)
      SHPXi(2,:)=(/ -aa*ee, -bb*ee, bb*ee, aa*ee, -aa*ff, -bb*ff, bb*ff, aa*ff /)
      SHPXi(3,:)=(/ -aa*cc, -bb*cc, -bb*dd, -aa*dd, aa*cc, bb*cc, bb*dd, aa*dd /)

      SHPXi(4,1)= aa*cc*ee
      SHPXi(4,2)= bb*cc*ee
      SHPXi(4,3)= bb*dd*ee
      SHPXi(4,4)= aa*dd*ee
      SHPXi(4,5)= aa*cc*ff
      SHPXi(4,6)= bb*cc*ff
      SHPXi(4,7)= bb*dd*ff
      SHPXi(4,8)= aa*dd*ff

    else if (nen==10) then !Tetraeder mit quadratischen Ansätzen
      lambda= (1.d0-xi(1)-xi(2)-xi(3))
     
      SHPXi(1,:)=(/ -4.d0*lambda+1.d0, 4.d0*xi(1)-1.d0, 0.d0, 0.d0, 4.d0*(lambda-xi(1)), 4.d0*xi(2), -4.d0*xi(2),&
                &-4.d0*xi(3), 4.d0*xi(3), 0.d0 /) 
      SHPXi(2,:)=(/ -4.d0*lambda+1.d0, 0.d0, 4.d0*xi(2)-1.d0, 0.d0, -4.d0*xi(1), 4.d0*xi(1), 4.d0*(-xi(2)+lambda),&
                &-4.d0*xi(3), 0.d0, 4.d0*xi(3) /) 
      SHPXi(3,:)=(/ -4.d0*lambda+1.d0, 0.d0, 0.d0, 4.d0*xi(3)-1.d0, -4.d0*xi(1), 0.d0, -4.d0*xi(2),&
                &4.d0*(lambda-xi(3)), 4.d0*xi(1), 4.d0*xi(2) /)
   
      SHPXi(4,1)= lambda*(2.d0*lambda-1.d0)
      SHPXi(4,2)= xi(1)*(2.d0*xi(1)-1.d0)
      SHPXi(4,3)= xi(2)*(2.d0*xi(2)-1.d0)
      SHPXi(4,4)= xi(3)*(2.d0*xi(3)-1.d0)
      SHPXi(4,5)= 4.d0*xi(1)*lambda
      SHPXi(4,6)= 4.d0*xi(1)*xi(2)
      SHPXi(4,7)= 4.d0*xi(2)*lambda
      SHPXi(4,8)= 4.d0*xi(3)*lambda
      SHPXi(4,9)= 4.d0*xi(1)*xi(3)
      SHPXi(4,10)= 4.d0*xi(2)*xi(3)

    else 
      write(*,*) 'ERROR.keine.Ahnung.was.fuer.ein.3DElement.mit.wievielen.Knoten.SHPXi'
      stop
    end if
  
  else  
    write(*,*) 'ERROR.keine.Ahnung.was.fuer.eine.Dimension.SHPXi' 
    stop
  end if
  
else  
  write(*,*) 'ERROR.keine.Ahnung.was.fuer.ein.elType.SHPXi' 
  stop
end if 


end subroutine getSHPXi

! getBEpsTMat(Btau,SHP(1:ndm,nen),tau_Voigt,nen,nTens,ndm,1)
! call getBEpsTMat(dblMat(1:ndm*nen,1:nPhi),SHP(1:ndm,1:nen),dTauDPhi,nen,nTens,ndm,nPhi)

subroutine getBEpsTMat(oo,BB,vv,nen,nTens,ndm,m)
implicit none 
!SHP(1:ndm,nen)=BB
integer, intent(in) :: nen,nTens,ndm,m
double precision, intent(in) :: BB(ndm,nen),vv(nTens,m)
double precision, intent(out) :: oo(ndm*nen,m)
integer :: i,j
double precision, parameter :: w22=.5d0*dsqrt(2.d0)
double precision :: vvv(nTens,m)

oo = 0.d0

vvv = vv

if (ndm==2) then 
  vvv(3,:)=w22*vvv(3,:)
  do j = 1,m
    do i = 1,nen
    oo(2*i-1,j) = BB(1,i)*vvv(1,j) + BB(2,i)*vvv(3,j)
    oo(2*i,j)   = BB(2,i)*vvv(2,j) + BB(1,i)*vvv(3,j)
    end do 
  end do 
else 
  vvv(4:6,:)=vvv(4:6,:)*w22
  do j = 1,m
    do i = 1,nen
      oo(3*i-2,j) = BB(1,i)*vvv(1,j) + BB(3,i)*vvv(5,j) + BB(2,i)*vvv(6,j)
      oo(3*i-1,j) = BB(2,i)*vvv(2,j) + BB(3,i)*vvv(4,j) + BB(1,i)*vvv(6,j)
      oo(3*i,j) = BB(3,i)*vvv(3,j) + BB(2,i)*vvv(4,j) + BB(1,i)*vvv(5,j)
    end do
  end do 
end if

end subroutine getBEpsTMat




subroutine getNPhiTMat(oo,NN,vv,nen,nPhi,m)
implicit none
integer, intent(in) :: nen,nPhi,m
double precision, intent(in) :: NN(nen),vv(nPhi,m)
double precision, intent(out) :: oo(nPhi*nen,m)
integer :: I,j,k,l

do k=1,m
  do I=1,nen
    l=nPhi*(I-1)
    do j=1,nPhi
      oo( l+j , k )=NN(I)*vv(j,k)
    end do
  end do
end do

end subroutine getNPhiTMat




subroutine getBPhiTMat(oo,BB,vv,nen,nPhi,ndm,m)
implicit none
integer, intent(in) :: nen,nPhi,ndm,m
double precision, intent(in) :: BB(ndm,nen),vv(ndm*nPhi,m)
double precision, intent(out) :: oo(nPhi*nen,m)
integer :: I,j,k,id,pos1, pos2

oo=0.d0
do k=1,m
  do I=1,nen
    pos1=nPhi*(I-1)
    do j=1,nPhi
      pos2=ndm*(j-1)
      do id=1,ndm
        oo(pos1+j,k)=oo(pos1+j,k)+BB(id,I)*vv(pos2+id,k)
      end do
    end do      
  end do
end do

end subroutine getBPhiTMat



! generate_K(KMat,SHP(1:ndm,nen),ca,nen,nTens,ndm)
subroutine generate_K(BTCB,BB,ca,nen,nTens,ndm)
implicit none

integer, intent(in) :: nen,nTens,ndm
double precision, intent(in) :: BB(ndm,nen),ca(nTens,nTens)
double precision :: BTCT(ndm*nen,nTens)!, CB(nTens,ndm*nen)
double precision, intent(out) :: BTCB(ndm*nen,ndm*nen)

    !getBEpsTMat(oo(ndm*nen,m),SHP(1:ndm,nen),vv(nTens,m),nen,nTens,ndm,m)
!     caT=transpose(ca)
    call getBEpsTMat(BTCT,BB,transpose(ca),nen,nTens,ndm,nTens)!SW2: I think you don't need the variables caT and CB. You can simply do 'call getBEpsTMat(BTCT,BB,transpose(ca),nen,nTens,ndm,nTens)', similar below
!     CB=transpose(BTCT) !CBT(nTens,ndm*nen)
    call getBEpsTMat(BTCB,BB,transpose(BTCT),nen,nTens,ndm,ndm*nen)!SW2: similar here
    
end subroutine generate_K



! UmatInterface(d,FF,ndm,nTens,tau_Voigt,ca,hr(pos1:pos1+nHistVars-1),hr(pos2:pos2+nHistVars-1),nHistVars)
subroutine UmatInterface(d,FF,ndm,nTens,tau,ca,hn,h1,nHistVars)
implicit none 

integer, intent(in) :: ndm, nTens, nHistVars 
double precision, intent(in) :: FF(ndm,ndm), d(*),hn(nHistVars)
double precision, intent(out) :: tau(nTens), ca(nTens, nTens),h1(nHistVars)

double precision :: FF_3D(3,3),E_Green_Lagr(3,3),E_Green_Lagr_Voigt(6),SS_Voigt_3D(6),FTBOXF(6,6),detF,CC(6,6),ca_3D(6,6), tau_3D(6)
double precision, parameter :: w2=dsqrt(2.d0)
integer :: i

if(ndm.eq.2)then 
  FF_3D(1:2,1:2)=FF(1:2,1:2)     !---plane strain (EVZ)---!
  FF_3D(1:2,3)=0.d0; FF_3D(3,1:2)=0.d0; FF_3D(3,3)=1.d0
else
  FF_3D=FF  
end if

E_Green_Lagr = 0.5d0*matmul(transpose(FF_3D),FF_3D) !E_Green_Lagrange 
do i=1,3;  E_Green_Lagr(i,i)=E_Green_Lagr(i,i)-0.5d0; end do

call trnsfrmMat2Vec(E_Green_Lagr, E_Green_Lagr_Voigt,3,6)
call SmallStrainUmat(d,E_Green_Lagr_Voigt,SS_Voigt_3D,CC)
call getATBOXSA(FF_3D, FTBOXF) 

tau_3D=matmul(transpose(FTBOXF),SS_Voigt_3D)

if(ndm.eq.2)then 
  tau(1:2)=tau_3D(1:2);  tau(3)=tau_3D(6)
else
  tau=tau_3D
end if

 ca_3D=matmul(matmul(transpose(FTBOXF),CC),FTBOXF)

if(ndm.eq.2)then 
   ca(1:2,1:2) = ca_3D(1:2,1:2); ca(1:2,3)   = ca_3D(1:2,6)
   ca(3,1:2)   = ca_3D(6,1:2); ca(3,3)     = ca_3D(6,6)
else
   ca=ca_3D
end if 
 
h1(2:nTens+1)=tau(1:nTens) !stress in History
end subroutine UmatInterface
!SW2: Here, you did a really good job!






subroutine UmatInterfaceGradExt(d,FF,phi,gradPhiVoigt,ndm,nTens,nPhi,tau,pi,xi,ca,cDPi,dTauDPhi,dPiDPhi,dXiDGPhi,hn,h1,nHistVars, isw)
implicit none 

integer, intent(in) :: ndm, nTens, nPhi, nHistVars, isw
double precision, intent(in) :: FF(ndm,ndm), phi(nPhi), d(*), gradPhiVoigt(ndm*nPhi), hn(nHistVars)
double precision, intent(out) :: tau(nTens), pi(nPhi), xi(ndm*nPhi), h1(nHistVars), ca(nTens, nTens), cDPi(nPhi,nTens)
double precision, intent(out) ::  dTauDPhi(nTens,nPhi), dPiDPhi(nPhi,nPhi), dXiDGPhi(ndm*nPhi,ndm*nPhi)

double precision :: FF_3D(3,3),E_Green_Lagr(3,3),E_Green_Lagr_Voigt(6),SS_Voigt_3D(6),FTBOXF(6,6),CC(6,6),ca_3D(6,6), tau_3D(6)
double precision :: gradPhiVoigt_3D(3*nPhi), xi_3D(3*nPhi), cDPi_3D(nPhi,6), dTauDPhi_3D(6,nPhi), dXiDGPhi_3D(3*nPhi,3*nPhi)
double precision, parameter :: w2=dsqrt(2.d0)
integer :: i,j,k

if(ndm.eq.2)then 
  FF_3D(1:2,1:2)=FF(1:2,1:2)     !---plane strain (EVZ)---!
  FF_3D(1:2,3)=0.d0; FF_3D(3,1:2)=0.d0; FF_3D(3,3)=1.d0
  do i=1,nPhi
    j=3*(i-1); k=2*(i-1)
    gradPhiVoigt_3D(j+1)=gradPhiVoigt(k+1)
    gradPhiVoigt_3D(j+2)=gradPhiVoigt(k+2)
    gradPhiVoigt_3D(j+3)=0.d0
  end do
else
  FF_3D=FF  
  gradPhiVoigt_3D=gradPhiVoigt
end if
! call schreibeb(FF,3,3,'Grad')
E_Green_Lagr = 0.5d0*matmul(transpose(FF_3D),FF_3D) !E_Green_Lagrange 
do i=1,3;  E_Green_Lagr(i,i)=E_Green_Lagr(i,i)-0.5d0; end do
! call schreibeb(E_Green_Lagr,3,3,'Green')

call trnsfrmMat2Vec(E_Green_Lagr, E_Green_Lagr_Voigt,3,6)

! call SmallStrainUmat(d,E_Green_Lagr_Voigt,SS_Voigt_3D,CC)

! subroutine umatIsoGradDam(eps,DChi,h1,hn,sig,dSigDeps,&
!     &dSigDDChi,dPiDDChi,gradDChi,Xi,dXiDGradDchi,nHistVars)

SS_Voigt_3D=0.d0; CC=0.d0; dTauDPhi_3D=0.d0; dPiDPhi=0.d0; xi_3D=0.d0; dXiDGPhi_3D=0.d0; cDPi_3D=0.d0; pi=0.d0
! call umatIsoGradDam(E_Green_Lagr_Voigt,phi(1),h1(1:17),hn(1:17),SS_Voigt_3D,CC,&
!     &dTauDPhi_3D(:,1),pi(1),dPiDPhi(1,1),gradPhiVoigt_3D(1:3),xi_3D(1:3),dXiDGPhi_3D(1:3,1:3),17)
! call umatIsoGradDamFibre(E_Green_Lagr_Voigt,phi(2),h1(18:19),hn(18:19),SS_Voigt_3D,CC,&
!     &dTauDPhi_3D(:,2),pi(2),dPiDPhi(2,2),gradPhiVoigt_3D(4:6),xi_3D(4:6),dXiDGPhi_3D(4:6,4:6),2)

call umatIsoGradDamFibre(E_Green_Lagr_Voigt,phi(1),h1,hn,SS_Voigt_3D,CC,&
    &dTauDPhi_3D(:,1),pi(1),dPiDPhi(1,1),gradPhiVoigt_3D(1:3),xi_3D(1:3),dXiDGPhi_3D(1:3,1:3),19,d)
if (d(30).lt.0.5) then
  !write(*,*) 'NoCC'
  call umatIsoGradDam(E_Green_Lagr_Voigt,phi(2),h1(1:17),hn(1:17),SS_Voigt_3D,CC,&
      &dTauDPhi_3D(:,2),pi(2),dPiDPhi(2,2),gradPhiVoigt_3D(4:6),xi_3D(4:6),dXiDGPhi_3D(4:6,4:6),17,d)

else if (d(30).gt.0.5) then
  !write(*,*) 'CC'
  call umatIsoGradDam_CC(E_Green_Lagr_Voigt,phi(2),h1(1:17),hn(1:17),SS_Voigt_3D,CC,&
      &dTauDPhi_3D(:,2),pi(2),dPiDPhi(2,2),gradPhiVoigt_3D(4:6),xi_3D(4:6),dXiDGPhi_3D(4:6,4:6),17,d)
      
!   write(*,*) 'Stiff:', CC(1,1), CC(1,2), CC(4,4) 
  !call schreibeb(CC,6,6,'Stif')
  
      
else
  write(*,*) 'ERROR: No Matrix material selected'
  
end if

! if (isw.eq.6) then
!   call schreibeb(xi_3D,6,1,'Xi :')
! !   call schreibeb(dXiDGPhi_3D,6,6,'DXidD')
!   call schreibeb(dPiDPhi(2,2),1,1,'dPiDDX')
!   call schreibeb(dTauDPhi_3D(:,2),6,1,'dSDDX')
!   call schreibeb(SS_Voigt_3D,6,1,'Sig:')
!    call schreibeb(CC(:,:),6,6,'dSDEps')
! endif

!     subroutine umatIsoGradDamFibre(eps,DChi,h1,hn,sig,dSigDeps,&
!     &dSigDDChi,pi,dPiDDChi,gradDChi,Xi,dXiDGradDChi,nHistVars)

call getATBOXSA(FF_3D, FTBOXF) 
tau_3D=matmul(transpose(FTBOXF),SS_Voigt_3D)
! dTauDPhi_3D(:,1)=matmul(transpose(FTBOXF),dTauDPhi_3D(:,1))
!  cDPi_3D(1,:)=dTauDPhi_3D(:,1)
! dTauDPhi_3D(:,2)=matmul(transpose(FTBOXF),dTauDPhi_3D(:,2))
!  cDPi_3D(2,:)=dTauDPhi_3D(:,2)
 
dTauDPhi_3D=matmul(transpose(FTBOXF),dTauDPhi_3D)
 cDPi_3D=transpose(dTauDPhi_3D)

if(ndm.eq.2)then 
  tau(1:2)=tau_3D(1:2);  tau(3)=tau_3D(6)
  do i=1,nPhi
    j=2*(i-1); k=3*(i-1)
    xi(j+1)=xi_3D(k+1)
    xi(j+2)=xi_3D(k+2)
  end do
else
  tau=tau_3D
  xi=xi_3D
end if

 ca_3D=matmul(matmul(transpose(FTBOXF),CC),FTBOXF)

if(ndm.eq.2)then 
   dXiDGPhi=0.d0 
   ca(1:2,1:2) = ca_3D(1:2,1:2); ca(1:2,3)   = ca_3D(1:2,6)
   ca(3,1:2)   = ca_3D(6,1:2); ca(3,3)     = ca_3D(6,6)
   cDPi(:,1:2) = cDPi_3D(:,1:2); cDPi(:,3)=cDPi_3D(:,6)
   dTauDPhi(1:2,:)=dTauDPhi_3D(1:2,:); dTauDPhi(3,:)=dTauDPhi_3D(6,:)
  do i=1,nPhi
    do j=1,nPhi
      dXiDGPhi(2*i-1:2*i, 2*j-1:2*j) = dXiDGPhi_3D(3*i-2:3*i-1, 3*j-2:3*j-1)
    end do
  end do
else
   ca=ca_3D; cDPi=cDPi_3D; dTauDPhi=dTauDPhi_3D; dXiDGPhi=dXiDGPhi_3D
end if 
 
! h1(2:nTens+1)=tau(1:nTens) !stress in History
end subroutine UmatInterfaceGradExt







subroutine getATBOXSA(aa, ATBOXSA)
implicit none

!     Declarations

integer::i,j
double precision, DIMENSION (3,3) :: aa
integer, DIMENSION (6,2) :: q

double precision,dimension(6)::p
double precision, DIMENSION (6,6) :: ATBOXSA

!----------------------------------------------------------------------------

q=reshape((/1,2,3,2,1,1,1,2,3,3,3,2/),(/6,2/))
!
p=(/1.d0,1.d0,1.d0,sqrt(2.d0),sqrt(2.d0),sqrt(2.d0)/)
!
do i=1,6
   do j=1,6
      ATBOXSA(i,j) = 0.5d0*(p(i)*p(j)*(aa(q(j,1),q(i,1))*aa(q(j,2),q(i,2))+aa(q(j,2),q(i,1))*aa(q(j,1),q(i,2))))
    end do
end do


!         write(*,*)'ATBOXSA is ', ATBOXSA
end subroutine


! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc]
      subroutine schreibeb(a,nz,ns,wort)
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc]
      implicit none
      integer n3s,ns,j,i,k,nz
      double precision a(nz,ns)
      character*6 wort
! c
! c.... Subroutine for outputting a matix on the screen
! c
      n3s=ns/3
      write(*,*)wort
      if(n3s.lt.1)goto 20
      do j=0,n3s-1
	write(*,*)'Spalten ',j*3+1,'-',j*3+3
	do i=1,nz
	  write(*,*)(a(i,j*3+k),k=1,3)
        end do
      end do
      if (ns.gt.3*n3s) then
	write(*,*)'Spalten ',3*n3s+1,'-',ns
	do i=1,nz
	  write(*,*)(a(i,k),k=3*n3s+1,ns)
        end do
      endif
      return
20    do i=1,nz
	  write(*,*)(a(i,j),j=1,ns)
      end do
! c
      return
      end

