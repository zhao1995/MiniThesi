      subroutine mate3d16(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
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

c     Outputs:
c         md = 8      - number of used data for control of d-array 
c         nh = 8      - number of history parameter at Gauss-Point
c         h2(nh)      - history array h2
c         Sig         - 2nd Piola-Kirchhoff Stress
c         Cmat        - algorithimic consistent tangent modulus
c
c     (c) f.k                                                   Feb,2012
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),Eps(*),Sig(*),Cmat(6,6),E_p(6),
     +          plout(10),xgp(3)
      if(isw.eq.1) then
c....   input data
        call mati3d16(d,md,nh)
c
      else
        xm   = d(1)/(2.d0*(1.d0+d(2))) 
        xd   = d(6)
        eta  = d(7)
c....   Get plastic variables from history variables at time n 
        do i = 1,6
           E_p(i) = h1(i)
        end do
        a = h1(7)
        ab= h1(8)
c
c       compute elastic-plastic tangent modulus and stresses
        call elas3d16(d,Eps,E_p,a,ab,dvp,sig,cmat)
c
c....   store history variables at time n+1 in  h2-array
        do i = 1,6
           h2(i) = E_p(i)
           plout(i)=h2(i)
        end do
        h2(7) = a
        h2(8) = ab
        plout(7)=h2(7)
        plout(8)=h2(8)
        valuse1 = a
        valuse2 = ab
      end if
c
      return
      end
c
      subroutine mati3d16(d,md,nh)
c-----------------------------------------------------------------------
c
c     mate3d16 small elastoplastic strains adhesive E = Ee + Ep
c     input material parameter 
c     
c-----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      dimension d(30)

      md1=9
      nh=8
      if(ior.lt.0) write(*,1001)
1001  format(
     + ' Input: 1')
      call dinput(d(1),md1)
c
      md2=9
      if(ior.lt.0) write(*,1002)
1002  format(
     + ' Input 2')
      call dinput(d(10),md2)
c
      md3=6
      if(ior.lt.0) write(*,1003)
1003  format(
     + ' Input 3')
      call dinput(d(19),md3)
c 
                  write(iow,1004) nh,(d(i),i=1,24)
      if(ior.lt.0)write(*  ,1004) nh,(d(i),i=1,24)
1004  format(5x,'Adhesive material elastoplastic material data',/,
     +  5x,'length nh of h1,h2 ...........',i12,/,
     +  5x,'Youngs modulus E .............',g12.5,/,
     +  5x,'Poissons ratio v .............',g12.4,/,
     +  5x,'Initial yield stress    y_0 ..',g12.5,/,
     +  5x,'Parameter a1 .................',g12.5,/,
     +  5x,'Parameter a2 .................',g12.5,/,
     +  5x,'Parameter a1star .............',g12.5,/,
     +  5x,'Parameter a2star .............',g12.5,/,
     +  5x,'Parameter q...................',g12.5,/,
     +  5x,'Parameter B...................',g12.5,/,
     +  5x,'Parameter H ..................',g12.4,/,
     +  5x,'Parameter evol................',g12.5,/,
     +  5x,'Parameter dvol ...............',g12.5,/,
     +  5x,'Parameter wvol ...............',g12.5,/,
     +  5x,'Parameter voln ...............',g12.5,/,
     +  5x,'Parameter edev................',g12.5,/,
     +  5x,'Parameter ddev ...............',g12.5,/,
     +  5x,'Parameter wdev ...............',g12.5,/,
     +  5x,'Parameter devn ...............',g12.5,/,
     +  5x,'Parameter cpdev...............',g12.5,/,
     +  5x,'Parameter betapdev............',g12.5,/,
     +  5x,'Parameter gammadev............',g12.5,/,
     +  5x,'Parameter cpvol...............',g12.5,/,
     +  5x,'Parameter betavol.............',g12.5,/,
     +  5x,'Parameter gammavol............',g12.5,/)
c
      md=md1+md2+md3
      return
      end 
c
      subroutine elas3d16(d,ec,ep,a,ab,dvp,sig,cmat)
c-----------------------------------------------------------------------
c
c     mate3d16 small elastoplastic strains isotropic E = Ee + Ep
c     calculate C, Sigma_V 
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),ep(6),sig(6),cmat(6,6),ee(6),sd(6),ec(6),cmat1(6,6)
     +          ,cmat2(6,6),ep1(6),sig1(6),xvar(6,6)
     +          ,testana(9,9),testloc(9,9),testit(9,9)
c    
      xchl  = dvp**(1.d0/3.d0)
      xE   = d(1)
      xn   = d(2)
      y_o  = d(3)
      a1   = d(4)
      a2   = d(5)
      a1s  = d(6)
      a2s  = d(7)
      q    = d(8)
      b    = d(9)
      xH   = d(10)
      evol = d(11)
      dvol = d(12)
      wvol = d(13)
      voln = d(14)
      edev = d(15)
      ddev = d(16)
      wdev = d(17)
      devn = d(18)
      br   = d(19)
      qr   = d(20)
      xhr  = d(21)
      epsg = d(22)
      expn = d(23)
      dk   = d(24)
      ep1=ep
      ah=a
      bh=ab
      an=ab+a
c
c     damage volumetric
      if (a.le.evol) then
        Wvolt=1.d0
      else
        Wvolt=exp(-wvol*(a-evol)**voln)       
      endif
c     damage deviatoric
      if (ab.le.edev) then
        Wdevt=1.d0
      else
        Wdevt=exp(-wdev*(ab-edev)**devn)       
      endif
      W = Wvolt*Wdevt
c      W=1.d0
c.....elasticity matrix
      call elas3d16s(Cmat,xE,xn)
c      
c.... trial strains  E^el = E - E^pl
      ee = ec - ep
c
c.....trial stresses   sig = cmat*ee
      sig = matmul(Cmat,ee)
      if (dabs(y_o/xE).lt.1.e-15) return
c
c.....yield stress 
cfk      y = y_o + xk*a + (y_i-y_o)*(1-dexp(-xd*a))
c
c.... compute g^trial
      trsig = (sig(1)+sig(2)+sig(3))
      sd(1) = sig(1) - trsig/3.d0 
      sd(2) = sig(2) - trsig/3.d0
      sd(3) = sig(3) - trsig/3.d0
      sd(4) = sig(4)
      sd(5) = sig(5)
      sd(6) = sig(6)
      trsig = trsig/W
      g_tr  =  (1.d0/2.d0*(sd(1)*sd(1)+sd(2)*sd(2)+sd(3)*sd(3))
     +       + (sd(4)*sd(4)+sd(5)*sd(5)+sd(6)*sd(6)))/W**2.d0 
c      
c.....Yield condition   f = g^trial - yz
      f = g_tr -3.d0**(-1)*((y_o+q*(1-exp(-b*an))+xH*an)**2.d0
     +         -a1*y_o*trsig-a2*trsig**2.d0)
      if (f.gt.1.d-10) then
         call plas3d16ana(d,sig,ec,a,ab,dvp,ep,cmat)      
      endif
      return
c
      end
     
      subroutine elas3d16s(cmat,e1,xnu)
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension cmat(6,6)
c      include 'tdata.h'
c.... linear isotropic elasticity matrix
c
      cc = e1/((1.0d0+xnu)*(1.0d0-2.0d0*xnu))
      tt = (1.0d0-2.0d0*xnu)/2.0d0
      cc1xnu = cc*(1.0d0-xnu)
      ccxnu = cc*xnu
      cctt  = cc*tt 
c
      cmat(1,1) = cc1xnu
      cmat(1,2) = ccxnu
      cmat(1,3) = ccxnu
      cmat(2,2) = cc1xnu
      cmat(2,3) = ccxnu
      cmat(3,3) = cc1xnu
      cmat(4,4) = cctt
      cmat(5,5) = cctt
      cmat(6,6) = cctt
c
c.... symmetry
      cmat(2,1) = cmat(1,2)
      cmat(3,1) = cmat(1,3)
      cmat(3,2) = cmat(2,3)
c
      return
      end

      subroutine plas3d16ana(d,asig,ae,a,ab,dvp,aep,cmat)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),ae(6),e(6,1),aep(6),ep(6,1)
     +          ,asig(6),sig(6,1),sd(6,1),sds(6,1),sdh(6,1)
     +          ,xR(6,1),cmat(6,6),cmati(6,6),xA(6,6),xAi(6,6)          
     +          ,xGa(9,9),xGai(9,9),xRa(9,1),dR(9,1),xQv(6,1),xQd(6,1)
     +          ,xdQds(6,1),xdFds(6,1),xQdevo(6,1),xQdede(6,1)
     +          ,xPo(6,6),xP(6,6),xIde(6,1),xtest(9,9),cmathlp(6,6)
     +          ,xtestr(9,9),aephlp(6)
      LOGICAL exist
c     
      xP=0.d0;xPo=0.d0;xGa=0.d0
      cmathlp=cmat
c
      do i=1,6
      sig(i,1) = asig(i)
      e(i,1)   = ae(i)
      ep(i,1)  = aep(i)
      xP(i,i)  = 1.d0
      xPo(i,i)  = 1.d0
      if (i>3) xIde(i,1)=0.d0
      if (i>3) xP(i,i)  =2.0d0
      if (i<4) then
      xIde(i,1)=1.d0
      do j=1,3
          xP(i,j)=xP(i,j)-1.d0/3.d0
          xPo(i,j)=xPo(i,j)-1.d0/3.d0  
      enddo
      endif
      enddo
c
      xk = d(1)/(3.d0*(1.d0-2.d0*d(2)))       ! bulk modulus 
      xm = d(1)/(2.d0*(1.d0+d(2)))            ! shear modulus
      xchl  = dvp**(1.d0/3.d0)
      xE   = d(1)
      xn   = d(2)
      y_o  = d(3)
      a1   = d(4)
      a2   = d(5)
      a1s  = d(6)
      a2s  = d(7)
      q    = d(8)
      b    = d(9)
      xH   = d(10)
      evol = d(11)
      dvol = d(12)
      wvol = d(13)
      voln = d(14)
      edev = d(15)
      ddev = d(16)
      wdev = d(17)
      devn = d(18)
      br   = d(19)
      qr   = d(20)
      xhr  = d(21)
      epsg = d(22)
      expn = d(23)
      dk   = d(24)
c.....start values
      da  = a+ab
      call pivot(cmat,6,6,cmati)
c.....start values
      gam = 0.d0
      dgam = 0.d0
      qa = a
      qb  = ab
      newton = 0
c
      xR=0.d0    
c.....begin local Newton iteration...................................
100   continue
      xdQds=0.d0;xdFds=0.d0;
      trsig = (sig(1,1)+sig(2,1)+sig(3,1))
      sd(1,1) = sig(1,1) - trsig/3.d0 
      sd(2,1) = sig(2,1) - trsig/3.d0
      sd(3,1) = sig(3,1) - trsig/3.d0
      sd(4,1) = sig(4,1)
      sd(5,1) = sig(5,1)
      sd(6,1) = sig(6,1)
      g_tr  =  1.d0/2.d0*(sd(1,1)*sd(1,1)+sd(2,1)*sd(2,1)
     +                   +sd(3,1)*sd(3,1))
     +       + (sd(4,1)*sd(4,1)+sd(5,1)*sd(5,1)+sd(6,1)*sd(6,1))
      newton = newton + 1 
      if (newton.gt.300) then
         write( * ,*) 'no convergence in local newton'
         write(iow,*) 'no convergence in local newton'
         call mate16abort('convnewton')
      end if
      if (a.le.evol) then
        Wvolt=1.d0
        wvolpkt=0.d0 
      else
        Wvolt=exp(-wvol*(a-evol)**voln)
        wvolpkt=-wvolt*voln*wvol*
     +         (a-evol)**(voln-1.d0)      
      endif
c     damage deviatoric
      if (ab.le.edev) then
        Wdevt=1.d0
        wdevpkt=0.d0
      else
        Wdevt=exp(-wdev*(ab-edev)**devn)
        Wdevpkt=-wdevt*devn*wdev*
     +         (ab-edev)**(devn-1.d0) 
      endif
      W = Wvolt*Wdevt
      xdQds(1,1)=(a1s*y_o/W+2.d0*a2s*trsig/W**2.d0)/3.d0
      xdQds(2,1)=(a1s*y_o/W+2.d0*a2s*trsig/W**2.d0)/3.d0
      xdQds(3,1)=(a1s*y_o/W+2.d0*a2s*trsig/W**2.d0)/3.d0
      xdFds(1,1)=(a1*y_o/W+2.d0*a2*trsig/W**2.d0)/3.d0
      xdFds(2,1)=(a1*y_o/W+2.d0*a2*trsig/W**2.d0)/3.d0
      xdFds(3,1)=(a1*y_o/W+2.d0*a2*trsig/W**2.d0)/3.d0
      do i=1,3
        xdFds(i,1)=xdFds(i,1)+sd(i,1)/W**2.d0
        xdQds(i,1)=xdQds(i,1)+sd(i,1)/W**2.d0
        sds(i,1)=sd(i,1)
        sdh(i,1)=sd(i,1)
      enddo
      do i=4,6
        sds(i,1)=sd(i,1)*2.d0
        sdh(i,1)=sd(i,1)*.5d0
        xdFds(i,1)=sd(i,1)/W**2.d0*2.d0
        xdQds(i,1)=sd(i,1)/W**2.d0*2.d0
      enddo
      xA=matmul(xIde,transpose(xIde))*2.d0*a2s/3.d0/W**2.d0
     +         +xP*1.d0/W**2.d0
      R=y_o+q*(1-exp(-b*da))+xH*da
      Rs=q*b*exp(-b*da)+xH
      dFdedev=-2.d0/3.d0*R*Rs+(-2.d0*g_tr/W**3.d0-1.d0/3.d0*
     +     (a1*y_o*trsig/W**2.d0+2.d0*a2*trsig**2.d0/W**3.d0))
     +      *Wvolt*Wdevpkt
      dFdevol=-2.d0/3.d0*R*Rs+(-2.d0*g_tr/W**3.d0-1.d0/3.d0*
     +     (a1*y_o*trsig/W**2.d0+2.d0*a2*trsig**2.d0/W**3.d0))
     +     *Wdevt*Wvolpkt
      xQv=gam/R*(matmul(matmul(xIde,transpose(xIde))*1.d0/3.d0,xdQds)
     +         +matmul(xA,trsig/3.d0*xIde))
      xQd=gam/R*(matmul(xPo,xdQds)+matmul(xA,sd))
      r2=-(g_tr/W**2.d0-(R**2.d0
     +        -a1*y_o*trsig/W-a2*trsig**2.d0/W**2.d0)/3.d0)
      xR=e-ep-matmul(cmati,sig)
     +   -gam*(sds/W**2.d0+(a1s*y_o/W
     +         +2.d0*a2s*trsig/W**2.d0)/3.d0*xIde)
      r3v =  -a+qa+gam/R*sum(matmul(trsig/3.d0*transpose(xIde),xdQds))
      r3d = -ab+qb+gam/R*sum(matmul(transpose(sd),xdQds))
      xQdevo=wvolpkt*Wdevt*(-2.d0/W**3.d0*sds
     +      -(a1s*y_o/W**2.d0+4.d0*a2s/W**3.d0*trsig)/3.d0*xIde)
      xQdede=wdevpkt*Wvolt*(-2.d0/W**3.d0*sds
     +      -(a1s*y_o/W**2.d0+4.d0*a2s/W**3.d0*trsig)/3.d0*xIde)
      xA=gam*xA+cmati
      do i=1,6
        do j=1,6
          xGa(i,j)=xA(i,j)
        enddo 
          xGa(7,i)=xdFds(i,1)
          xGa(8,i)=-xQv(i,1)
          xGa(9,i)=-xQd(i,1)
          xGa(i,7)=xdQds(i,1)
          xGa(i,8)=gam*xQdevo(i,1)
          xGa(i,9)=gam*xQdede(i,1)
          xRa(i,1)=xR(i,1)
      enddo
      xGa(7,7)=0.d0
      xGa(7,8)=dFdevol
      xGa(7,9)=dFdedev
      xGa(8,7)=-sum(matmul(transpose(1.d0/R*xIde*trsig/3.d0),xdQds))
      xGa(9,7)=-sum(matmul(transpose(1.d0/R*sd),xdQds))
      xGa(8,8)=1.d0
     +        -gam/R*sum(matmul(transpose(xIde*trsig/3.d0),xQdevo))
     +        +gam/R**2.d0*sum(matmul(trsig/3.d0*transpose(xIde),xdQds))
     *         *Rs
      xGa(9,8)=-gam/R*sum(matmul(transpose(sd),xQdevo))+
     +        +gam/R**2.d0*sum(matmul(transpose(sd),xdQds))*Rs
      xGa(8,9)=-gam/R*sum(matmul(transpose(xIde*trsig/3.d0),xQdede))
     +        +gam/R**2.d0*sum(matmul(trsig/3.d0*transpose(xIde),xdQds))
     *        *Rs 
      xGa(9,9)=1.d0-gam/R*sum(matmul(transpose(sd),xQdede))+
     +        +gam/R**2.d0*sum(matmul(transpose(sd),xdQds))*Rs 
      xRa(7,1)=r2
      xRa(8,1)=r3v
      xRa(9,1)=r3d
c      call inv(xGa,9,9,xGai)
      call pivot(xGa,9,9,xGai)
      dR=matmul(xGai,xRa)
      do i=1,6
        sig(i,1)=sig(i,1)+dR(i,1)
      enddo
      gam=gam+dR(7,1)
      a=a+dR(8,1)
      ab=ab+dR(9,1)
      da  = a+ab
      resi = sqrt(sum(xRa**2))
      if (resi.gt.1.d-10) goto 100
      continue
      do i=1,6
          aephlp(i)=aep(i)+gam*xdQds(i,1)
c          asig(i)=sig(i,1)
          do j=1,6
              cmat(i,j)=xGai(i,j)
          enddo
      enddo
      call plas3d16iter(d,sig,ae,a,ab,qa,qb,aep,cmathlp,gam,xtest)
      do i=1,9
        do j =1,9
        if (abs(xGa(i,j)).ge.1.d-10) then
            xtestr(i,j)=1.d0-xtest(i,j)/xGa(i,j)
        else   
            xtestr(i,j)=xtest(i,j)-xGa(i,j)
        endif
        enddo
      enddo
      do i=1,6
          asig(i)=sig(i,1)
      enddo
      aep=ae-matmul(cmati,asig)
      if (gam.lt.0.d0) then
         write( * ,*) 'gam decreases'
         write(iow,*) 'no convergence in local newton'
         call mate16abort('gammaerror')
      endif
      return
      end     
      
      subroutine cmat3d16a(d,sig,ec,a,c,ep,cmat1,cmat,dvp,xAi)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),sig(6,1),ep(6,1),ec(6,1),cmat(6,6),cmathlp(6,6),
     +         sighlp(6,1),ephlp(6,1),ehlp(6,1),cmat1(6,6),xAi(9,9),
     +         xAorg(9,9),xA(9,9)
c
      xAi=0.0d0
      delta=1.d-12
      gam1=0.d0
      ephlp=ep
      ehlp=ec
      ehlp(i,1)=ec(i,1)
      sighlp=matmul(cmat1,(ehlp-ep))
      cmat=cmat1
      b=a
      f=c
      gam=gam1
      call check(d,sighlp,ehlp,b,f,gam1,cmat,dvp,xAorg)
c           
c.....start values
      do i=1,6
        ephlp=ep
        ehlp=ec
        ehlp(i,1)=ec(i,1)+delta;
        sighlp=matmul(cmat1,(ehlp-ep))
        cmat=cmat1
        b=a
        f=c
        gam=0.d0
        call check(d,sighlp,ehlp,b,f,gam,cmat,dvp,xA)
        do j=1,6
          xAi(j,i)=(sighlp(j,1)-sig(j,1))/delta
        enddo
        xAi(i,7)=(gam-gam1)/delta
        xAi(i,8)=(b-a)/delta
        xAi(i,9)=(f-c)/delta
      enddo
        ephlp=ep
        ehlp=ec
        sighlp=matmul(cmat1,(ehlp-ep))
        cmat=cmat1
        b=a
        f=c
        gam=delta
        call check(d,sighlp,ehlp,b,f,gam,cmat,dvp,xA)
        do j=1,6
          xAi(j,7)=(sighlp(j,1)-sig(j,1))/delta
        enddo
        xAi(8,7)=(gam-gam1)/delta
        xAi(8,8)=(b-a)/delta
        xAi(8,9)=(f-c)/delta
        ephlp=ep
        ehlp=ec
        sighlp=matmul(cmat1,(ehlp-ep))
        cmat=cmat1
        b=a+delta
        f=c
        call check(d,sighlp,ehlp,b,f,gam,cmat,dvp,xA)
        do j=1,6
          xAi(j,8)=(sighlp(j,1)-sig(j,1))/delta
        enddo
        xAi(8,7)=(gam-gam1)/delta
        xAi(8,8)=(b-a)/delta
        xAi(8,9)=(f-c)/delta
        ephlp=ep
        ehlp=ec
        sighlp=matmul(cmat1,(ehlp-ep))
        cmat=cmat1
        b=a
        f=c+delta
        call check(d,sighlp,ehlp,b,f,gam,cmat,dvp,xA)
        do j=1,6
          cmathlp(j,i)=(sighlp(j,1)-sig(j,1))/delta
          xAi(i,j)=cmathlp(j,i)
        enddo
        xAi(9,7)=(gam-gam1)/var
        xAi(9,8)=(b-a)/var
        xAi(9,9)=(f-c)/var

c
      do i=1,9
      do j=1,9
          xAi(i,j)=xAi(i,j)-xAorg(i,j)
      enddo
      enddo
c      
      return
      end
      
      subroutine plas3d16iter(d,asig,ae,a,ab,qa,qb,aep,cmat,gam,xMa)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),sig(6,1),sd(6,1),ep(6,1),cmat(6,6),xA(6,6),xdQ(6,1)
     +          ,xdF(6,1),xR(6,1),cmati(6,6),xAi(6,6),dsig(6,1),e(6,1)
     +          ,xIde(6,1),cmath(6,6),dsd(6,1), asig(6),ae(6),aep(6)
     +          ,xdds(6,6),sig1(6,1),xdds2(6,6),xdD(6,1),xMa(9,9)
     +          ,xRa(9,1),xdelR(9,1),sigv(6,1),xR0(6,1),xMr(7,7)
     +          ,res(1,1)
c
      xIde=0.d0
      do i=1, 6
      sig(i,1) = asig(i)
      e(i,1)   = ae(i)
      ep(i,1)  = aep(i)
      if (i<4) xIde(i,1)=1.d0 
      enddo
c
      xk = d(1)/(3.d0*(1.d0-2.d0*d(2)))       ! bulk modulus 
      xm = d(1)/(2.d0*(1.d0+d(2)))            ! shear modulus
      xE   = d(1)
      xn   = d(2)
      y_o  = d(3)
      a1   = d(4)
      a2   = d(5)
      a1s  = d(6)
      a2s  = d(7)
      q    = d(8)
      b    = d(9)
      xH   = d(10)
      evol = d(11)
      dvol = d(12)
      wvol = d(13)
      voln = d(14)
      edev = d(15)
      ddev = d(16)
      wdev = d(17)
      devn = d(18)
      br   = d(19)
      qr   = d(20)
      xhr  = d(21)
      epsg = d(22)
      expn = d(23)
      dk   = d(24)
c.....start values
c
c.....begin local Newton iteration...................................
100   continue
      newton = newton + 1 
      if (newton.gt.100) then
         write( * ,*) 'no convergence in local newton'
         write(iow,*) 'no convergence in local newton'
         stop
      end if
c     damage volumetric
      var=1.d-12
      varg=1.d-13
      call plasloc3d16(d,sig,ae,a,ab,qa,qb,aep,gam,
     +                       cmat,xR0,r20,r30v,r30d)
c      xr0=0.d0
c      r20=0.d0
c      r30v=0.d0
c      r3d0=0.d0
      do i=1,6
        sigv=sig
        if (abs(sigv(i,1)).gt.1.d-2) then
          var=1.d-10*abs(sigv(i,1))
        else
          var=1.d-12
        endif
        sigv(i,1)=sigv(i,1)+var
        call plasloc3d16(d,sigv,ae,a,ab,qa,qb,aep,gam
     +         ,cmat,xR,r2,r3v,r3d)
        do j=1,6
          xMa(j,i)=(xR(j,1)-xR0(j,1))/var
          xRa(j,1)=-xR0(j,1)
        enddo
          xMa(8,i)=(r3v-r30v)/var
          xMa(9,i)=(r3d-r30d)/var
          xMa(7,i)=(r2-r20)/var  
      enddo
      if (a.gt.1.d-2) then
          var=1.d-10*a
        else
          var=1.d-12
      endif
      avar=a+var
      call plasloc3d16(d,sig,ae,avar,ab,qa,qb,aep,gam,
     +                cmat,xR,r2,r3v,r3d)
      do j=1,6
        xMa(j,8)=(xR(j,1)-xR0(j,1))/var
      enddo
      xMa(8,8)=(r3v-r30v)/var
      xMa(9,8)=(r3d-r30d)/var
      xMa(7,8)=(r2-r20)/var
      if (ab.gt.1.d-2) then
          var=1.d-10*ab
      else
          var=1.d-12
      endif
      abvar=ab+var
      call plasloc3d16(d,sig,ae,a,abvar,qa,qb,aep,gam,
     +                       cmat,xR,r2,r3v,r3d)
      do j=1,6
        xMa(j,9)=(xR(j,1)-xR0(j,1))/var
      enddo
      xMa(8,9)=(r3v-r30v)/var
      xMa(9,9)=(r3d-r30d)/var
      xMa(7,9)=(r2-r20)/var
      if (gam.gt.1.d-2) then
          varg=1.d-10*gam
      else
          varg=1.d-12
      endif
      gamvar=gam+varg
      call plasloc3d16(d,sig,ae,a,ab,qa,qb,aep,gamvar,
     +                       cmat,xR,r2,r3v,r3d)
      do j=1,6
        xMa(j,7)=(xR(j,1)-xR0(j,1))/varg
      enddo
      xMa(8,7)=(r3v-r30v)/varg
      xMa(9,7)=(r3d-r30d)/varg
      xMa(7,7)=(r2-r20)/varg
c      do i=1,7
c        do j=1,7  
c        xMr(i,j)=xMa(i,j)c
c        enddo
c        xMr(7,j)=xMa(7,j)+xMa(8,j)/xMa(8,8)*xMa(7,8)+
c     +        xMa(9,j)*xMa(7,9)/xMa(9,9)
c      enddo
c       r2=r2-r3v/xMa(8,8)*xMa(7,8)-r3d*xMa(7,9)/xMa(9,9)        
c      call pivot(xMa,9,9,xMa)
c.....end local Newton iteration.....................................
      return
      end

      subroutine cmat3d16(d,sig,ec,a,c,ep,cmat1,cmat)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),sig(6,1),ep(6,1),ec(6,1),cmat(6,6),cmathlp(6,6),
     +         sighlp(6,1),ephlp(6,1),ehlp(6,1),cmat1(6,6)   
c
      delta=1.d-13  
c.....start values
      do i=1,6
        ephlp=ep
        ehlp=ec
        ehlp(i,1)=ec(i,1)+delta;
        sighlp=matmul(cmat1,(ehlp-ep))
        cmat=cmat1
        b=a
        f=c
        call plas3d16iter(d,sighlp,ehlp,b,f,ephlp,cmat)
        do j=1,6
          cmathlp(j,i)=(sighlp(j,1)-sig(j,1))/delta
        enddo
      enddo
c
      cmat=cmathlp
c      
      return
      end      
     
      subroutine plasloc3d16(d,sig,ae,a,ab,qa,qb,aep,gam,cmat,xR,r2,
     +                                                       r3v,r3d)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),sig(6,1),sd(6,1),dsd(6,1)
     +          ,e(6,1),ae(6),ep(6,1),aep(6),cmat(6,6),cmati(6,6)
     +          ,xdQ(6,1),xdQv(6,1),xdQd(6,1),xR(6,1),xIde(6,1)
      xIde=0.d0
      do i=1,6
        e(i,1)   = ae(i)
        ep(i,1)  = aep(i)
        if (i<4) xIde(i,1)=1.d0 
      enddo
c
      xdQv=0.d0;xdQd=0.d0;
c      
      xk = d(1)/(3.d0*(1.d0-2.d0*d(2)))       ! bulk modulus 
      xm = d(1)/(2.d0*(1.d0+d(2)))            ! shear modulus
      xE   = d(1) 
      xn   = d(2)
      y_o  = d(3)
      a1   = d(4)
      a2   = d(5)
      a1s  = d(6)
      a2s  = d(7)
      q    = d(8)
      b    = d(9)
      xH   = d(10)
      evol = d(11)
      dvol = d(12)
      wvol = d(13)
      voln = d(14)
      edev = d(15)
      ddev = d(16)
      wdev = d(17)
      devn = d(18)
      br   = d(19)
      qr   = d(20)
      xhr  = d(21)
      epsg = d(22)
      expn = d(23)
      dk   = d(24)
      da  = a+ab
      call pivot(cmat,6,6,cmati)
      xdQds=0.d0;xdFds=0.d0;
      trsig = (sig(1,1)+sig(2,1)+sig(3,1))
      sd(1,1) = sig(1,1) - trsig/3.d0 
      sd(2,1) = sig(2,1) - trsig/3.d0
      sd(3,1) = sig(3,1) - trsig/3.d0
      sd(4,1) = sig(4,1)
      sd(5,1) = sig(5,1)
      sd(6,1) = sig(6,1)
      g_tr  =  1.d0/2.d0*(sd(1,1)*sd(1,1)+sd(2,1)*sd(2,1)
     +                   +sd(3,1)*sd(3,1))
     +       + (sd(4,1)*sd(4,1)+sd(5,1)*sd(5,1)+sd(6,1)*sd(6,1))
      if (a.le.evol) then
        Wvolt=1.d0
      else
        Wvolt=exp(-wvol*(a-evol)**voln)    
      endif
c     damage deviatoric
      if (ab.le.edev) then
        Wdevt=1.d0
      else
        Wdevt=exp(-wdev*(ab-edev)**devn)
      endif
      W = Wvolt*Wdevt 
      trsigeff=trsig/W
      xdQv(1,1)=(a1s*y_o/W+2.d0*a2s*trsigeff/W)/3.d0
      xdQv(2,1)=(a1s*y_o/W+2.d0*a2s*trsigeff/W)/3.d0
      xdQv(3,1)=(a1s*y_o/W+2.d0*a2s*trsigeff/W)/3.d0
      do i=1,3
        xdQd(i,1)=sd(i,1)/W**2.d0
        dsd(i,1)=sd(i,1)
      enddo
      do i=4,6
        dsd(i,1)=sd(i,1)*2.d0
        xdQd(i,1)=sd(i,1)/W**2.d0*2.d0
      enddo
      xdQ=xdQd+xdQv
c.....start values
      r2=g_tr/W**2.d0-((y_o+q*(1-exp(-b*da))+xH*da)**2.d0
     +        -a1*y_o*trsigeff-a2*trsigeff**2.d0)/3.d0
      xR=-e+ep+matmul(cmati,sig)+gam*(dsd/W**2.d0+3.d0**(-1.d0)
     +   *(a1s*y_o/W+2.d0*a2s*trsig/W**2.d0)*xIde)
      r3v =  a-qa-gam/(y_o+q*(1-exp(-b*da))+xH*da)
     *      *sum(matmul(trsig/3.d0*transpose(xIde),xdQv)) 
      r3d = ab-qb-gam/(y_o+q*(1-exp(-b*da))+xH*da)
     *      *sum(matmul(transpose(sd),xdQd))
      return
      end
   
      subroutine check(d,asig,ae,a,ab,gam,cmat,dvp,xGai)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),ae(6),e(6,1),aep(6),ep(6,1)
     +          ,asig(6),sig(6,1),sd(6,1),sds(6,1),sdh(6,1)
     +          ,xR(6,1),cmat(6,6),cmati(6,6),xA(6,6),xAi(6,6)          
     +          ,xGa(9,9),xGai(9,9),xRa(9,1),dR(9,1),xQv(6,1),xQd(6,1)
     +          ,xdQds(6,1),xdFds(6,1),xQdevo(6,1),xQdede(6,1)
     +          ,xPo(6,6),xP(6,6),xIde(6,1)
c     
      xP=0.d0;xPo=0.d0;xGa=0.d0
c
      do i=1,6
      sig(i,1) = asig(i)
      e(i,1)   = ae(i)
      ep(i,1)  = aep(i)
      xP(i,i)  = 1.d0
      xPo(i,i)  = 1.d0
      if (i>3) xIde(i,1)=0.d0
      if (i>3) xP(i,i)  =2.0d0
      if (i<4) then
      xIde(i,1)=1.d0
      do j=1,3
          xP(i,j)=xP(i,j)-1.d0/3.d0
          xPo(i,j)=xPo(i,j)-1.d0/3.d0  
      enddo
      endif
      enddo
c
      xk = d(1)/(3.d0*(1.d0-2.d0*d(2)))       ! bulk modulus 
      xm = d(1)/(2.d0*(1.d0+d(2)))            ! shear modulus
      xchl  = dvp**(1.d0/3.d0)
      xE   = d(1)
      xn   = d(2)
      y_o  = d(3)
      a1   = d(4)
      a2   = d(5)
      a1s  = d(6)
      a2s  = d(7)
      q    = d(8)
      b    = d(9)
      xH   = d(10)
      evol = d(11)
      dvol = d(12)
      wvol = d(13)
      voln = d(14)
      edev = d(15)
      ddev = d(16)
      wdev = d(17)
      devn = d(18)
      br   = d(19)
      qr   = d(20)
      xhr  = d(21)
      epsg = d(22)
      expn = d(23)
      dk   = d(24)
c.....start values
      da  = a+ab
      call pivot(cmat,6,6,cmati)
c.....start values
      dgam = 0.d0
      qa = a
      qb  = ab
      newton = 0
c
      xR=0.d0    
c.....begin local Newton iteration...................................
100   continue
      xdQds=0.d0;xdFds=0.d0;
      trsig = (sig(1,1)+sig(2,1)+sig(3,1))
      sd(1,1) = sig(1,1) - trsig/3.d0 
      sd(2,1) = sig(2,1) - trsig/3.d0
      sd(3,1) = sig(3,1) - trsig/3.d0
      sd(4,1) = sig(4,1)
      sd(5,1) = sig(5,1)
      sd(6,1) = sig(6,1)
      g_tr  =  1.d0/2.d0*(sd(1,1)*sd(1,1)+sd(2,1)*sd(2,1)
     +                   +sd(3,1)*sd(3,1))
     +       + (sd(4,1)*sd(4,1)+sd(5,1)*sd(5,1)+sd(6,1)*sd(6,1))
      newton = newton + 1 
      if (newton.gt.300) then
         write( * ,*) 'no convergence in local newton'
         write(iow,*) 'no convergence in local newton'
         stop
      end if
      if (a.le.evol) then
        Wvolt=1.d0
        wvolpkt=0.d0 
      else
        Wvolt=exp(-wvol*(a-evol)**voln)
        wvolpkt=-exp(-wvol*(a-evol)**voln)*voln*wvol*
     +         (a-evol)**(voln-1.d0)      
      endif
c     damage deviatoric
      if (ab.le.edev) then
        Wdevt=1.d0
        wdevpkt=0.d0
      else
        Wdevt=exp(-wdev*(ab-edev)**devn)
        Wdevpkt=-exp(-wdev*(ab-edev)**devn)*devn*wdev*
     +         (ab-edev)**(devn-1.d0) 
      endif
      W = Wvolt*Wdevt
      xdQds(1,1)=3.d0**(-1.d0)*(a1s*y_o/W+2.d0*a2s*trsig/W**2.d0)
      xdQds(2,1)=3.d0**(-1.d0)*(a1s*y_o/W+2.d0*a2s*trsig/W**2.d0)
      xdQds(3,1)=3.d0**(-1.d0)*(a1s*y_o/W+2.d0*a2s*trsig/W**2.d0)
      xdFds(1,1)=(a1*y_o/W+2.d0*a2*trsig/W**2.d0)/3.d0
      xdFds(2,1)=(a1*y_o/W+2.d0*a2*trsig/W**2.d0)/3.d0
      xdFds(3,1)=(a1*y_o/W+2.d0*a2*trsig/W**2.d0)/3.d0
      do i=1,3
        xdFds(i,1)=xdFds(i,1)+sd(i,1)/W**2.d0
        xdQds(i,1)=xdQds(i,1)+sd(i,1)/W**2.d0
        sds(i,1)=sd(i,1)
        sdh(i,1)=sd(i,1)
      enddo
      do i=4,6
        sds(i,1)=sd(i,1)*2.d0
        sdh(i,1)=sd(i,1)*.5d0
        xdFds(i,1)=sd(i,1)/W**2.d0*2.d0
        xdQds(i,1)=sd(i,1)/W**2.d0*2.d0
      enddo
      xA=matmul(xIde,transpose(xIde))*2.d0*a2s/3.d0/W**2.d0
     +         +xP*1.d0/W**2.d0
      R=y_o+q*(1-exp(-b*da))+xH*da
      Rs=q*b*exp(-b*da)+xH
      dFdedev=-2.d0/3.d0*R*Rs+(-2.d0*g_tr/W**3.d0-1.d0/3.d0*
     +     (a1*y_o*trsig/W**2.d0+2.d0*a2*trsig**2.d0/W**3.d0))
     +      *Wvolt*Wdevpkt
      dFdevol=-2.d0/3.d0*R*Rs+(-2.d0*g_tr/W**3.d0-1.d0/3.d0*
     +     (a1*y_o*trsig/W**2.d0+2.d0*a2*trsig**2.d0/W**3.d0))
     +     *Wdevt*Wvolpkt
      xQv=gam/y_o*(matmul(matmul(xIde,transpose(xIde))/3.d0,xdQds)
     +         +matmul(xA,trsig/3.d0*xIde))
      xQd=gam/y_o*(matmul(xPo,xdQds)+matmul(xA,sd))
      r2=-(g_tr/W**2.d0-(R**2.d0
     +        -a1*y_o*trsig/W-a2*trsig**2.d0/W**2.d0)/3.d0)
      xR=e-ep-matmul(cmati,sig)
     +   -gam*(sds/W**2.d0+3.d0**(-1.d0)*(a1s*y_o/W
     +         +2.d0*a2s*trsig/W**2.d0)*xIde)
      r3v =  -a+qa+gam/y_o*sum(matmul(trsig/3.d0*transpose(xIde),xdQds))
      r3d = -ab+qb+gam/y_o*sum(matmul(transpose(sd),xdQds))
      xQdevo=wvolpkt*Wdevt*(-2.d0/W**3.d0*sds
     +      -(a1s*y_o/W**2.d0+4.d0*a2s/W**3.d0*trsig)/3.d0*xIde)
      xQdede=wdevpkt*Wvolt*(-2.d0/W**3.d0*sds
     +      -(a1s*y_o/W**2.d0+4.d0*a2s/W**3.d0*trsig)/3.d0*xIde)
      xA=gam*xA+cmati
      do i=1,6
        do j=1,6
          xGa(i,j)=xA(i,j)
        enddo 
          xGa(7,i)=xdFds(i,1)
          xGa(8,i)=-xQv(i,1)
          xGa(9,i)=-xQd(i,1)
          xGa(i,7)=xdQds(i,1)
          xGa(i,8)=gam*xQdevo(i,1)
          xGa(i,9)=gam*xQdede(i,1)
          xRa(i,1)=xR(i,1)
      enddo
      xGa(7,7)=0.d0
      xGa(7,8)=dFdevol
      xGa(7,9)=dFdedev
      xGa(8,7)=-sum(matmul(transpose(1.d0/y_o*xIde*trsig/3.d0),xdQds))
      xGa(9,7)=-sum(matmul(transpose(1.d0/y_o*sd),xdQds))
      xGa(8,8)=1.d0
     +         -gam/y_o*sum(matmul(transpose(xIde*trsig/3.d0),xQdevo))
      xGa(9,8)=-gam/y_o*sum(matmul(transpose(sd),xQdevo))
      xGa(8,9)=-gam/y_o*sum(matmul(transpose(xIde*trsig/3.d0),xQdede))
      xGa(9,9)=1.d0-gam/y_o*sum(matmul(transpose(sd),xQdede))
      xRa(7,1)=r2
      xRa(8,1)=r3v
      xRa(9,1)=r3d
c     call DPOTRI( UPLO, N, A, LDA, INFO )
      call pivot(xGa,9,9,xGai)      
      dR=matmul(xGai,xRa)
      do i=1,6
        sig(i,1)=sig(i,1)+dR(i,1)
      enddo
      gam=gam+dR(7,1)
      a=a+dR(8,1)
      ab=ab+dR(9,1)
      da  = a+ab
      resi=abs(r2)
      if (resi.gt.1.d-10) goto 100
      continue
      do i=1,6
          aep(i)=aep(i)+gam*xdQds(i,1)
          asig(i)=sig(i,1)
          do j=1,6
              cmat(i,j)=xGai(i,j)
          enddo
      enddo
      return
      end          
      
      subroutine mate16abort(b)
      character(len=10) b
      integer nmbu
!$OMP CRITICAL
      call fexit()
      open(newunit=nmbu,file='kaputt.txt',status='replace')
      write(nmbu,*) b
      close(nmbu)
      stop
!$OMP END CRITICAL
      return
      end
             