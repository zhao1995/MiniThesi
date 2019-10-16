      subroutine elmt24(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)

c   geht nicht fuer mehr als 1 Material wg. SAVE??
c   eingebaut, dass auch ein Verbindungselement möglich ist
c   Eingabe: E=0, G=G, h=-groß Zahl  nur für T_Q  

c-----------------------------------------------------------------------
c.....shear stresses in thin walled beams from Q or M_T  
c     with Ei, Gi formulated for E_ref=1 -> ni = Ei and G_ref=1
c
c     ityp = 1 (default)
c.....from Q  
c     to find M two calculations are necessary
c     M: yM mit Q_y = 0 und Q_z = 1
c     M: zM mit Q_y = 1 und Q_z = 0  
c
c     ityp = 2
c.....from M_T (Bredt)             
c     M_T = M_T0 + Qy*zqm - Qz*yqm (Qy, Qz in S) 
c     omega_s scaled to 1 and is given with respect to M
c     I_T open and closed cross sections
c     I_T für TauG oder nur I_TG??
c     idw = 6 reac?
c     verwoelbung offen? 
c     graphische Darstellung Tau_rand offen + Richtung(=M_T)?
c 
c
c.....(c) w.wagner f. gruttmann 07/01           
c-----------------------------------------------------------------------
c..... d-field 
c      d(1) = Ei
c      d(2) = Gi
c      d(3) = bi
c      d(4) = Qy only for Q 
c      d(5) = Qz only for Q
c      d(6) = ns  number of points for stre      (def=4)  only for Q
c      d(7) = nps number of points for plot,stre (def=10) 
c      d(8) = ityp
c      d(9) = M_T0  only for M_T
c      d(10)= ic = 0 I_T=only I_T_closed
c                = 1 I_T=I_T_open + I_T_closed
c-----------------------------------------------------------------------
      USE bdata       
      USE cdata 
      USE eldata 
      USE fornam 
      USE iofile 
      USE pdata1 
      USE pdata2 
      USE pdata6 
      USE pdata10 
      USE strnam 
      implicit double precision(a-h,o-z)
      logical fltau          
      dimension ix(*),xl(ndm,*),tl(*),d(*),s(nst,*),ul(ndf,*),p(nst),
     +          ixl(5),ixld(4),xll(3,4),xdp(2)
      dimension h1(*),h2(*),h3(*)
      character typtau(2)*6
      save fltau,cphi,ityp,typtau
      save area,area0,az,azz,ay,ayy,ayz,ys,zs,ayq,azq,tm,awy,awz,ait
      save ym,zm,cmt,aito
      data ixl/1,2,3,4,1/,ixld /1,2,3,1/,eps/1.0e-9/,eps1/1.e-7/ 

      data typtau /' Tau_Q','Tau_MT'/
      ielno = 24
c.... go to correct array processor
      go to(1,2,3,4,2,6,2,2,2,2,2,2,13) isw
c.... input material properties
1     if(ior.lt.0) write(*,1001)
1001  format(' Input: E G b Qy Qz ns nps ityp M_T0 I_Toc'/3x,'>',$)
      call dinput(d,10)
      ityp = d(8) 
      if(ityp.le.0) ityp = 1
      if(ityp.gt.2) ityp = 2
      ns  = d(6)
      nsp = d(7)
      if(ns .eq.0) d(6)=  4.d0
      if(nsp.eq.0) d(7)= 10.d0
                   write(iow,1002) typtau(ityp),(d(i),i=1,7),ityp,d(9),
     +                              d(10)
      if(ior.lt.0) write(  *,1002) typtau(ityp),(d(i),i=1,7),ityp,d(9),
     +                              d(10)        
1002  format(5x,'Materialdata for ',a,' in thin walled beams',/,
     + 5x,'E ...  = ',e12.6,/,
     + 5x,'G ...  = ',e12.6,/,
     + 5x,'b ...  = ',e12.6,/
     + 5x,'Qy...  = ',e12.6,/
     + 5x,'Qz...  = ',e12.6,/
     + 5x,'ns...  = ',f4.0,'  only for Q' /
     + 5x,'nps..  = ',f4.0,/    
     + 5x,'ityp.  = ',i4,' (1=Q,2=M_T)',/   
     + 5x,'M_T0.  = ',e12.6,/ 
     + 5x,'Typ I_T= ',e12.6,' (0= only I_Tclosed, 1= I_Tclosed+open)')  
c.... define node numbering for plot mesh routine, see pltord
      inord(ielno)   = 4
      ipord(1,ielno) = 1
      ipord(2,ielno) = 2
      ipord(3,ielno) = 2
      ipord(4,ielno) = 1
c.... description of stresses  
      do i = 1,11         
         forsus(i) =  ' '
      end do
      forsus( 1) =  '  STRESS Tau_s '
      forsus( 2) =  '  FLUX   T_s   '
      if(ityp.eq.1) forsus( 3) =  '  Warping Phi_s'
      if(ityp.eq.2) forsus( 3) =  'Warping Omega_s' ! with respect to M
      fltau = .true.
      area   = 0.d0
      area0  = 0.d0
      az     = 0.d0
      ay     = 0.d0
      azz    = 0.d0
      ayy    = 0.d0
      ayz    = 0.d0
      tm     = 0.d0
      awy    = 0.d0
      awz    = 0.d0
      ait    = 0.d0
      aito   = 0.d0
      return
c.... calculate cross section values
2     dz = xl(2,2)-xl(2,1)
      dy = xl(1,2)-xl(1,1)
      sl = dsqrt(dy*dy + dz*dz)
c                   
      y1  = xl(1,1)    
      z1  = xl(2,1)    
      al  = d(1)*d(3)*sl  
      hl = d(3)

c.... thickness no influence on area0
      if(hl.lt.0) hl=0
      area0= area0 + hl*sl 

      area= area+ al 
      az  = az  + al*(z1+0.5d0*dz)
      ay  = ay  + al*(y1+0.5d0*dy)
      azz = azz + al*(z1*z1+z1*dz+dz*dz/3.d0)
      ayy = ayy + al*(y1*y1+y1*dy+dy*dy/3.d0)
      ayz = ayz + al*(y1*z1+0.5d0*(y1*dz+z1*dy)
     +                +dy*dz/3.d0)

      if(n.eq.numel) then
        ys    = ay/area
        zs    = az/area
        azzs  = azz - zs*zs*area
        ayys  = ayy - ys*ys*area
        ayzs  = ayz - ys*zs*area
c??     dayyzz= azz - ayy
        xad  = 0.5d0*(azzs-ayys)
        if (dabs(ayzs).lt.1.e-5 .and. dabs(xad).lt.1.e-5) then
          phi = 0.0d0
        else
          phi = datan2(-ayzs,xad) 
        end if
        a1i =  0.5d0*(azzs+ayys)
        a2i =  0.5d0*(azzs-ayys)
        aetet  = a1i + a2i*dcos(phi) - ayzs*dsin(phi)
        azeze  = a1i - a2i*dcos(phi) + ayzs*dsin(phi)
        pi     = datan(1.0d0)*4.d0
        phi    = phi*90.d0/pi
                     write(iow,4003)area,az,ay,azz,ayy,ayz,ys,zs,
     +                              azzs,ayys,ayzs,aetet,azeze,phi      
        if(ior.lt.0) write(  *,4003)area,az,ay,azz,ayy,ayz,ys,zs,        
     +                              azzs,ayys,ayzs,aetet,azeze,phi      
c....   copy cross section values:  with respect to center of gravity
        ayy = ayys
        azz = azzs
        ayz = ayzs
        if(ityp.eq.1) then 
c...      basic values
          dn  = azz*ayy-ayz*ayz
          if(ayy.eq.0.d0.and.ayz.eq.0.d0) then
            azq = d(5)/azz
            if(d(4).eq.0.d0) then
              ayq = 0.d0
            else
              stop 'Last Q_z nicht aufnehmbar, I_y=0'
            end if   
          else if(azz.eq.0.d0.and.ayz.eq.0.d0) then
            ayq = d(4)/ayy
            if(d(5).eq.0.d0) then
              azq = 0.d0
            else
              stop 'Last Q_y nicht aufnehmbar, I_z=0'
            end if   
          else
            ayq  = (d(4)*azz-d(5)*ayz)/dn                 
            azq  = (d(5)*ayy-d(4)*ayz)/dn                 
          end if
cww     write(*,*) ayq,azq
        end if
      end if
4003  format(5x,'cross section values in thin walled beams',/,
     + 5x,'n*A ....... ',e14.8,/,
     + 5x,'n*S_y...... ',e14.8,/,
     + 5x,'n*S_z...... ',e14.8,/
     + 5x,'n*I_y...... ',e14.8,/
     + 5x,'n*I_z...... ',e14.8,/
     + 5x,'n*I_yz..... ',e14.8,/
     + 5x,'n*ys....... ',e14.8,/
     + 5x,'n*zs....... ',e14.8,/ 
     + 5x,'n*I_ys..... ',e14.8,/
     + 5x,'n*I_zs..... ',e14.8,/
     + 5x,'n*I_yzs.... ',e14.8,/
     + 5x,'n*I_etas... ',e14.8,/
     + 5x,'n*I_zetas.. ',e14.8,/
     + 5x,'phi........ ',e14.8) 
      return
c
c.... load vector and stiffness matrix
3     dz = xl(2,2)-xl(2,1)
      dy = xl(1,2)-xl(1,1)
      sl = dsqrt(dy*dy + dz*dz)
c.... loadvector
      if(ityp.eq.1) then
        y1q  = xl(1,1)-ys
        z1q  = xl(2,1)-zs
        c3   = -d(1)*sl*sl/2.0d0/d(2)*(ayq*y1q+azq*z1q)
        c4   = -d(1)*sl*sl/6.0d0/d(2)*(ayq*dy +azq*dz) 
c....   thickness 
        hl = abs(d(3))
c
        p(1) = -d(2)*hl/sl*(c3+c4)
        p(2) = -d(2)*hl/sl*(c3+2.0d0*c4)
      else if(ityp.eq.2) then
        xsi_n = (xl(1,1)*dy+xl(2,1)*dz)/sl
        rn1   =  xl(1,1)-xsi_n*dy/sl
        rn2   =  xl(2,1)-xsi_n*dz/sl
        rn    =  dsqrt(rn1*rn1+rn2*rn2)
c....   sign for rn
        vrn   =  (-rn1*dz+rn2*dy)
        vr1   =  1.d0
        vr    =  dsign(vr1,vrn)
        rn    =  rn*vr
c
        p(1) = -d(2)*d(3)*rn        
        p(2) =  d(2)*d(3)*rn               
      end if
c.... stiffness matrix
      f      =  d(2)*hl/sl
      s(1,1) =  f
      s(1,2) = -f
      s(2,1) = -f
      s(2,2) =  f
      return
c.... output stress values etc 
4     dz = xl(2,2)-xl(2,1)
      dy = xl(1,2)-xl(1,1)
      sl = dsqrt(dy*dy + dz*dz) 
      y1q  = xl(1,1)-ys
      z1q  = xl(2,1)-zs
c.... thickness = 0 for connection elements 
      hl = d(3)
      if(hl.lt.0) hl=0
      
      if(ityp.eq.1) then
        c3   = -d(1)*sl*sl/2.0d0/d(2)*(ayq*y1q+azq*z1q)
        c4   = -d(1)*sl*sl/6.0d0/d(2)*(ayq*dy +azq*dz) 
      end if
      if(fltau) then
        if(n.eq.1) then
          cphi   = 0.d0
        end if
c.....  sum values for cphi
        if(ityp.eq.1) then
          cphi=cphi+((ul(1,1)+ul(1,2))/2.0d0-c3/6.0d0-c4/4.0d0)*hl*sl
        else if(ityp.eq.2) then
          c1 = ul(1,1)
          c2 = ul(2,1)-c1
          cphi=cphi+(c1+0.5d0*c2)*d(3)*sl
c.....    and warping moments
          awy=awy+d(1)*d(3)*sl*(y1q*(c1+0.5d0*c2)+dy*(c1/2.d0+c2/3.d0))
          awz=awz+d(1)*d(3)*sl*(z1q*(c1+0.5d0*c2)+dz*(c1/2.d0+c2/3.d0))
c.....    and GIT torsional stiffness
c....     rn
          xsi_n = (xl(1,1)*dy+xl(2,1)*dz)/sl
          rn1   =  xl(1,1)-xsi_n*dy/sl
          rn2   =  xl(2,1)-xsi_n*dz/sl
          rn    =  dsqrt(rn1*rn1+rn2*rn2)
c....     sign for rn
          vrn   =  (-rn1*dz+rn2*dy)
          vr1   =  1.d0
          vr    =  dsign(vr1,vrn)
          rn    =  rn*vr
c....     I_T ST.Venant
          h=d(3) 
          ait   = ait+d(2)*h*sl*(rn-c2/sl)*rn
c....     I_T open           
          aito  = aito+d(2)*sl*h**3/3.d0
        end if
        if(n.eq.numel) then
          fltau = .false.
          cphi  = cphi/area0
cww       write(*,*) cphi
          if(ityp.eq.2) then
            iopen = d(10)
            if(iopen.eq.1) ait = ait+aito 
c.....      print torsional stiffness
                         write(iow,4008) ait
            if(ior.lt.0) write(  *,4008) ait
c.....      calculate center of shear
            dn  = azz*ayy-ayz*ayz
            if(ayy.eq.0.d0.and.ayz.eq.0.d0) then !(Abfragen nicht vollständig geprüft!!ww)
              ym  = -awz/azz
              zm  = zs 
            else if(azz.eq.0.d0.and.ayz.eq.0.d0) then
              ym = ys 
              zm  = awy/ayy
            else
              ym = -(awz*ayy-awy*ayz)/dn
              zm =  (awy*azz-awz*ayz)/dn
            end if
                         write(iow,4006) ym,ym-ys
            if(ior.lt.0) write(*  ,4006) ym,ym-ys
                         write(iow,4007) zm,zm-zs
            if(ior.lt.0) write(*  ,4007) zm,zm-zs
c....       scaling factor theta = M_T /G I_T = verdrillung = twist
            if(ait.ne.0.d0) then ! Absturz, wenn offenes Profil ohne offene Anteile gerechnet wird.
              cmt = (d(9) + d(4)*(zm-zs) - d(5)*(ym-ys))/ait  
                           write(iow,4009) cmt     
              if(ior.lt.0) write(*  ,4009) cmt       
            end if 
          end if
        end if
      else            
c....   print shear stress    
        mct = mct - 1
        if(mct.le.0) then
                       write(iow,4000)o,head,typtau(ityp)
          if(ior.lt.0) write(  *,4000)o,head,typtau(ityp)
          mct = 50
        end if
                     write(iow,4001) n,ma,ix(1),ix(2)
        if(ior.lt.0) write(  *,4001) n,ma,ix(1),ix(2)
c
        if(ityp.eq.1) then
c...      shear flux t = g*b/l*(c2+2*c3*xsi+3*c4*xsi^2)
          c1   =   ul(1,1)
          c2   =   ul(1,2)-c1-c3-c4

          xsi  = 0.d0
          ns   = d(6)
          dxsi = 1.0d0/float(ns) 
          do ii = 1,ns+1
            xsi    = (ii-1)*dxsi
            pxsi   = c1+c2*xsi+c3*xsi**2+c4*xsi**3
            txsi   = d(2)*hl/sl*(c2+2.0d0*c3*xsi+3.0d0*c4*xsi*xsi)
            tauxsi = 0.d0
            if(hl.ne.0.d0) tauxsi = txsi/hl
                         write(iow,4004) xsi,tauxsi,txsi,pxsi-cphi
            if(ior.lt.0) write(  *,4004) xsi,tauxsi,txsi,pxsi-cphi
          end do
c...      extremum
          if(c4.eq.0.d0) goto 41
          xsie = -c3/(3.d0*c4)
          if(xsie.gt.0.d0.and.xsie.lt.1.0d0) then
           txsi   = d(2)*hl/sl*(c2+2.0d0*c3*xsie+3.0d0*c4*xsie*xsie)
           tauxsi = 0.d0
           if(hl.ne.0.d0) tauxsi = txsi/hl
           pxsi   = c1+c2*xsie+c3*xsie**2+c4*xsie**3
                        write(iow,4005) xsie,tauxsi,txsi,pxsi-cphi
           if(ior.lt.0) write(  *,4005) xsie,tauxsi,txsi,pxsi-cphi
          end if
        else if(ityp.eq.2) then
c...      shear flux t = g*b*c2/l = constant.
          xsi_n = (xl(1,1)*dy+xl(2,1)*dz)/sl
          rn1   =  xl(1,1)-xsi_n*dy/sl
          rn2   =  xl(2,1)-xsi_n*dz/sl
          rn    =  dsqrt(rn1*rn1+rn2*rn2)
c....     sign for rn
          vrn   =  (-rn1*dz+rn2*dy)
          vr1   =  1.d0
          vr    =  dsign(vr1,vrn)
          rn    =  rn*vr
          wks   = (ul(1,2) - ul(1,1))/sl
          txsi  = d(2)*d(3)*(wks-rn) * cmt  ! scale by cmt
c.....    open cross section (iopen=1)
          if(abs(txsi).lt.eps1.and.iopen.eq.1) then
            txsi = d(2)*cmt*d(3)*d(3)
          end if 
c
          tauxsi = 0.d0
          if(d(3).ne.0.d0) tauxsi = txsi/d(3)
c.....    omega -> omega bar -> omega tilde           
          pxsi1 = (ul(1,1)-cphi + ym* z1q     - zm* y1q    )*cmt
          pxsi2 = (ul(1,2)-cphi + ym*(z1q+dz) - zm*(y1q+dy))*cmt   
                       write(iow,4004) 0.d0,tauxsi,txsi,pxsi1
          if(ior.lt.0) write(  *,4004) 0.d0,tauxsi,txsi,pxsi1
                       write(iow,4004) 1.d0,tauxsi,txsi,pxsi2
          if(ior.lt.0) write(  *,4004) 1.d0,tauxsi,txsi,pxsi2
        end if
41      continue
        if(ityp.eq.1) then
c....     center of shear
          c5   =  c2+c3+c4
          c6   =  0.5d0*c2+2.0d0/3.0d0*c3+0.75d0*c4
          dtyz =  d(2)*hl*(xl(2,1)*c5 + dz*c6)
          dtzy =  d(2)*hl*(xl(1,1)*c5 + dy*c6)
          tm = tm - dy/sl*dtyz + dz/sl*dtzy
          if(n.eq.numel) then
            qy = d(4)
            qz = d(5)
            ym = 999.d0 
            zm = 999.d0 
            if(qz.ne.0.d0.and.qy.eq.0.d0)then ! print only if one direction is loaded
              ym =  tm/qz
                           write(iow,4006) ym,ym-ys
              if(ior.lt.0) write(*  ,4006) ym,ym-ys
            end if
            if(qy.ne.0.d0.and.qz.eq.0.d0) then
              zm = -tm/qy
                           write(iow,4007) zm,zm-zs
              if(ior.lt.0) write(*  ,4007) zm,zm-zs
            end if
          end if
        end if  
      end if  
c
4000  format(/,a1,20a4,//,2x,'ELEMENT SHEAR STRESSES from ',a,/,
     1 2x,'El',2x,'Mat',1x,
     2 '  Node ','       Tau_s        T_s          Phiq')
4001  format(1x,i4,i4,i4,i4)
4004  format(12x,f5.3,2x,3(1x,e12.6))
4005  format(6x,'max',3x,f5.3,2x,3(1x,e12.6))
4002  format(1x,'Calculate stresses first!')
4006  format(5x,'center of shear:',/,
     + 5x,' y_M  ',e14.8,/, 
     + 5x,'yq_M  ',e14.8)   
4007  format(5x,'center of shear:',/,
     + 5x,' z_M  ',e14.8,/, 
     + 5x,'zq_M  ',e14.8)   
4008  format(5x,'torsional stiffness',/,
     + 5x,'n*I_T...... ',e14.8)   
4009  format(5x,'twist:',/,
     + 5x,' theta = M_T/n* I_T...... ',e14.8)   
      return
c
c.... reaction (flux)
6     dz = xl(2,2)-xl(2,1)
      dy = xl(1,2)-xl(1,1)
      sl = dsqrt(dy*dy + dz*dz)
c.... thickness = 0 for connection elements 
      hl = d(3)
      if(hl.lt.0) hl=0
      if(ityp.eq.1) then
        y1q  = xl(1,1)-ys
        z1q  = xl(2,1)-zs
        c3   = -d(1)*sl*sl/2.0d0/d(2)*(ayq*y1q+azq*z1q)
        c4   = -d(1)*sl*sl/6.0d0/d(2)*(ayq*dy +azq*dz) 
c
        f1   = -d(2)*hl/sl*(c3+c4)
        f2   = -d(2)*hl/sl*(c3+2.0d0*c4)
      else if(ityp.eq.2) then
        xsi_n = (xl(1,1)*dy+xl(2,1)*dz)/sl
        rn1   =  xl(1,1)-xsi_n*dy/sl
        rn2   =  xl(2,1)-xsi_n*dz/sl
        rn    =  dsqrt(rn1*rn1+rn2*rn2)
c....   sign for rn
        vrn   =  (-rn1*dz+rn2*dy)
        vr1   =  1.d0
        vr    =  dsign(vr1,vrn)
        rn    =  rn*vr
c
        f1   = -d(2)*d(3)*rn        
        f2   =  d(2)*d(3)*rn               
      end if
      tq   = d(2)*hl/sl*(ul(1,1)-ul(1,2)) 
      p(1) =  tq - f1
      p(2) = -tq - f2
      return
c
c.... plot values
13    if(nfp.gt.3) return
      ipglsave = ipgl
c.... filled plot 
cww      ipgl = 3   ! plottet not filled
c
c.... thickness = 0 for connection elements 
      hl = d(3)
      if(hl.lt.0) hl=0

      if(fltau.and.nfp.eq.1) then
        if(n.eq.1) write(*,4002)
        goto 131
      end if   
      y1q  = xl(1,1)-ys
      z1q  = xl(2,1)-zs
      dz = xl(2,2)-xl(2,1)
      dy = xl(1,2)-xl(1,1)
      sl = dsqrt(dy*dy + dz*dz)
      sn = dz/sl
      cs = dy/sl
c
      if(ityp.eq.1) then
        c3   = -d(1)*sl*sl/2.0d0/d(2)*(ayq*y1q+azq*z1q)
        c4   = -d(1)*sl*sl/6.0d0/d(2)*(ayq*dy +azq*dz) 
c
c...    shear flux t = g*b/l*(c2+2*c3*xsi+3*c4*xsi^2)
        c1   =   ul(1,1)
        c2   =   ul(1,2)-c1-c3-c4
      else if(ityp.eq.2) then
        iopen = d(10)
c...    shear flux t = g*b*c2/l
        xsi_n = (xl(1,1)*dy+xl(2,1)*dz)/sl
        rn1   =  xl(1,1)-xsi_n*dy/sl
        rn2   =  xl(2,1)-xsi_n*dz/sl
        rn    =  dsqrt(rn1*rn1+rn2*rn2)
c....   sign for rn
        vrn   =  (-rn1*dz+rn2*dy)
        vr1   =  1.d0
        vr    =  dsign(vr1,vrn)
        rn    =  rn*vr
        wks   = (ul(1,2) - ul(1,1))/sl
c...    warping function phi = c1+c2*xsi omega -> omega bar -> omega tilde
        c1    = ul(1,1)-cphi + ym* z1q     - zm* y1q    
        c2    = ul(1,2)-cphi + ym*(z1q+dz) - zm*(y1q+dy) - c1   
      end if 
      nsp  = d(7)
      y0 = xl(1,1)
      z0 = xl(2,1)
      y2 = xl(1,2)
      z2 = xl(2,2)
      dyn= (y2-y0)/nsp 
      dzn= (z2-z0)/nsp
      dxsi = 1.0d0/float(nsp) 
c.... loop over all plot increments
      do ii = 1,nsp
c....   coordinates
        y2 = y0 + ii*dyn
        z2 = z0 + ii*dzn
        y1 = y2 - dyn
        z1 = z2 - dzn
c....   flux and stress
        xsi1 = (ii-1)*dxsi
        xsi2 =  ii   *dxsi
        if(ityp.eq.1) then
          t1   = d(2)*hl/sl*(c2+2.0d0*c3*xsi1+3.0d0*c4*xsi1*xsi1)
          t2   = d(2)*hl/sl*(c2+2.0d0*c3*xsi2+3.0d0*c4*xsi2*xsi2)
          p1   = c1+c2*xsi1+c3*xsi1**2+c4*xsi1**3-cphi
          p2   = c1+c2*xsi2+c3*xsi2**2+c4*xsi2**3-cphi
        else if(ityp.eq.2) then
c.....    closed cross section
          t1   = d(2)*d(3)*(wks-rn)*cmt
c.....    open cross section
          if(abs(t1).lt.eps1.and.iopen.eq.1) then
            t1 = d(2)*cmt*d(3)*d(3)
          end if 
          t2 = t1
c
          p1   = (c1+c2*xsi1)*cmt
          p2   = (c1+c2*xsi2)*cmt
        end if 
        tau1 = 0.d0
        tau2 = 0.d0
        if(hl.ne.0.d0) then
          tau1 = t1/hl
          tau2 = t2/hl
        end if
c....   values for increment
        if(nfp.eq.1) then
          s1 = tau1
          s2 = tau2
        else if(nfp.eq.2) then
          s1 = t1
          s2 = t2
        else if(nfp.eq.3) then
          s1 = p1
          s2 = p2
        end if
        klayf = 1
        if(flfp) then
          ccfp = max(abs(s1),abs(s2))
          ccfp1 = max(s1,s2)
          ccfp2 = min(s1,s2)
          xmaxf = max(xmaxf,ccfp1)
          xminf = min(xminf,ccfp2)
          cfp  = max(cfp,ccfp)
        else
          call pzero(xll,12)
          if(abs(s1).lt.eps.and.abs(s2).lt.eps) goto 131
          if(abs(s1).lt.eps) s1 = 0.d0
          if(abs(s2).lt.eps) s2 = 0.d0
c
          if((s1.ge.0.d0.and.s2.ge.0.d0).or.
     +       (s1.le.0.d0.and.s2.le.0.d0)) then
c....       plotte trapez
            sm = 0.5d0*(s1+s2)
            call pppcolf(sm)
            xll(1,1) = y1     
            xll(2,1) = z1     
            xll(1,2) = y2     
            xll(2,2) = z2     
            xll(1,3) = y2 - sn*s2*cfp
            xll(2,3) = z2 + cs*s2*cfp
            xll(1,4) = y1 - sn*s1*cfp
            xll(2,4) = z1 + cs*s1*cfp
            call plot9s(ixl,xll,3,4)
          else
c.....      Vorzeichenwechsel -> rechne zwischen (0,x1) und (sl,x2) lokal
            sld = dxsi*sl
            x1d = 0.d0
            x2d = sld 
            y1d = s1 
            y2d = s2 
c.....      gerade y = y1d + a1 * x   mit a1 = (y2d-y1d)/sld 
            a1  = (y2d-y1d)/sld 
c.....      Durchstosspunkt
            xd = -y1d/a1
            if(xd.le.x1d.or.xd.ge.x2d) goto 131  ! ausserhalb
            xdp(1) =  y1 * (1.d0-xd/sld) + y2 * xd/sld
            xdp(2) =  z1 * (1.d0-xd/sld) + z2 * xd/sld
c....     plotte dreieck 1
            sm = 0.5d0*s1
            call pppcolf(sm)
c           call pzero(xll,12)
            xll(1,1) = y1     
            xll(2,1) = z1     
            xll(1,2) = xdp(1)
            xll(2,2) = xdp(2)
            xll(1,3) = y1 - sn*s1*cfp
            xll(2,3) = z1 + cs*s1*cfp 
            call plot9s(ixld,xll,3,3)
c....     plotte dreieck 2
            sm = 0.5d0*s2
            call pppcolf(sm)
            call pzero(xll,12)
            xll(1,1) = xdp(1)
            xll(2,1) = xdp(2)
            xll(1,2) = y2     
            xll(2,2) = z2     
            xll(1,3) = y2 - sn*s2*cfp
            xll(2,3) = z2 + cs*s2*cfp 
            call plot9s(ixld,xll,3,3)
          end if
c...      Beginn Schreibe Anfangs-und Endwert
          dx1 = .002/scale
          if(ii.eq.1) then          
            nt = s1            
            yt  = y1 -12*dx1                            
            zt  = z1                                            
            call plotl(yt,zt,0.0d0,3)
            call pppcol(1)
            call plablx(s1)
          else if(ii.eq.nsp) then          
            nt = s2            
            yt  = y1 -12*dx1                      
            zt  = z1                               
            call plotl(yt,zt,0.0d0,3)
            call pppcol(4)
            call plablx(s2)
          end if
c...      End Schreibe Anfangs-und Endwert
        end if
      end do
131   ipgl = ipglsave
      return
      end
c
      subroutine plablx(xn)
c-----------------------------------------------------------------------
c     plot number on screen                                             |
c-----------------------------------------------------------------------
      USE hpgl1
      USE pftn77      
      USE plotter
      implicit double precision (a-h,o-z)
      character yy*8
      write(yy,'(f8.2)') xn
      call draw_text(yy,ixa,iya,icc)
      if(nexte.gt.0.and.iprin.eq.1) then
        call hptext (xxp(2),yyp(2),yy,ihpgl)
      end if
      return
      end
