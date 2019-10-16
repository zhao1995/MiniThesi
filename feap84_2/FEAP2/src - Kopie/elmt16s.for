      subroutine elmt16(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c     3D-Seilelement 
c-----------------------------------------------------------------------
c     Eingabe: 
c     Karte 1: 
c     d( 1) = E
c     d( 2) = A
c     d( 3) = Wichte gamma -> Eigengewicht in -z-Ri=gamma*A bez. auf Länge 
c     d( 4) = Temperaturausdehnungskoeffizient
c     d( 5) = Verhältnis Seillänge/Seilsehnenlänge im Ausgangszustand
c     d( 6) = Momentankofiguration durch Nachziehen des Ausgangszustands 0=ja, 1=nein
c
c     Karte 2: 
c     d(11) = Vorspannkraft(in Seilsehnenrichtung)
c     d(12) = Streckenlast x-global/Projektion
c     d(13) = Streckenlast y-global/Projektion
c     d(14) = Streckenlast z-global/Projektion
c     d(15) = Streckenlast x-lokal/Seillänge
c     d(16) = Streckenlast y-lokal/Seillänge
c     d(17) = Streckenlast z-lokal/Seillänge
c     d(18)= Temperaturänderung
c-----------------------------------------------------------------------
c.....2D Seil+FW       : ndm=2,ndf=2 -> P(4)
c.....2D Seil+FW+Balken: ndm=2,ndf=3 -> P(6) 

c.....3D Seil+FW       : ndm=3,ndf=3 -> P(6)
c.....3D Seil+FW+Balken: ndm=3,ndf=6 -> P(12) 
c-----------------------------------------------------------------------
c     forc,1 Seilkraft
c     [forc,2 lokales KOS in Elementmitte] nicht aktiv, da Konflikt bei
c             Mischung mit Stab-Elementen!!  
c-----------------------------------------------------------------------
c     OFFEN
c     Lagerkräfte aus lokalen Lasten
c
c     Seilkraft:
c     # 1 Newton-Verfahren geht,   Startwerte??
c     # 2 Pegasus, fehlt  
c     # 3 direkte Lösung  
c     # im Moment Vergleich 1+3 und Ausgabe von Differenzen 
c     Abfangen von negativer vorspannung derzeit: c
c
c     Plot FORC auf verformtes Netz: geht nur auf x=X+U, nicht auf x=X+scal*u
c
c     geht das allein 2D?? oder nur 3D?? 
c
c
c-----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE eldata
      USE fornam
      USE iofile
      USE pdata6 
      USE pdata10 
      USE pltran
      USE prlod 
      implicit double precision(a-h,o-z)
      logical zoom
      dimension d(*),ul(ndf,*),xl(ndm,*),yp(3,2),
     +          s(nst,nst),p(nst),tm(3,3),rl(6),pl(6),g1(7),g0(7),
     +          res(2),tl(*),ix(*),ixl(5),xll(3,4),dy(3),yl(3,2)
      dimension h1(*),h2(*),h3(*)
      
      data ixl /1,2,3,4,1/,eps/1.0e-5/
c----------------      
      ielno = 16
c----------------      
c.... go to correct array processor
      go to(1,2,3,3,5,3,7,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2), isw
c.... input material properties
1     if(ior.lt.0) write(*,1002)
1002  format(' Input: E, A,gam,alphaT,s_0/l_0,igeo//,
     *         H_0,qx_g,qy_g,qz_g,qx_l,qy_l,qz_l,dT'/)

      call dinput(d,6)

      call dinput(d(11),8)

                   write(iow,1001) (d(i),i=1,6),(d(i),i=11,18)
      if(ior.lt.0) write(  *,1001) (d(i),i=1,6),(d(i),i=11,18)
1001  format(5x,'materialdata cable element:',/,
     1      5x,'elastic modulus......',e12.5,/,
     2      5x,'area.................',e12.5,/,
     3      5x,'gamma................',e12.5,/,
     4      5x,'alphaT...............',e12.5,/,
     5      5x,'relation.. s_0/l_0...',e12.5,/,
     6      5x,'igeo 0=NZ...1=AZ.....',f4.0,/,
     7      5x,'pre-stress force.....',e12.5,/,
     8      5x,'qx_global............',e12.5,/,
     9      5x,'qy_global............',e12.5,/,     
     +      5x,'qz_global............',e12.5,/,
     1      5x,'qx_lokal.............',e12.5,/,
     2      5x,'qy_lokal.............',e12.5,/,     
     3      5x,'qz_lokal.............',e12.5,/,
     4      5x,'delta T..............',e12.5,/)
      if(d(11).lt.0.d0) 
     + stop 'negative pre-stress force in cable element not possible!'

c.... define node numbering for plot mesh routine, see pltord
      inord(ielno)   = 3
      ipord(1,ielno) = 1
      ipord(2,ielno) = 2
      ipord(3,ielno) = 1

c.... description of stresses  
      do i = 1,11         
        forsus(i) =  ' '
      end do
      forsus( 1) =  'Cable Force'
c     forsus( 2) =  'Local COS'

2     return

c.... stiffness matrix
3     continue
c.... yl: aktuelle Knotenkoordinaten
      igeo = d(6)
      do i=1,nel
        do j=1,ndm
          yl(j,i)=xl(j,i)+ul(j,i)
        end do
      end do

c.... g1: verformte Elementgeometrie für x+u
      g1(1) = yl(1,2) - yl(1,1)                          ! dx
      g1(2) = yl(2,2) - yl(2,1)                          ! dy
      g1(3) = yl(3,2) - yl(3,1)                          ! dz  
      g1(4) = dsqrt(g1(2)*g1(2)+g1(3)*g1(3))             ! dy^2+dz^2	
      g1(5) = dsqrt(g1(1)*g1(1)+g1(3)*g1(3))             ! dx^2+dz^2
      g1(6) = dsqrt(g1(1)*g1(1)+g1(2)*g1(2))             ! dx^2+dy^2
      g1(7) = dsqrt(g1(1)*g1(1)+g1(2)*g1(2)+g1(3)*g1(3)) ! ds

c.....g0:...unverformt Elementgeometrie       
      g0(1) = xl(1,2) - xl(1,1)                          ! dX
      g0(2) = xl(2,2) - xl(2,1)                          ! dY
      g0(3) = xl(3,2) - xl(3,1)                          ! dZ
      g0(4) = dsqrt(g0(2)*g0(2)+g0(3)*g0(3))             ! dy^2+dz^2	
      g0(5) = dsqrt(g0(1)*g0(1)+g0(3)*g0(3))             ! dx^2+dz^2
      g0(6) = dsqrt(g0(1)*g0(1)+g0(2)*g0(2))             ! dx^2+dy^2
      g0(7) = dsqrt(g0(1)*g0(1)+g0(2)*g0(2)+g0(3)*g0(3)) ! dS       
 

c.....Anfangsseillänge       
      s0 = g0(7)*d(5)          

c.... Sehnenlänge
      if(igeo.eq.0) dl = g1(7)  
      if(igeo.eq.1) dl = g0(7)  

c.... verformte Sehnenlänge
      dll=g1(7)
                                                                     
c.... Transformationmatrix                             
      if(igeo.eq.0) call trafo16(tm,g1)                                          
      if(igeo.eq.1) call trafo16(tm,g0)                                          

c.... Lastvektor P(12) global qx,qy,qz-g  immer verformt
      p(1) = prop*0.5d0 * d(12)*g1(4)
      p(2) = prop*0.5d0 * d(13)*g1(5)
      p(3) = 0.d0
      if(ndm.gt.2)
     *p(3)=prop*0.5d0 *(d(14)*g1(6) - d(2)*d(3)*g1(7))
      p(ndf+1) = p(1)
      p(ndf+2) = p(2)
      p(ndf+3) = p(3)


c.... Lastvektor lokal qx,qy,qz immer verformt
      pl(1) = prop*0.5d0 * d(15)*g1(7)
      pl(2) = prop*0.5d0 * d(16)*g1(7)
      pl(3) = 0.d0
      if(ndm.gt.2)
     +pl(3) = prop*0.5d0 * d(17)*g1(7) 
      pl(ndm+1) = pl(1)
      pl(ndm+2) = pl(2)
      pl(ndm+3) = pl(3)
      call trafo16p(pl,tm)

c.....Berechnung der Streckenlast senkrecht zur Seilsehne 
      call qresult(g1,qt,p,pl,ndm,ndf)

c.... Seilkraft immer verformt
      call seilkraft(d,h,qt,s0,dll,n,iow)

      if(isw.eq.4.or.isw.eq.13) goto 4

c.... Rücktransformation für igeo=1 
      if(igeo.eq.1) then
         ukx=(ul(1,2)-ul(1,1))/dl  
         h=h/(1.d0+ukx)
      end if

c.... Innere Kräfte
      call pzero(rl,6)
      rl(1) = -h
      rl(4) =  h
      if(igeo.eq.1) then
        vkx=(ul(2,2)-ul(2,1))/dl  
        wkx=(ul(3,2)-ul(3,1))/dl  
        rl(2) = -h*vkx 
        rl(3) = -h*wkx 
        rl(5) =  h*vkx 
        rl(6) =  h*wkx 
      end if
      
      call trafo16p(rl,tm)

c.... Residuum
      do i=1,ndm  !ndm ???
        p(    i)=p(    i) + pl(i)   - rl(i)
        p(ndf+i)=p(ndf+i) + pl(i+3) - rl(i+3)
      end do
     
      if (isw.eq.6) goto 30

c.... Tang. Steifigkeitsmatrix
      call ktang16(s,d,h,qt,s0,dl,dll,vkx,wkx,ndf,nst,igeo) 
      call trafo16k(s,tm,nst,ndf)
30    continue
      return

c.... Berechnung der Auflagerkräfte
4     rx=-prop*0.5d0*d(12)*g1(4)
      ry=-prop*0.5d0*d(13)*g1(5)
      rz=-prop*0.5d0*(d(14)*g1(6)-d(2)*d(3)*g1(7))

c.... Umrechnung der Seilkraft (2.PK -> Cauchy)
      if(igeo.eq.1) then
         ukx=(ul(1,2)-ul(1,1))/dl  
         h=h*(1.d0+ukx)
      end if

      if(igeo.eq.1) h=h*(1.d0+ukx)

c.... Berechnung der Seilkraft am Anfang/Ende
      res(1)=dsqrt((rx-h*g1(1)/g1(7))*(rx-h*g1(1)/g1(7))
     +               +(ry-h*g1(2)/g1(7))*(ry-h*g1(2)/g1(7))
     +               +(rz-h*g1(3)/g1(7))*(rz-h*g1(3)/g1(7)))
      res(2)=dsqrt((rx+h*g1(1)/g1(7))*(rx+h*g1(1)/g1(7))
     +               +(ry+h*g1(2)/g1(7))*(ry+h*g1(2)/g1(7))
     +               +(rz+h*g1(3)/g1(7))*(rz+h*g1(3)/g1(7)))
      if (isw.eq.13)  goto 13
      mct = mct - 2
c.... Print Results
      if (mct .gt. 0) go to 20
      write (iow,*) 'Cable No.  S_chord     S_node_1    S_node_2'
      write (  *,*) 'Cable No.  S_chord     S_node_1    S_node_2'
      mct = 50
20    write(iow,1000)n,h,res(1),res(2)
      write(  *,1000)n,h,res(1),res(2)
      return

c.... mass matrix 
5     return

c.... surface loads
7     return

c.... plot forces
13    continue
c.... normal vector
c.... compute dir. cosin terms for tr(3,3) and member length
      if (scal.eq.0.d0) then
        dx1 = g0(1)
        dx2 = g0(2)
        dx3 = g0(3)
        dl  = g0(7)
      else
        dx1 = g1(1)
        dx2 = g1(2)
        dx3 = g1(3)
        dl  = g1(7)
        do i=1,nel
          do j=1,ndm
            xl(j,i)=yl(j,i)
          end do
        end do
      end if

c.... unit vector e_xL
      dx1 = dx1/dl
      dx2 = dx2/dl
      dx3 = dx3/dl

c.....Hilfsvektor e_Z
      if(nel.eq.2) then ! 2 Knoten e_z = 0,0,1
        dz1 = 0.d0
        dz2 = 0.d0
        dz3 = 1.d0
      end if

c.....e_yL =  e_xL x e_Z = -e_Z x e_xL
      fac = -1.d0      
      dy1 = fac*(dz2*dx3 - dz3*dx2)
      dy2 = fac*(dz3*dx1 - dz1*dx3)
      dy3 = fac*(dz1*dx2 - dz2*dx1)

c.... length of e_yL
      dn  = dsqrt(dy1*dy1+dy2*dy2+dy3*dy3)

c.... special case e_xL = +-e_Z
      deps = 1.e-7*dl
      if(dn.lt.deps) then
        dn  =  1.d0
        dy1 =  0.d0
        dy2 = -1.d0
        dy3 =  0.d0
      end if

c.... unit vector e_yL
      dy1 = dy1/dn
      dy2 = dy2/dn
      dy3 = dy3/dn

c.... unit vector e_zL = e_xL x e_yL
      dz1 = dx2*dy3 - dx3*dy2
      dz2 = dx3*dy1 - dx1*dy3
      dz3 = dx1*dy2 - dx2*dy1
c
      klayf = 1
      if(abs(nfp).gt.1) return 

      s1  = res(1) ! Kraft am Anfang
      s2  = res(2) ! Kraft am Ende

      if(flfp) then 
        ccfp  = max(abs(s1),abs(s2))
        ccfp1 = max(s1,s2)
        ccfp2 = min(s1,s2)
        xmaxf = max(xmaxf,ccfp1)
        xminf = min(xminf,ccfp2)
        cfp  = max(cfp,ccfp)
      else 
cc....   plot local axis
c        if(nfp.eq.2) then
c          call pppcol(2)
c          call plloco16(tm,xl,g1(7),ndm)
c          return
c        end if
c....   Randwerte
        if(abs(s1).lt.eps) s1 = 0.d0
        if(abs(s2).lt.eps) s2 = 0.d0
        if(abs(s1).lt.eps.and.abs(s2).lt.eps) return

        if(ifor.eq.12) then !x-y
          dp1 =  dy1
          dp2 =  dy2
          dp3 =  dy3
        elseif(ifor.eq.-12) then
          dp1 = -dy1
          dp2 = -dy2
          dp3 = -dy3
        elseif(ifor.eq.13) then !x-z
          dp1 =  dz1
          dp2 =  dz2
          dp3 =  dz3
        elseif(ifor.eq.-13) then
          dp1 = -dz1
          dp2 = -dz2
          dp3 = -dz3
        end if
c.....  schleife fuer alle Schnittgroessen
        nk = 6
c.....  Weginkrement für veränderliche Werte
        do kk = 1,3
          dy(kk) = (xl(kk,2) - xl(kk,1))/float(nk)
        end do
        dl1 =  dl/float(nk)
c
        do 135 ii = 1,nk
c.....  Randkoordinaten
        do kk = 1,3
          yp(kk,1) = xl(kk,1) + dy(kk) *(ii-1)
          yp(kk,2) = xl(kk,1) + dy(kk) * ii
        end do
c.....  Randschnittgrössen
        dla = dl1 * (ii-1)
        dle = dl1 *  ii
        if(nfp.eq.1) then
         rxa=-prop*(0.5d0-dla/g1(7))*d(12)*g1(4)
         rya=-prop*(0.5d0-dla/g1(7))*d(13)*g1(5)
         rza=-prop*(0.5d0-dla/g1(7))*(d(14)*g1(6)-d(2)*d(3)*g1(7))
         rxe=-prop*(0.5d0-dle/g1(7))*d(12)*g1(4)
         rye=-prop*(0.5d0-dle/g1(7))*d(13)*g1(5)
         rze=-prop*(0.5d0-dle/g1(7))*(d(14)*g1(6)-d(2)*d(3)*g1(7))
         resa=dsqrt((rxa-h*g1(1)/g1(7))*(rxa-h*g1(1)/g1(7))
     +               +(rya-h*g1(2)/g1(7))*(rya-h*g1(2)/g1(7))
     +               +(rza-h*g1(3)/g1(7))*(rza-h*g1(3)/g1(7)))
         rese=dsqrt((-rxe+h*g1(1)/g1(7))*(-rxe+h*g1(1)/g1(7))
     +               +(-rye+h*g1(2)/g1(7))*(-rye+h*g1(2)/g1(7))
     +               +(-rze+h*g1(3)/g1(7))*(-rze+h*g1(3)/g1(7)))
         s1  = resa
         s2  = rese
        end if
        if(abs(s1).lt.eps) s1 = 0.d0
        if(abs(s2).lt.eps) s2 = 0.d0
c
c.....  Berechne Farbe                                
        sm = (s1+s2)*0.5d0
        call pppcolf(sm)
c.....  plotte Schnittgroesse          
c....   plotte trapez
        xll(1,1) = yp(1,1)
        xll(2,1) = yp(2,1)
        xll(3,1) = yp(3,1)
        xll(1,2) = yp(1,2)
        xll(2,2) = yp(2,2)
        xll(3,2) = yp(3,2)
        xll(1,3) = yp(1,2) + dp1*s2*cfp
        xll(2,3) = yp(2,2) + dp2*s2*cfp 
        xll(3,3) = yp(3,2) + dp3*s2*cfp
        xll(1,4) = yp(1,1) + dp1*s1*cfp
        xll(2,4) = yp(2,1) + dp2*s1*cfp 
        xll(3,4) = yp(3,1) + dp3*s1*cfp
        call plxtrn(xll,tra,vr,3,4) ! drehe bei rot
        if(zoom(xll,3,4)) call plot9s(ixl,xll,ndm,4) 
135     continue
      end if
      return

1000  format(i5,4(2x,f10.4))
      end
c
      subroutine qresult(g1,qt,p,pl,ndm,ndf)
c-----------------------------------------------------------------------
c     Berechnung der Streckenlast senkrecht zur (verformten) Seilsehne 
c     qt      : Streckenlast senkrecht zur Seilsehne
c     g1      : Geometrie für x+u
c     p(12)   : globaler Lastvektor aus globalen Lasten
c     pl(6)   : globaler Lastvektor aus lokalen  Lasten
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension p(12),pl(6),g1(7),q(3)
      
c.....Resultierende Streckenlast
      q(1) = (p(1)+p(ndf+1)+pl(1)+pl(ndm+1))/g1(7)
      q(2) = (p(2)+p(ndf+2)+pl(2)+pl(ndm+2))/g1(7)
      q(3)  = 0.d0
      if(ndm.gt.2) 
     +q(3) = (p(3)+p(ndf+3)+pl(3)+pl(ndm+3))/g1(7)

      qres=dsqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))

c.....Winkel zwischen Streckenlastvektor und Seilsehnenvektor
      eps=1e-12
      if (abs(qres).lt.eps) then 
        sn=0.d0
      else if (abs(qres).gt.eps) then           
        cs=(q(1)*g1(1)+q(2)*g1(2)+q(3)*g1(3))/(qres*g1(7))
        sn=dsqrt(1-cs*cs)
      end if       

c.....Streckenlast senkrecht zur Seilsehne
      qt=qres*sn
      end
c
      subroutine seilkraft(d,h,qt,s0,dl,n,iow)
c-----------------------------------------------------------------------
c     Berechnung der Seilkraft
c
c     H^3+H^2[EA(1-L/S0+alpt*delT)-H0] - EAqt^2*L^3/(24*S0)=0
c     H^3 + H^2*b - c = 0 
c
c     H=h,h1 - Seilkraft
c     qt     - Streckenlast senkrecht zur Seilsehne
c     s0     - Anfangsseillänge
c     dl     - Sehne
c     n      - Elementnummer
c-----------------------------------------------------------------------
c     iterativ: Newton-Verfahren    zum Test
c     iterativ: Pegasus-Verfahren   nicht programmiert
c     direkt:   Kardanische Formeln sollte finale Lösung sein.
c
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension d(*)

c.....Koeffizienten der Seilgleichung H^3+H^2*b-c=0

      ea =d(1)*d(2) 
      b=ea*(1.d0-dl/s0)+ea*d(4)*d(18)-d(11) 
      c=ea*qt*qt*dl**3/24.d0/s0

c.....Startwert H für Seilsehnenkraft, wieso ww?? 
      ze =1e-12
      if (b.lt.ze) then
        h=-(2.d0/3.d0)*b+10.d0
      else if (b.gt.ze) then   
        h=10.d0
      end if      

c.....1 - Newtonverfahren, alternativ Pegasus möglich WW       
      ic= 0
100   ic=ic+1

      if (ic.gt.100) then
        write(iow,*)'no convergence in cable element 16, el= ',n
        write(  *,*)'no convergence in cable element 16, el= ',n
        stop
      end if

      fs=h*h*h+b*h*h-c
      fsd=3*h*h+2*b*h
      h1=h-(fs/fsd)
      delta=(h-h1)
      h=h1

      if(abs(delta).gt.1e-12) go to 100


c.... 2 - Direkte Lösung: Handbuch der Mathematik Cardanische Formel p 109 WW
c.... Umbennung
      r = b
      t =-c
      s = 0  
c.... Substitution x->y-r/3
      p=s-r**2/3.d0
      q=2.d0/27.d0*r**3-r*s/3.d0+t
c.... Anzahl der Wurzeln
      crit=(q/2.d0)**2+(p/3.d0)**3
c.... Lösung
      pi = 4.d0*datan(1.d0)
      if(crit.ge.0.d0) then
c....   Cardanische Formel       
        cc = sqrt(crit)
        u1 = (-q/2.d0+cc)**(1.d0/3.d0)  
        v1 = (-q/2.d0-cc)**(1.d0/3.d0)  
        y1 = u1+v1
c....   1 reelle + 2 konjugiert komplexe Lösungen (hier=0 gesetzt!)         
        x1 = y1-r/3.d0
        x2 = 0.d0
        x3 = 0.d0
      else
c....   'Casus Irreducibilis'
        rr   = sqrt(-p**3/27.d0) 
        cphi = -q/2.d0/rr 
        if(cphi.gt. 1.d0) cphi= 1.d0
        if(cphi.lt.-1.d0) cphi=-1.d0
        if(dabs(cphi).lt.1.e-12) then ! phi=90
          phi3  = pi/6.d0
        else 
          phi3  = (dacos(cphi))/3.d0
        end if 
        c1   = 2.d0*rr**(1.d0/3.d0)
        c2   = 2.d0*pi/3.d0
        y1   = c1*cos(phi3)
        y2   = c1*cos(phi3+c2)
        y3   = c1*cos(phi3+2.d0*c2)
c....   3 reelle Lösungen (Annahme x1=massgebend!)         
        x1 = y1-r/3.d0
        x2 = y2-r/3.d0
        x3 = y3-r/3.d0
      end if  
c
c.... Ergebnisvergleich (TEST)
      if(dabs(h-x1).gt.1.e-6) then
        write(iow,1000) n,ic,h,x1,x2,x3
        write(  *,1000) n,ic,h,x1,x2,x3
1000    format(1x,'Unterschied Seilkraft>1.e-6 in Element: ',i4,/,
     +  1x,'H aus Newton nach: ',i3,' Iterationen:',f12.5,/,       
     +  1x,'H=X1 direkt, sowie X2,X3       :',3x,3(f12.5,2x))
      end if
      
      end

c
      subroutine ktang16(s,d,h,qt,s0,dl,dll,vkx,wkx,ndf,nst,igeo) 
c-----------------------------------------------------------------------
c     tangentielle Steifigkeitsmatrix 
c     s       - KT = K0 + KG
c     K0      - linearer Anteil
c     KU      - materieller Anteil
c     KG      - geometrische Matrix
c     h       - Seilkraft
c     qt      - Streckenlast senkrecht zur Seilsehne
c     s0      - Seillänge unverformt
c     dl      - Sehne 
c     dll     - Sehne verformt 
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension s(nst,*),d(*)

c.... Materieller Anteil 1:         
      ea = d(1)*d(2) 

      phim1= s0/dll + qt**2*dll**2*ea/12.d0/h**3 ! red. E-Modul

      fm=ea/phim1/dl 

      s(    1,    1) =  fm
      s(    1,ndf+1) = -fm    
      s(ndf+1,ndf+1) =  fm


      if(igeo.eq.1) then
c....   Materieller Anteil 2:         
        vkx2=vkx*vkx  
        wkx2=wkx*wkx  
        
        s(    1,    2) =  fm*vkx
        s(    1,    3) =  fm*wkx
        
        s(    1,ndf+2) = -fm*vkx
        s(    1,ndf+3) = -fm*wkx
        
        s(    2,    2) =  fm*vkx2
        s(    2,    3) =  fm*vkx*wkx 
        
        s(    2,ndf+1) = -fm*vkx
        s(    2,ndf+2) = -fm*vkx2
        s(    2,ndf+3) = -fm*vkx*wkx 
        
        s(    3,    3) =  fm*wkx2
        
        s(    3,ndf+1) = -fm*wkx
        s(    3,ndf+2) = -fm*vkx*wkx
        s(    3,ndf+3) = -fm*wkx2
        
        s(ndf+1,ndf+2) =  fm*vkx 
        s(ndf+1,ndf+3) =  fm*wkx 
        
        s(ndf+2,ndf+2) =  fm*vkx2
        s(ndf+2,ndf+2) =  fm*vkx*wkx
        
        s(ndf+3,ndf+3) =  fm*wkx2

      end if
             

c.... Geometrischer Anteil  
      fg = h/dl


      s(    2,    2) = s(    2,    2) + fg
      s(    3,    3) = s(    3,    3) + fg

      s(    2,ndf+2) = s(    2,ndf+2) - fg
      s(    3,ndf+3) = s(    3,ndf+3) - fg

      s(ndf+2,ndf+2) = s(ndf+2,ndf+2) + fg
      s(ndf+3,ndf+3) = s(ndf+3,ndf+3) + fg


c.....Symmetrie
      do i = 1,nst
        do  j = 1,i
          s(i,j) =s(j,i)
        end do
      end do  
      end
c  
      subroutine trafo16(tm,g)
c-----------------------------------------------------------------------
c     Berechnung der Transformationsmatrix
c-----------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
      dimension g(7),tm(3,3)

c.... compute dir. cosin terms for tm(3,3) and member length
      call pzero(tm,9)

c.... unit vector e_x
      tm(1,1) = g(1)/g(7)
      tm(1,2) = g(2)/g(7)
      tm(1,3) = g(3)/g(7)
      
      dly = dsqrt(tm(1,1)*tm(1,1)+tm(1,2)*tm(1,2))

c.... vector e_y parallel to global x-y-Ebene
      if (abs(dly).gt.1.d-9) then
         tm(2,1) = -tm(1,2)/dly
         tm(2,2) =  tm(1,1)/dly
      else 
         tm(2,1) = 0.d0
         tm(2,2) = 1.d0
      end if
      tm(2,3) = 0.d0
      
c.... vector e_z = e_x x e_y
      tm(3,1) = -tm(1,3)*tm(2,2)
      tm(3,2) =  tm(1,3)*tm(2,1)
      tm(3,3) =  tm(1,1)*tm(2,2) - tm(1,2)*tm(2,1)
      end

      subroutine trafo16k(s,tm,nst,ndf)
c-----------------------------------------------------------------------
c     transform local tang. stiffness matrix
c      to global stiffness matrix: Kg = Tt * Kl * T
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension s(nst,*), tm(3,3), zw(6,6)
      call pzero(zw,6*6)
c.... first multiplication (zw=Tt*Kl)
      do i=1,3
        do j=1,3
          zw(i,j)     = tm(1,i)*s(1,j)         + tm(2,i)*s(2,j) 
     +                + tm(3,i)*s(3,j)
          zw(i,j+3)   = tm(1,i)*s(1,j+ndf)     + tm(2,i)*s(2,j+ndf)
     +                + tm(3,i)*s(3,j+ndf)
          zw(i+3,j)   = tm(1,i)*s(ndf+1,j)     + tm(2,i)*s(ndf+2,j)
     +                + tm(3,i)*s(ndf+3,j)
          zw(i+3,j+3) = tm(1,i)*s(ndf+1,j+ndf) + tm(2,i)*s(ndf+2,j+ndf)
     +                + tm(3,i)*s(ndf+3,j+ndf)
        end do     
      end do   
c.... second multiplication (Kg=a*T)
      do i=1,3
        do j=1,3
          s(i,j)         = zw(i,1)*tm(1,j)   + zw(i,2)*tm(2,j)
     +                   + zw(i,3)*tm(3,j)
          s(i,j+ndf)     = zw(i,4)*tm(1,j)   + zw(i,5)*tm(2,j)
     +                   + zw(i,6)*tm(3,j)
          s(i+ndf,j)     = zw(i+3,1)*tm(1,j) + zw(i+3,2)*tm(2,j)
     +                   + zw(i+3,3)*tm(3,j)
          s(i+ndf,j+ndf) = zw(i+3,4)*tm(1,j) + zw(i+3,5)*tm(2,j)
     +                   + zw(i+3,6)*tm(3,j)
        end do    
      end do   

      return
      end

      subroutine trafo16p(rl,tm)
c-----------------------------------------------------------------------
c     transform local vector P to global vector: Pg = Tt * P
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension rl(6), tm(3,3) ,rz(6)
      do i=1,3
        rz(i)   = tm(1,i)*rl(1) + tm(2,i)*rl(2) + tm(3,i)*rl(3)
        rz(i+3) = tm(1,i)*rl(4) + tm(2,i)*rl(5) + tm(3,i)*rl(6)
      end do  

      do i=1,6
         rl(i) = rz(i)
      end do

      return
      end
c
      subroutine trafo16u(vl,ul,tm,ndf)
c-----------------------------------------------------------------------
c     transform global to local displacement vector v =T*u
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension ul(ndf,*), vl(3,2), tm(3,3)
      call pzero(vl,6)
      do l = 1,2
        do i= 1,3
          do k = 1,3
            vl(i,k) = vl(i,k) + tm(i,k)*ul(k,l)
          end do  
        end do  
      end do  

      return
      end
c
      subroutine plloco16(tr,xl,dl,ndm)
c--------------------------------------------------------+
c.... Plot local basis
c--------------------------------------------------------+
      implicit double precision (a-h,o-z)
      dimension tr(3,3),x0(3),xl(ndm,*),d1(3),d2(3)
c...  plot at midpoint of element (global coordinates!)
      x0(1) = 0.5*(xl(1,1) + xl(1,2))
      x0(2) = 0.5*(xl(2,1) + xl(2,2))
      x0(3)=0.d0
      if(ndm.gt.2) 
     +x0(3) = 0.5*(xl(3,1) + xl(3,2))
c.... length for plot
      xm = 0.25d0 * dl
c...  plot axis 
      call pltaxs16(tr,x0,xm)
      return
      end
c
      subroutine pltaxs16(tr,x0,xm)
c----------------------------------------------------------------------
c.... draw vectors for axes, including rot-macro 
c.... t_1=black
c.... t_2=red
c.... t_3=blue
c----------------------------------------------------------------------
      USE pdata2
      USE pltran
      implicit double precision (a-h,o-z)
      logical zoom
      dimension xx(3,5),x0(3),x0z(3),tr(3,3)
c.... test if vectors are in plot-region
      x0z(1) = x0(1)
      x0z(2) = x0(2)
      x0z(3) = x0(3)
      call plxtrn(x0z,tra,vr,3,1)
      if(zoom(x0z,3,1)) then
c....   perspective projecton of axes
        do 120 m = 1,3  
	        call pppcol(m)
          call pzero(xx,15)
          do 100 n = 1,3  
            fac1    = tr(m,n)*xm
            xx(n,1) = x0(n)
            xx(n,2) = xx(n,1) + fac1
100       xx(n,5) = xx(n,2)
          fac1 = tr(m,1)*xm
          fac2 = tr(m,2)*xm
          fac3 = tr(m,3)*xm
          xx(1,3) = xx(1,2) -.3*fac1 - .1*(fac2+fac3)
          xx(2,3) = xx(2,2) -.3*fac2 + .1*(fac1+fac3)
          xx(3,3) = xx(3,2) -.3*fac3 + .1*(fac1+fac2)
          xx(1,4) = xx(1,2) -.3*fac1 + .1*(fac2+fac3)
          xx(2,4) = xx(2,2) -.3*fac2 - .1*(fac1+fac3)
          xx(3,4) = xx(3,2) -.3*fac3 - .1*(fac1+fac2)
c.......  transform vector if necessary
          call plxtrn(xx,tra,vr,3,5)
c....     plot the vector
          call plotl(xx(1,1),xx(2,1),xx(3,1),3)
          call plotl(xx(1,2),xx(2,2),xx(3,2),2)
          call plotl(xx(1,2),xx(2,2),xx(3,2),ipgl)
          call plotl(xx(1,3),xx(2,3),xx(3,3),2)
          call plotl(xx(1,4),xx(2,4),xx(3,4),2)
          call plotl(xx(1,2),xx(2,2),xx(3,2),2)
          if(ipgl.eq.1) call clpan
          call plotl(xx(1,2),xx(2,2),xx(3,2),3)
          call plabl(m)
120     continue
      end if
      return
      end
c
