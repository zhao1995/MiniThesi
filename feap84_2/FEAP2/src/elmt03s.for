      subroutine elmt03(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-------------------------------------------------------+
c     lineares 2D-Stabelement                           |
c-------------------------------------------------------+
c
c   * Trapezlast in  n,q                                |
c   * Temperaturlast (delta ueber Hoehe)                |
c   * elastische Bettung                                |
c   * Massenmatrix (lumped/konsistent)                  |
c   * Loads can be defined via qloa                     |
c   * Theorie 2. Ordnung                                |
c   * Stabilitaet (K+lambda*K_G)*Phi = 0                |
c   * Plot Schnittgroessen  FORC,...                    |
c     1 = N, 2 = Q, 3 = M, 4 = Press(c)                 |
c   * el.no 7 used for plot, see inord,ipord            |
c-------------------------------------------------------+
c      Freiheitsgrade  u,w,phi in Richtung x,y,z        |
c      y,w   phi   F_y   M_z                            |
c      |   <-|      |  <-|                              |
c      . ->x,u      .-> F_x                             |
c-------------------------------------------------------+
c       Schnittgroessen (lokal)  nach Statik            |
c           Q                            Q              |
c     |->   ^                              <-| M        |
c     | M   | -> N ---------------  ->N  |   |          |
c                                        v              |
c-------------------------------------------------------+
c 1.Karte d(1)  = E                                     |
c           d(2)        = A                             |
c           d(3)        = I                             |
c           d(4)        = h                             |
c           d(5)        = q_1                           |
c           d(6)        = q_2                           |
c           d(7)        = n_1                           |
c           d(8)        = n_2                           |
c           d(9)        = 0/1/2 Theorie 2. Ordnung      |
c                      0=nein                           |
c                      1=ja, normal                     |
c                      2=ja  rechte Seite               |
c 2.Karte d(10)   = rho                                 |
c           d(11)   = c                                 |
c           d(12)   = b                                 |
c           d(13)   = a_t                               |
c           d(14)   = T_N  old... to                    |
c           d(15)   = T_M  old ...tu                    |
c           d(16)   = 0/1 Schnittgr. ohne/mit Zahlen    |
c             
c d-feld    d(17)                                       |
c-------------------------------------------------------+
c     Loads from QLOA
c
c      q01 = qy1  
c      q02 = qy2  
c      q03 = qx1  
c      q04 = qx2  
c      keine Temperatur
c
c-------------------------------------------------------+
c offen: Lasten global/lokal                            |
c        Spannungen plotten                             |
c        Spannungen > zul Sigma plotten                 |
c        Einheiten  spezialisieren                      |
c        Profildatei                                    |
c        Schubspannungsberechnung ist nur 1.5 Q/A       |
c        schnittgroessen und spannungen zwischen Knoten |
c        forc,-i fehlt                                  |
c-------------------------------------------------------+
      USE bdata
      USE cdata
      USE eldata
      USE evdata
      USE fornam
      USE iofile
      USE pdata6
      USE pdata7
      USE pdata10
      USE prlod
      USE qload
      USE strnam
      implicit double precision (a-h,o-z)
      dimension  xll(3,4),xl(ndm,*),tl(*),ix(*),ixl(5),ixl1(4),
     1         d(*),ul(ndf,*),s(nst,nst),p(nst),vl(3,2),sig(8),
     2         press(2),str(6),smat(nst,nst)
      dimension h1(*),h2(*),h3(*)
      character comp*25   
      data ixl/1,2,3,4,1/,ixl1/1,2,3,1/
      ielno = 3
c.... Sprung zu gewuenschtem Programmteil
      go to (1,2,3,3,5,3,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,3), isw
      return
      
1     if(ior.lt.0) write(*,3000)
      call dinput(d,9)
      if(ior.lt.0) write(*,3001)
      call dinput(d(10),7)

c.... elast. foundation always/only for compression = 0/1
      if(d(11).ge.0.d0)then
         comp='(tension/compression)'
      else
         comp='(  only  compression)'
      end if   
c
                   write(iow,2000) (d(i),i=1,11),comp,(d(i),i=12,16)
      if(ior.lt.0) write(*  ,2000) (d(i),i=1,11),comp,(d(i),i=12,16)

c.... drehe lasten q, damit sie lokal wirken       
      d(5) = -d(5)
      d(6) = -d(6)

c.... define node numbering for plot mesh routine, see pltord
      inord(ielno)   = 4
      ipord(1,ielno) = 1
      ipord(2,ielno) = 2
      ipord(3,ielno) = 2
      ipord(4,ielno) = 1
      ipla = 2
c.... description of stresses  
      forsus( 1) =  'N-FORCE N_x/S_x'
      forsus( 2) =  'Q-FORCE Q_z/R_z'
      forsus( 3) =  '  MOMENT  M_y  '
      forsus( 4) =  'FOUND. PRESSURE'
      do i = 5,11         
         forsus(i) =  ' '
      end do

2     return

c.... Steifigkeitsmatrix
3     cb = d(11)

      cs = xl(1,2) - xl(1,1)
      sn = xl(2,2) - xl(2,1)


      dl = dsqrt(cs*cs + sn*sn)
      cs = cs/dl
      sn = sn/dl
      
    
c.... Transformation der Verschiebungen auf lokale Richtungen
      call tranb7(cs,sn,ul,vl,3,ndf,1)


c....   Trapez--Lastvektor
        call qload03(qx1,qx2,qy1,qy2,d,aqloa,propq,prop,numel,n,
     +               mqloa,isw)

        p(1) = dl/6.d0 * (2.d0*qx1 +      qx2 )
        p(4) = dl/6.d0 * (     qx1 + 2.d0*qx2 )
        p(2) = 0.05d0*dl*(7.d0*qy1 + 3.d0*qy2 )
        p(5) = 0.05d0*dl*(3.d0*qy1 + 7.d0*qy2 )
        p(3) = 0.05d0*dl*(           qy1 + 2.d0/3.d0*qy2 )*dl 
        p(6) = 0.05d0*dl*(-2.d0/3.d0*qy1 -           qy2 )*dl 
c
      sigm = 0.d0 
      if(d(13).gt.0.0) then
c....     Temperaturlast fuer t = tm + tb/h
cww          tm = 0.5d0*(d(15)+d(14))
cww          tb =       (d(15)-d(14))
          tm = d(14)
          tb = d(15)
c....     Normalkraft N_t = sigm, Biegemoment M_t = sigb
c         Biegevorzeichen gedreht
          sigm = d(1)*d(2)*d(13)*tm
          sigb = d(1)*d(3)*d(13)*tb/d(4)
          p(1) = p(1) - sigm*prop
          p(4) = p(4) + sigm*prop
          p(3) = p(3) - sigb*prop
          p(6) = p(6) + sigb*prop
      end if
  
      if(isw.eq.22) go to 33
      
c.... Steifigkeitsmatrix
        eal = d(1)*d(2)/dl
        eil = d(1)*d(3)/(dl*dl*dl)
        s(1,1) =  eal
        s(1,4) = -eal
        s(4,4) =  eal
        s(2,2) =  12.d0* eil
        s(2,3) =   6.d0* eil * dl
        s(2,5) = -12.d0* eil
        s(2,6) =   6.d0* eil * dl
        s(3,3) =   4.d0* eil * dl * dl
        s(3,5) =  -6.d0* eil * dl
        s(3,6) =   2.d0* eil * dl * dl
        s(5,5) =  12.d0* eil
        s(5,6) =  -6.d0* eil * dl
        s(6,6) =   4.d0* eil * dl * dl

c.... Anteil K elastische Bettung
      if(dabs(cb).gt.0.0) then
c....   add term always/only for negative values
        press(1) = dabs(cb)*vl(2,1)
        press(2) = dabs(cb)*vl(2,2)
        pressm=press(1)+press(2)
        if(cb.gt.0.d0.or.(cb.lt.0.d0.and.pressm.le.0.0d0)) then
          scb     = dabs(cb)*d(12)*dl/420.d0
          s(2,2) = s(2,2) + scb * 156.d0
          s(2,3) = s(2,3) + scb *  22.d0 * dl
          s(2,5) = s(2,5) + scb *  54.d0
          s(2,6) = s(2,6) - scb *  13.d0 * dl
          s(3,3) = s(3,3) + scb *   4.d0 * dl * dl
          s(3,5) = s(3,5) + scb *  13.d0 * dl
          s(3,6) = s(3,6) - scb *   3.d0 * dl * dl
          s(5,5) = s(5,5) + scb * 156.d0
          s(5,6) = s(5,6) - scb *  22.d0 * dl
          s(6,6) = s(6,6) + scb *   4.d0 * dl * dl
        end if
      end if
c.... Anteil K_G Theorie 2. Ordnung
      ith2 = d(9)

      if(ith2.eq.2) then ! Iteration nur Residuum
        smat = 0
        smat = s 
      end if  
        
      if(ith2.gt.0) then
         xn = eal * (vl(1,2) - vl(1,1) ) /30.d0/dl ! aus Verschiebung 
         xn = xn - sigm*prop /30.d0/dl             ! aus Temperatur 
          s(2,2) = s(2,2) + xn * 36.d0
          s(2,3) = s(2,3) + xn *  3.d0 * dl
          s(2,5) = s(2,5) - xn * 36.d0
          s(2,6) = s(2,6) + xn *  3.d0 * dl
          s(3,3) = s(3,3) + xn *  4.d0 * dl * dl
          s(3,5) = s(3,5) - xn *  3.d0 * dl
          s(3,6) = s(3,6) - xn *         dl * dl
          s(5,5) = s(5,5) + xn * 36.d0
          s(5,6) = s(5,6) - xn *  3.d0 * dl
          s(6,6) = s(6,6) + xn *  4.d0 * dl * dl

      end if

c....   Symmetrie
        s(4,1) = s(1,4)
        s(3,2) = s(2,3)
        s(5,2) = s(2,5)
        s(6,2) = s(2,6)
        s(5,3) = s(3,5)
        s(6,3) = s(3,6)
        s(6,5) = s(5,6)

        if(ith2.eq.2) then
          smat(4,1) = smat(1,4)
          smat(3,2) = smat(2,3)
          smat(5,2) = smat(2,5)
          smat(6,2) = smat(2,6)
          smat(5,3) = smat(3,5)
          smat(6,3) = smat(3,6)
          smat(6,5) = smat(5,6)
        end if

c     Residuum R = P - (K_e+K_G) * v_e 
      do 32 i = 1,6
         do 32 k = 1,3
            p(i) = p(i) - s(i,k)*vl(k,1) - s(i,k+3)*vl(k,2)
32    continue

      if(ith2.eq.2) s = smat ! Iteration nur Residuum

      
      if(isw.eq.4.or.isw.eq.13) goto 4
c....   Transformation von K und P auf globale Richtungen
33      call tranb7(cs,sn,p,s,nst,ndf,2)
        return
        
4     inum = d(16)
c.... Berechnung der Schnittgroessen lokal x-z aus Knotenkraeften  S = - R = -(K*v-P) für  x-y
      sig(1) =  p(1)
      sig(2) = -p(2)
      sig(3) =  p(3)
      sig(4) = -p(4)
      sig(5) =  p(5)
      sig(6) = -p(6)

c.... Pressung bei elastischer Bettung
      press(1) = dabs(cb)*vl(2,1)
      press(2) = dabs(cb)*vl(2,2)
      if(cb.lt.0.d0) then
        if(press(1).gt.0.d0) press(1) = 0.d0
        if(press(2).gt.0.d0) press(2) = 0.d0
      end if
c
      if(isw.eq.4) then
        mct = mct - 3
        if(mct.lt.0) then
          if(ith2.eq.0) then
            if(d(4).gt.0.0d0) then
c...          mit  Spannungen
              if(ior.lt.0) write(*  ,2005) head
                           write(iow,2005) head
            else
c...          ohne Spannungen
              if(ior.lt.0) write(*  ,2006) head
                           write(iow,2006) head
            end if
          else
c...        Th. 2. Ordnung
            if(ior.lt.0) write(*  ,2007) head
                         write(iow,2007) head
          end if
          mct = 50
        end if

        if(ith2.eq.0 .and. d(4).gt.0.0d0) then
c.....    Spannungsberechnung   S = N/A +- M/W   T = 1.5 Q/A
          str(1)  = sig(1)/d(2) + sig(3)/d(3)*0.5d0*d(4)
          str(2)  = sig(1)/d(2) - sig(3)/d(3)*0.5d0*d(4)
          str(3) = 1.5d0*sig(2)/d(2)
          str(4) = sig(4)/d(2) + sig(6)/d(3)*0.5d0*d(4)
          str(5) = sig(4)/d(2) - sig(6)/d(3)*0.5d0*d(4)
          str(6) = 1.5d0*sig(5)/d(2)
          if(ior.lt.0) write(*  ,2002) 
     1    n,ma,' L ',(sig(i),i=1,3),press(1),(str(i),i=1,3),
     2         ' R ',(sig(i),i=4,6),press(2),(str(i),i=4,6)
                         write(iow,2002) 
     1    n,ma,' L ',(sig(i),i=1,3),press(1),(str(i),i=1,3),
     2         ' R ',(sig(i),i=4,6),press(2),(str(i),i=4,6)
        end if

        if(ith2.eq.0 .and. d(4).eq.0.0d0) then
          if(ior.lt.0) write(*  ,2003) 
     1    n,ma,' L ',(sig(i),i=1,3),press(1),
     2         ' R ',(sig(i),i=4,6),press(2)
                         write(iow,2003) 
     1    n,ma,' L ',(sig(i),i=1,3),press(1),
     2         ' R ',(sig(i),i=4,6),press(2)
        end if

        if(ith2.ne.0) then
c.....    Berechnung Stabkennzahl eps = L*sqrt(S/EI)
          sl = dabs(sig(1))
          sr = dabs(sig(4))
          ei = d(1)*d(3)
          epsl = dl*dsqrt(sl/ei)
          epsr = dl*dsqrt(sr/ei)
          epsm = max(epsl,epsr)
          if(ior.lt.0) write(*  ,2004) 
     1    n,ma,' L ',(sig(i),i=1,3),press(1),epsl,
     2         ' R ',(sig(i),i=4,6),press(2),epsr
                         write(iow,2004) 
     1    n,ma,' L ',(sig(i),i=1,3),press(1),epsl,
     2         ' R ',(sig(i),i=4,6),press(2),epsr
          if(epsm.gt.2.5d0) then
            if(ior.lt.0) write(*  ,2008) n
                         write(iow,2008) n
          end if
        end if

      else if(isw.eq.13) then

        if(iplma(ma).eq.0)       return ! only if MATN
        if(nfp.lt.1.or.nfp.gt.4) return
c....   Plotte  N,Q,M,press 
        nd2 = ndf+nfp
        if(nfp.eq.4) then
          if(dabs(d(11)).eq.0.0) return
          sig(4) = press(1)
          sig(8) = press(2)
        end if 
        klayf= 1
        if(flfp) then
          ccfp = max(abs(sig(nfp)),abs(sig(nd2)))
          ccfp1 = max(sig(nfp),sig(nd2))
          ccfp2 = min(sig(nfp),sig(nd2))
          xmaxf = max(xmaxf,ccfp1)
          xminf = min(xminf,ccfp2)
          cfp  = max(cfp,ccfp)
        else
          s1 = sig(nfp)
          s2 = sig(nd2)
          fac = -0.5
          if(nfp.eq.4) then 
            fac = 1.3
            if(dabs(d(11)).eq.0.0) return 
          end if 
c.....    Abbruch bei kleinen Werten
          sq = dsqrt(s1*s1+s2*s2)
          if(sq.lt.1.e-8) return
          if(abs(s1).lt.1.e-8) s1 = 0.0d0
          if(abs(s2).lt.1.e-8) s2 = 0.0d0
c.....    Randkkordinaten
          x1 = xl(1,1)
          y1 = xl(2,1)
          x2 = xl(1,2)
          y2 = xl(2,2)
c......   Lasten
          call qload03(qx1,qx2,qy1,qy2,d,aqloa,propq,prop,numel,n,
     +               mqloa,isw)
          q1 = qy1
          q2 = qy2
          p1 = qx1
          p2 = qx2
c.....    Auswertung fuer  N,Q,M, inkrementweise Berechnung ns = 6
          ns = 6
          x0 = x1
          y0 = y1
          r0 = 0.
          dx = (x2-x1)/ns
          dy = (y2-y1)/ns
          dr = sqrt(dx*dx+dy*dy)
          dr0 = dr*ns
c.....    Schleife ueber alle Plotinkremente
          do 135 is = 1,ns
c.....      Randkoordinaten 1,2
            x2 = x0 + is*dx
            y2 = y0 + is*dy
            r2 = r0 + is*dr
            x1 = x2 - dx
            y1 = y2 - dy
            r1 = r2 - dr
c.....      Randschnittgroessen an Stelle 1,2
c.....      dn=p             
            dp1 = (p2-p1)*r1/dr0
            dp2 = (p2-p1)*r2/dr0
c.....      dq-press*b             
            ps1= press(1) 
            ps2= press(2)
            ps1b=ps1*d(12)
            ps2b=ps2*d(12)
            dq1 = (q2-ps2b-q1+ps1b)*r1/dr0 
            dq2 = (q2-ps2b-q1+ps1b)*r2/dr0 
c.....      N
            sn1 = sig(1) - p1*r1 - dp1*r1/2.  
            sn2 = sig(1) - p1*r2 - dp2*r2/2. 
c.....      Q
            sq1 = sig(2) + (q1-ps1b)*r1 + dq1*r1/2.  
            sq2 = sig(2) + (q1-ps1b)*r2 + dq2*r2/2. 
c.....      M            
            sm1 = sig(3) + sig(2)*r1 + (q1-ps1b)*r1*r1/2. + dq1*r1*r1/6. 
            sm2 = sig(3) + sig(2)*r2 + (q1-ps1b)*r2*r2/2. + dq2*r2*r2/6.
c.....      pressure             
            ps1g= ps1+(ps2-ps1)*r1/dr0
            ps2g= ps1+(ps2-ps1)*r2/dr0
c
            if(nfp.eq.1) then
              s1 = sn1
              s2 = sn2
            else if(nfp.eq.2) then
              s1 = sq1
              s2 = sq2
            else if(nfp.eq.3) then
              s1 = sm1
              s2 = sm2
            else if(nfp.eq.4) then
              s1 = ps1g
              s2 = ps2g
            end if
            s1   = s1*cfp
            s2   = s2*cfp
c.....      Punkte 3,4 (für Trapez 1,2,3,4)
            x3   = x2 + sn*s2
            y3   = y2 - cs*s2
            x4   = x1 + sn*s1
            y4   = y1 - cs*s1
c......     Zahlen
            if(is.eq.1.and.inum.eq.1) then
              x6 = x1 + sn*s1*fac
              y6 = y1 - cs*s1*fac
              call plotl(x6,y6,0.0d0,3)
              call pppcol(1)
              ns1 = s1/cfp
              call plabl(ns1)
            end if
            if(is.eq.ns.and.inum.eq.1) then
              x5 = x2 + sn*s2*fac
              y5 = y2 - cs*s2*fac
              call plotl(x5,y5,0.0d0,3)
              call pppcol(1)
              ns2 = s2/cfp
              call plabl(ns2)
            end if
c.....      Berechne Farbe                                
            sm = 0.5d0*(s1+s2)/cfp
            call pppcolf(sm)
c.....      plotte Schnittgroesse 
            if(s1.ge.0.0.and.s2.ge.0.0) then
              goto 136
            else if(s1.le.0.0.and.s2.le.0.0) then
              goto 136
            else
c.....        Vorzeichenwechsel -> Durchstosspunkt D + sonderfall x1=x2
              if(dabs(x1-x2).gt.1.e-8) then
                a1 =(y1-y2)/(x1-x2)
                b1 =y2-a1*x2
                a2=(y3-y4)/(x3-x4)
                b2=y4-a2*x4
                xd=(b2-b1)/(a1-a2)
              else
                xd=x1
                a2=(y3-y4)/(x3-x4)
                b2=y4-a2*x4
              end if
              yd=a2*xd+b2
c....         test ob D zwischen 1 und 2
c....         berechen  Vektoren P1 = D-1 und P2 = 2-D
c....         P1*P2 negativ -> D ausserhalb, dann nicht plotten
              p1x = xd - x1
              p1y = yd - y1
              p1l = sqrt(p1x*p1x+p1y*p1y)
              p2x = x2 - xd
              p2y = y2 - yd
              p2l = sqrt(p2x*p2x+p2y*p2y)
c.....        nicht moeglich
              if(p1l.eq.0.0.and.p2l.eq.0.0) goto 135
c.....        P1 = 0 -> D = 1  P2 = 0 -> D = 2
              if(p1l.eq.0.0.or.p2l.eq.0) goto 138
              cosphi = (p1x*p2x + p1y*p2y)/p1l/p2l
              if(cosphi.lt.0.0) goto 135
c....         plotte dreieck 1
c.....        Berechne Farbe                                
138           sm = 0.5d0*s1/cfp
              call pppcolf(sm)
              call pzero(xll,12)
              xll(1,1) = x1
              xll(2,1) = y1
              xll(1,2) = xd
              xll(2,2) = yd
              xll(1,3) = x4
              xll(2,3) = y4
              call plot9s(ixl1,xll,3,3)
c....         plotte dreieck 2
c.....        Berechne Farbe                                
              sm = 0.5d0*s1/cfp
              call pppcolf(sm)
              call pzero(xll,12)
              xll(1,1) = xd
              xll(2,1) = yd
              xll(1,2) = x2
              xll(2,2) = y2
              xll(1,3) = x3
              xll(2,3) = y3
              call plot9s(ixl1,xll,3,3)
              goto 135
            end if
c....       plotte trapez
136         call pzero(xll,12)
            xll(1,1) = x1
            xll(2,1) = y1
            xll(1,2) = x2
            xll(2,2) = y2
            xll(1,3) = x3
            xll(2,3) = y3
            xll(1,4) = x4
            xll(2,4) = y4
            call plot9s(ixl,xll,3,4)
135       continue
        end if
      end if
      return
c------------------------------------
c.... Massenmatrizen und Beulmatrizen
5     cs = xl(1,2) - xl(1,1)
      sn = xl(2,2) - xl(2,1)
      dl = dsqrt(cs*cs + sn*sn)
      cs = cs/dl
      sn = sn/dl
      if(imtyp.eq.1) then
c....   lumped Massenmatrix s. unten
c....   konsistente Massenmatrix
        cm  = d(10)*d(2)*dl /420.d0
          s(1,1) =   cm * 140.d0
          s(1,4) =   cm *  70.d0
          s(4,4) =   cm * 140.d0
          s(2,2) =   cm * 156.d0
          s(2,3) =   cm *  22.d0* dl
          s(2,5) =   cm *  54.d0
          s(2,6) =  -cm *  13.d0* dl
          s(3,3) =   cm *   4.d0* dl * dl
          s(3,5) =   cm *  13.d0* dl
          s(3,6) =  -cm *   3.d0* dl * dl
          s(5,5) =   cm * 156.d0
          s(5,6) =  -cm *  22.d0* dl
          s(6,6) =   cm *   4.d0* dl * dl
      else if(imtyp.ne. 1) then
c....   Beulmatrix
c....     Transformation der Verschiebungen auf lokale Richtungen
          call tranb7(cs,sn,ul,vl,3,ndf,1)
c....     Normalkraft/-30l
          xn = -d(1)*d(2)/dl * (vl(1,2) - vl(1,1) ) /30.d0/dl
          s(2,2) =    xn *  36.d0
          s(2,3) =    xn *   3.d0 * dl
          s(2,5) =  - xn *  36.d0
          s(2,6) =    xn *   3.d0 * dl
          s(3,3) =    xn *   4.d0 * dl * dl
          s(3,5) =  - xn *   3.d0 * dl
          s(3,6) =  - xn *        dl * dl
          s(5,5) =    xn *  36.d0
          s(5,6) =  - xn *   3.d0 * dl
          s(6,6) =    xn *   4.d0 * dl * dl
      end if
c....   Symmetrie
        s(4,1) = s(1,4)
        s(3,2) = s(2,3)
        s(5,2) = s(2,5)
        s(6,2) = s(2,6)
        s(5,3) = s(3,5)
        s(6,3) = s(3,6)
        s(6,5) = s(5,6)
c....   Transformation auf globale Richtungen
        call tranb7(cs,sn,p,s,nst,ndf,2)
      if(imtyp.eq.1) then
c....   lumped Massenmatrix ohne Transformation, sonst evtl. neg. Massen
          p(1) = d(10)*d(2)*0.5d0*dl
          p(2) = p(1)
          p(3) = p(1)*d(3)/d(2)
          p(4) = p(1)
          p(5) = p(2)
          p(6) = p(3)
      end if
      return
c.... format statements
2000  format(5x,'Bernoulli Stabelement'//
     +   10x,'E ',e14.5/10x,'A ',e14.5/10x,'I ',e14.5/
     +   10x,'h ',e14.5/10x,'q1',e14.5/10x,'q2',e14.5/
     +   10x,'n1',e14.5/10x,'n2',e14.5/
     +   10x,'Theor. 2.Ordnung 0=nein/1=ja',e14.5/
     +   10x,'rho = gamma/g',e14.5/
     +   10x,'Bettungsziffer',e14.5,a/
     +   10x,'Breite  ',e14.5/
     +   10x,'alpha_T ',e14.5/
     +   10x,'Temp. to',e14.5/10x,'Temp. tu',e14.5/
     +   10x,'Zahlen f. N,M,Q  0=nein/1=ja',e14.5)
2002  format(2i3,a3,7(1x,g9.3),/,6x,a3,7(1x,g9.3))
2003  format(2i3,a3,4(1x,g10.4),/,6x,a3,4(1x,g10.4))
2004  format(2i3,a3,5(1x,g10.4),/,6x,a3,5(1x,g10.4))
2005  format(1x,20a4/21x,'Beam Resultants for x-z(Statik)'/
     1'El. Q K-NR. N_x      Q_z       M_y       Press(c)',
     +'  S_(h/2)   S_(-h/2)  T_xz')
2006  format(1x,20a4/21x,'Beam Resultants for x-z(Statik)'/
     1 ' El. Q K-NR. ',' N_x     Q_z        M_y        Press(c)')
2007  format(1x,20a4/21x,'Beam Resultants for x-z(Statik)'/
     1 ' El. Q K-NR. ',' S_x     R_z        M_y        Press(c)',
     2 '   eps(<2.5!)')
2008  format(1x,'Stabkennzahl fuer Element ',i5,' zu gross !')
3000  format(' Input: E,A,I,h,q1,q2,n1,n2,2.Ordn.'/3x,'>',$)
3001  format(' Input: rho,c,b,a_T,to,tu,Numbers '/3x,'>',$)
      end
      subroutine tranb7(cs,sn,p,s,nst,ndf,isw)
      implicit double precision (a-h,o-z)
      dimension p(ndf,2),s(nst,nst),ss(2)
c.... transform the displacements to local form
      if(isw.eq.1) then
        do 90 k = 1,2
          s(1,k) = cs*p(1,k) + sn*p(2,k)
          s(2,k) =-sn*p(1,k) + cs*p(2,k)
          s(3,k) = p(3,k)
 90     continue
      else
c.... skip for an identity transformation
c        if(cs.ge.0.999999d0) return
c.... postmultiply local stiffness by transformation array
        j1 = 0
        nsiz = ndf + ndf
        do 130 k = 1,2
          do 120 i = 1,nsiz
            ss(1) = s(i,j1+1)
            ss(2) = s(i,j1+2)
            s(i,j1+1) = ss(1)*cs - ss(2)*sn
            s(i,j1+2) = ss(1)*sn + ss(2)*cs
120       continue
          j1 = j1 + ndf
130     continue
c.... premultiply result by the transpose of the transformation array
        j1 = 0
        do 190 k = 1,2
          do 160 i = 1,nsiz
            ss(1) = s(j1+1,i)
            ss(2) = s(j1+2,i)
            s(j1+1,i) = cs*ss(1) - sn*ss(2)
            s(j1+2,i) = sn*ss(1) + cs*ss(2)
160       continue
c.... transform the load vector
          ss(1)  = cs*p(1,k) - sn*p(2,k)
          p(2,k) = sn*p(1,k) + cs*p(2,k)
          p(1,k) = ss(1)
          j1 = j1 + ndf
190     continue
      end if
      return
      end
c
      subroutine qload03(qx1,qx2,qy1,qy2,d,q,propq,prop,numel,n,mqloa,
     + isw)
c----------------------------------------------------------
c.... add loads from macro qloa
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*)
      if(isw.eq.4.or.isw.eq.13) then 
        qx1 = d(7)*prop 
        qx2 = d(8)*prop 
        qy1 = d(5)*prop 
        qy2 = d(6)*prop 
        if(mqloa.ne.1) then  
          qx1 = qx1+q(n,3)*propq 
          qx2 = qx2+q(n,4)*propq 
          qy1 = qy1-q(n,1)*propq 
          qy2 = qy2-q(n,2)*propq 
        end if
      else if(isw.eq.22) then 
        qx1 = 0.d0      
        qx2 = 0.d0      
        qy1 = 0.d0      
        qy2 = 0.d0        
        if(mqloa.ne.1) then  
          qx1 = q(n,3)*propq 
          qx2 = q(n,4)*propq 
          qy1 = -q(n,1)*propq 
          qy2 = -q(n,2)*propq 
        end if
      else
        qx1 = d(7)*prop 
        qx2 = d(8)*prop 
        qy1 = d(5)*prop 
        qy2 = d(6)*prop 
      end if
      
      return
      end
