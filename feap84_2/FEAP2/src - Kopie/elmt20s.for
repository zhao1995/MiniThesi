      subroutine elmt20(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-------------------------------------------------------+
c     lineares Traegerrostelement                       |
c-------------------------------------------------------+
c   * dofs: w,phix,phiy, analog Platte                  |
c   * Gleichlast                       )                |
c   * Plot Schnittgroessen  FORC,...                    |
c     1 = Q_z, 2 = M_T, 3 = M_y                         |
c-------------------------------------------------------+
c 1.Karte   d(1)  = E                                   |
c           d(2)  = G                                   |
c           d(3)  = I_T                                 |
c           d(4)  = I_y                                 |
c           d(5)  = q_1                                 |
c           d(6)  = q_2                                 |
c           d(7)  = 0/1 Schnittgr. ohne/mit Zahlen      |
c-------------------------------------------------------+
      USE bdata
      USE cdata
      USE eldata
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
     1         d(*),ul(ndf,*),s(nst,*),p(nst),vl(3,2),sig(6)
      dimension h1(*),h2(*),h3(*)
      data ixl/1,2,3,4,1/,ixl1/1,2,3,1/
      ielno = 20
c.... Sprung zu gewuenschtem Programmteil
      go to (1,2,3,3,2,3,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,3), isw
      return
1     if(ior.lt.0) write(*,3000)
      call dinput(d,7)
                   write(iow,2000) (d(i),i=1,7)
      if(ior.lt.0) write(*  ,2000) (d(i),i=1,7)
      
c.... define node numbering for plot mesh routine, see pltord
        inord(ielno)   = 4
        ipord(1,ielno) = 1
        ipord(2,ielno) = 2
        ipord(3,ielno) = 2
        ipord(4,ielno) = 1
        ipla = 1
c.... description of stresses  
      forsus( 1) =  '  FORCE   Q_z  '
      forsus( 2) =  '  MOMENT  M_T  '
      forsus( 3) =  '  MOMENT  M_y  '
      do i = 4,11         
         forsus(i) =  ' '
      enddo
2     return
c....
c.... Transformationsmatrix
3     cs = xl(1,2) - xl(1,1)
      sn = xl(2,2) - xl(2,1)
      dl = dsqrt(cs*cs + sn*sn)
      cs = cs/dl
      sn = sn/dl
c.... Transformation der Verschiebungen auf lokale Richtungen
      call trans20(cs,sn,ul,vl,3,ndf,1)
c.... Steifigkeitsmatrix
      gitl   = d(2)*d(3)/dl
      eiyl   = d(1)*d(4)/(dl*dl*dl)
c.... Torsion
      s(2,2) =  gitl
      s(2,5) = -gitl
      s(5,5) =  gitl
c...  Biegung
      s(1,1) =  12.d0* eiyl
      s(1,3) =  -6.d0* eiyl * dl
      s(1,4) = -12.d0* eiyl
      s(1,6) =  -6.d0* eiyl * dl
      s(3,3) =   4.d0* eiyl * dl * dl
      s(3,4) =   6.d0* eiyl * dl
      s(3,6) =   2.d0* eiyl * dl * dl
      s(4,4) =  12.d0* eiyl
      s(4,6) =   6.d0* eiyl * dl
      s(6,6) =   4.d0* eiyl * dl * dl
c.... Symmetrie
      s(3,1) = s(1,3)
      s(4,1) = s(1,4)
      s(6,1) = s(1,6)
      s(4,3) = s(3,4)
      s(6,3) = s(3,6)
      s(6,4) = s(4,6)
      s(5,2) = s(2,5)
c
c....   Trapez--Lastvektor
        call qload20(q1,q2,d,aqloa,numel,n,mqloa,propq,prop,isw)
        p(1) =-0.05d0 * dl *      (     7.d0*q1 +      3.d0*q2)
        p(4) =-0.05d0 * dl *      (     3.d0*q1 +      7.d0*q2)
        p(3) = 0.05d0 * dl * dl * (          q1 + 2.d0/3.d0*q2)
        p(6) =-0.05d0 * dl * dl * (2.d0/3.d0*q1 +           q2)
c
c     Residuum R = K * v - P
      do 32 i = 1,6
         do 32 k = 1,3
            p(i) = p(i) - s(i,k)*vl(k,1) - s(i,k+3)*vl(k,2)
32    continue
      if(isw.eq.4.or.isw.eq.13) goto 4
c.... Transformation von K und P auf globale Richtungen
      call trans20(cs,sn,p,s,nst,ndf,2)
      return
c.... Berechnung der Schnittgroessen   -R = K*v-P 
4     inum = d(8)
      do i = 1,3
        sig(i) =     p(i)
        sig(i+3) = - p(i+3)
      enddo
      if(isw.eq.4) then
c....   drucke Schnittgroessen
        mct = mct - 3
        if(mct.lt.0) then
          if(ior.lt.0) write(*  ,2001) head
                       write(iow,2001) head
          mct = 50
        endif
        if(ior.lt.0) write(*  ,2002)
     1     n,ma,' L ',(sig(i),i=1,3),' R ',(sig(i),i=4,6)
                       write(iow,2002) 
     1     n,ma,' L ',(sig(i),i=1,3),' R ',(sig(i),i=4,6)

      elseif(isw.eq.13) then
c....   plotte Schnittgroessen
        if(nfp.lt.1.or.nfp.gt.3) return
        nd2 = ndf+nfp
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
c.....    Lasten
          call qload20(q1,q2,d,aqloa,numel,n,mqloa,propq,prop,isw)
c.....    Auswertung fuer  Q,MT,Mz,inkrementweise Berechnung ns
          ns = 10
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
c.....      dq             
            dq1 = (q2-q1)*r1/dr0 
            dq2 = (q2-q1)*r2/dr0
c.....      Q
            sq1 = sig(1) + q1*r1 + dq1*r1/2.  
            sq2 = sig(1) + q1*r2 + dq2*r2/2. 
c.....      M_T   
            smt1= sig(2)+(sig(5)-sig(2))*r1/dr0
            smt2= sig(2)+(sig(5)-sig(2))*r2/dr0
c.....      M            
            sm1 = sig(3) + sig(1)*r1 + q1*r1*r1/2. + dq1*r1*r1/6. 
            sm2 = sig(3) + sig(1)*r2 + q1*r2*r2/2. + dq2*r2*r2/6.
c
            if(nfp.eq.1) then
              s1 = sq1
              s2 = sq2
            elseif(nfp.eq.2) then
              s1 = smt1
              s2 = smt2
            elseif(nfp.eq.3) then
              s1 = sm1
              s2 = sm2
            endif
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
            endif
            if(is.eq.ns.and.inum.eq.1) then
              x5 = x2 + sn*s2*fac
              y5 = y2 - cs*s2*fac
              call plotl(x5,y5,0.0d0,3)
              call pppcol(1)
              ns2 = s2/cfp
              call plabl(ns2)
            endif
c.....      Berechne Farbe                                
            sm = 0.5d0*(s1+s2)/cfp
            call pppcolf(sm)
c.....      plotte Schnittgroesse 
            if(s1.ge.0.0.and.s2.ge.0.0) then
              goto 136
            elseif(s1.le.0.0.and.s2.le.0.0) then
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
              endif
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
            endif
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
        endif
      endif
      return
c.... format statements
2000  format(5x,'Traegerrost, Bernoulli+St.Venant'//
     +   10x,'E  ',e14.5/10x,'G  ',e14.5/10x,'I_T',e14.5/
     +   10x,'I_y',e14.5/10x,'q1 ',e14.5/10x,'q2 ',e14.5/
     +   10x,'Zahlen f. N,M,Q  0=nein/1=ja',e14.5)
2001  format(1x,20a4/21x,'Schnittgroessen im lokalen KOS'/
     1'El. Q K-NR. Q_z         M_T          M_y')
2002  format(2i3,a3,3(1x,g12.5),/,6x,a3,3(1x,g12.5))
3000  format(' Input: E,G,I_T,I_y,q1,q2'/3x,'>',$)
      end
c
      subroutine trans20(cs,sn,p,s,nst,ndf,isw)
c-----------------------------------------------
c.... transformations                           
c-----------------------------------------------
      implicit double precision (a-h,o-z)
      dimension p(ndf,2),s(nst,nst),ss(3)
c.... transform the displacements to local form
      if(isw.eq.1) then
        do 90 k = 1,2
          s(1,k) =   -p(1,k)
          s(2,k) = cs*p(2,k) + sn*p(3,k)
          s(3,k) = sn*p(2,k) - cs*p(3,k)
 90     continue
      else
c.... skip for an identity transformation
c        if(cs.ge.0.999999d0) return
c.... postmultiply local stiffness by transformation array
        j1 = 0
        nsiz = ndf + ndf
        do 130 k = 1,2
          do 120 i = 1,nsiz
            ss(2) = s(i,j1+2)
            ss(3) = s(i,j1+3)
            s(i,j1+1) = -s(i,j1+1)
            s(i,j1+2) = ss(2)*cs + ss(3)*sn
            s(i,j1+3) = ss(2)*sn - ss(3)*cs
120       continue
          j1 = j1 + ndf
130     continue
c.... premultiply result by the transpose of the transformation array
        j1 = 0
        do 190 k = 1,2
          do 160 i = 1,nsiz
            ss(2) = s(j1+2,i)
            ss(3) = s(j1+3,i)
            s(j1+1,i) = -s(j1+1,i)
            s(j1+2,i) = cs*ss(2) + sn*ss(3)
            s(j1+3,i) = sn*ss(2) - cs*ss(3)
160       continue
c.... transform the load vector
          ss(2)  = cs*p(2,k) + sn*p(3,k)
          p(1,k) =   -p(1,k) 
          p(3,k) = sn*p(2,k) - cs*p(3,k)
          p(2,k) = ss(2)
          j1 = j1 + ndf
190     continue
      endif
      return
      end
c
      subroutine qload20(q1,q2,d,q,numel,n,mqloa,propq,prop,isw)
c----------------------------------------------------------
c.... set loads from macro qloa/mate
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*)
      if(isw.eq.22) then
        q1 = 0.d0
        q2 = 0.d0
        if(mqloa.ne.1) then 
          q1 = q(n,1)*propq 
          q2 = q(n,2)*propq 
        end if
      else if(isw.eq.13) then
        q1 = d(5)*prop
        q2 = d(6)*prop
        if(mqloa.ne.1) then 
          q1 = q1 + q(n,1)*propq 
          q2 = q2 + q(n,2)*propq 
        end if
      else
        q1 = d(5)*prop
        q2 = d(6)*prop
      end if
      return
      end
