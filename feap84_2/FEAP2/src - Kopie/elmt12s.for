      subroutine elmt12(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c----------------------------------------------------------------------+
c
c.... two dim. laplace equation element,  3-9 node isoparametric       
c.....    ityp  =  1  -  heat transfer                                 
c.....          =  2  -  flow in porous media - horizontal slab        
c.....          =  3  -  flow in porous media - vertical slab          
c.....          =  4  -  flow in porous media - simulated free surface 
c.....          =  5  -  warping function for arbitrary cross sections 
c                        including: all cross section properties       
c                        (with or without loads possible)              
c                        symmetry with respect to x, y, x/y axis       
c.....          =  6  -  torsion function for arbitrary cross sections 
c                        boundary conditions by edge                   
c                        symmetry with respect to x, y, x/y axis       
c.....          =  7  -  warping function for arbitrary cross sections 
c                        from shear forces                             
c                        symmetry with respect to x, y, x/y axis       
c.....          =  8  -  stress  function for arbitrary cross sections 
c                        from shear forces  with  nue                  
c                        symmetry with respect to x, y, x/y axis       
c
c     plot flow/shear stress as arrows with macro flux                 
c
c     D-Array 
c     ITYP  1                 2                3               4
c     d(1)  1-conductivity    1-permeability   1-permeability  1-permeability 
c     d(2)  2-conductivity    2-permeability   2-permeability  2-permeability
c     d(3)  1-x angle deg.    1-x angle deg.   1-x angle deg.  1-x angle deg.
c     d(4)  heat source*rho   flow source*rho  flow source*r   flow source*r 
c     d(5)  spec. heat *rho   spec. flow *rho  spec. flow *r   spec. flow *r
c     d(6)  not used          not used         rho*g           not used
c     d(7)  not used          not used         surface elev.   not used 
c     d(8-17) not used        not used         not used        not used 
c           
c           
c     ITYP  5               6    
c     d(1)  1-shear modulus 1.d0  
c     d(2)  2-shear modulus 1.d0  
c     d(3)  0.d0            0.d0   
c     d(4)  0.d0            2.d0   
c     d(5)  0.d0            0.d0   
c     d(5)  0.d0            0.d0   
c     d(7-17) not set       not set 
            
c     ITYP  7/8                       
c     d(1)  1.d0   
c     d(2)  1.d0   
c     d(3)  0.d0                  
c     d(4)  0.d0                  
c     d(5)  0.d0                  
c     d(6)  0.d0                  
c     d(7)  area .............. A     
c     d(8)  1-moment of inertia I_y   
c     d(9)  2-moment of inertia I_z   
c     d(10) 2-moment of inertia I_yz  
c     d(11) center of gravity   y_s   
c     d(12) center of gravity   z_s   
c     d(13) 1-shear force       Q_y   
c     d(14) 2-shear force       Q_z   
c     d(15) Poisson ratio       nue   
c     d(16) center of stress f. y_0   
c     d(17) center of stress f. z_0   
c     elmt no. used for plot, see inord, ipord                      
c                                                                   
c     Baustelle
c     stand 1.8.00 Einbau: unterschiedliche Gs, offen: g12 und 
c                          drehwinkel im Typ 5, noch keine doku
c                          el geht mit g1=g2=1 fuer Eingabe 0
c     stand davor in elmt12s.old           
c----------------------------------------------------------------------+
      USE bdata
      USE cdata
      USE eldata
      USE errin2
      USE errin3
      USE iofile
      USE pdata6
      USE pdata7
      USE strnam
      implicit double precision (a-h,o-z)
      logical fflg4, flgsec(4)
      dimension xl(ndm,*),tl(*),d(*),ul(*),s(nst,*),p(*),ix(*)
     1         ,shp(3,9),sg(16),tg(16),wg(16),td(3),ipord3(4),xa(19)
      dimension h1(*),h2(*),h3(*)
      character wlab(2)*12
      save ityp,kat,fflg4,flgsec,xa,ksym
cww      data  sg/-1.0, 1.0,1.0,-1.0,0.0/
cww      data  tg/-1.0,-1.0,1.0, 1.0,0.0/
      data wlab/'  p l a n e ','axisymmetric'/
      data ipord3 /1,2,3,1/
      if(ndm.eq.3) ipla = 1
      ielno = 12
c.... transfer to correct processor
      go to (1,2,3,4,5,3,2,8,9,2,2,2,2),isw
      return
c.... input material properties
1     if(ior.lt.0) write(*,1000)
      call dinput(td,3)
      ityp = td(1)
      kat  = td(2)
      ityp = max(1,min(8,ityp))
      if (ityp.ge.5) fflg4 = .true.
      if (ityp.ge.5) ksym = td(3)
      if (kat.ne.2) kat=1
      if (ior.lt.0.and.ityp.eq.1) write(*,1001)
      if (ior.lt.0.and.ityp.eq.2) write(*,1002)
      if (ior.lt.0.and.ityp.eq.3) write(*,1003)
      if (ior.lt.0.and.ityp.eq.4) write(*,1004)
      if (ior.lt.0.and.ityp.eq.5) write(*,1005)
      if (ityp.le.5) call dinput(d(1),7)
      if (ityp.eq.5.and.dabs(d(1)).lt.1.0d-10) d(1) = 1.d0 ! no input for ityp=5 poss.
      if (ityp.gt.5) then
        d(1) = 1.0d0
        d(2) = 1.0d0
        d(3) = 0.0d0
        d(4) = 0.0d0
        d(5) = 0.0d0
        d(6) = 0.0d0
  
        if (ityp.eq.6) d(4) = 2.0d0
        if (ityp.eq.7) call dinput(d(7),11)!  7-A,8-Iy,9-Iz,10-Iyz
c                                          ! 11-ys,12-zs,13-Qy,14-Qz
c                                          ! 15-nue,16-y0,17-z0
        if (ityp.eq.8) call dinput(d(7),11)!  7-A,8-Iy,9-Iz,10-Iyz
c                                          ! 11-ys,12-zs,13-Qy,14-Qz
c                                          ! 15-nue,16-y0,17-z0
      end if 
      if (ityp.ge.5) then
        do i=1,4
         flgsec(i) = .true.
        end do
      end if
      if (dabs(d(2)).lt.1.0d-10) then
        d(2) = d(1)
        d(3) = 0.0d0
      end if
      if     (ityp.eq.1) then
                        write(iow,2001) (d(i),i=1,5),wlab(kat)
        if   (ior.lt.0) write(*  ,2001) (d(i),i=1,5),wlab(kat)
      else if (ityp.eq.2) then
                        write(iow,2002) (d(i),i=1,5),wlab(kat)
        if   (ior.lt.0) write(*  ,2002) (d(i),i=1,5),wlab(kat)
      else if (ityp.eq.3) then
                        write(iow,2003) (d(i),i=1,7),wlab(kat)
        if   (ior.lt.0) write(*  ,2003) (d(i),i=1,7),wlab(kat)
      else if (ityp.eq.4) then
                        write(iow,2004) (d(i),i=1,5),wlab(kat)
        if   (ior.lt.0) write(*  ,2004) (d(i),i=1,5),wlab(kat)
      else if (ityp.eq.5) then
                        write(iow,2005) (d(i),i=1,2),wlab(kat)
        if   (ior.lt.0) write(*  ,2005) (d(i),i=1,2),wlab(kat)
      else if (ityp.eq.6) then
                        write(iow,2006) wlab(kat)
        if   (ior.lt.0) write(*  ,2006) wlab(kat)
      else if (ityp.eq.7) then
                        write(iow,2007) (d(i),i=7,17),wlab(kat)
        if   (ior.lt.0) write(*  ,2007) (d(i),i=7,17),wlab(kat)
      else if (ityp.eq.8) then
                        write(iow,2008) (d(i),i=7,17),wlab(kat)
        if   (ior.lt.0) write(*  ,2008) (d(i),i=7,17),wlab(kat)
      end if
c.... description of 'stresses'  
      do i = 1,24
         strsus(i) = ' '                    
      end do
      if(ityp.le.4) then
        strsus( 1) = '  1-FLOW       '
        strsus( 2) = '  2-FLOW       '
        strsus( 3) = '  MAX FLOW     '
        strsus( 4) = '  PRESSURE     '
        if(ityp.eq.1) strsus( 4) = '  TEMPERATURE  '
        if(ityp.eq.4) strsus( 4) = '  HEAD         '
      else if(ityp.ge.5) then
        strsus( 1) = '  STRESS T_13  '
        strsus( 2) = '  STRESS T_23  '
        strsus( 3) = '  STRESS T_max '
        strsus( 4) = '  WARPING WQ   '
        strsus( 5) = '  WARPING WT   '
        if(ityp.eq.6) strsus( 4) = 'STRESS FUNCTION'
        if(ityp.eq.7) strsus( 4) = '  WARPING WQ   '
        if(ityp.ge.6) strsus( 5) = '               '
      end if               
c
c.... define node numbering for plot mesh routine, see pltord
      if(nen.eq.3) then
        inord(ielno) = 4
        do ii = 1,4
          ipord(ii,ielno) = ipord3(ii)
        end do
      end if
c
2     return
c.... compute conductivity (stiffness) matrix
3     ll=2
      if(nel.gt.4) ll = 3
      call pgauss(ll,lint,sg,tg,wg)
      do 303 l = 1,lint
        call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
        xsj = xsj*wg(l)
        if (kat.eq.2) then
          rr = 0.0d0
          do 301 i = 1,nel
            rr = rr+shp(3,i)*xl(1,i)
301       continue
          xsj = xsj*rr
        end if
cww        if (ityp.eq.7 .or. ityp.eq.8) then
        if (ityp.eq.5 .or. ityp.eq.7 .or. ityp.eq.8) then
          yy = 0.d0
          zz = 0.d0
          do  i = 1,nel
            yy = yy+shp(3,i)*xl(1,i)
            zz = zz+shp(3,i)*xl(2,i)
          end do
c          if(ityp.eq.7) then
            y7 = yy - d(11)      ! coordinates   yq = y - ys
            z7 = zz - d(12)
c          else if(ityp.eq.8) then
            y8 = yy - d(16)      ! coordinates   yq = y - y0
            z8 = zz - d(17)
c          end if
        end if
        j1 = 1
        do 303 j = 1,nel
c...      heat/flow internally generated (source) 
          if (ityp.lt.5) then
            p(j1) = p(j1)+d(4)*shp(3,j)*xsj
          else if (ityp.eq.5) then
            p(j1) = p(j1)+(d(1)*zz*shp(1,j)-d(2)*yy*shp(2,j))*xsj
          else if (ityp.eq.7) then
c...        load terms warping Q  
c           p=int[Q_z/D*(z*I_z-y*Iyz)+Q_y/D*(y*I_y-z*I_yz)]*N_K da
            det = d(8)*d(9)-d(10)*d(10)
            p(j1) = p(j1) + ( d(14)/det*(d(9)*z7-d(10)*y7)
     +                       +d(13)/det*(d(8)*y7-d(10)*z7))
     +                    * shp(3,j)*xsj
c....       statt eloa   
            a1  =  (d(13)*d(8)-d(14)*d(10))/det      
            a2  =  (d(14)*d(9)-d(13)*d(10))/det      
            g1  =  -d(15)/(1.d0+d(15))*0.5d0*a1*z8*z8
            g2  =   d(15)/(1.d0+d(15))*0.5d0*a2*y8*y8

            p(j1) = p(j1)+(g1*shp(1,j)-g2*shp(2,j))*xsj
c
          else if (ityp.eq.8) then
c..         load term stress function Q
c           p=int[nu/(1+nu)1/D
c           [(Q_z*I_z-Q_y*Iyz)*(y-y0)-(Q_y*I_y-Q_z*I_yz)*(z-z_0)]*N_K da
            det = d(8)*d(9)-d(10)*d(10)
            p(j1) = p(j1) + d(15)/(1.d0+d(15))/det*
     +      ((d(14)*d(9)-d(13)*d(10))*y8-(d(13)*d(8)-d(14)*d(10))*z8)
     +                   * shp(3,j)*xsj
          end if
c...      stiffness matrix
          a1 = d(1)*shp(1,j)*xsj
          a2 = d(2)*shp(2,j)*xsj
          i1 = 1
          do 302 i = 1,nel
            s(i1,j1) = s(i1,j1) + (a1*shp(1,i) + a2*shp(2,i))
            i1 = i1 + ndf
302       continue
          j1 = j1 + ndf
303   continue
c.... compute the residual
      do 304 i = 1,nst
        do 304 j = 1,nst
          p(i) = p(i) - s(i,j)*ul(j)
304   continue
      return
c.... print the flows(stresses) in each element
4     ll = 1
      if(nel.gt.4) ll = 2
      call pgauss(ll,lint,sg,tg,wg)
      if(ityp.le.4) then
      do 410 l = 1,lint
        call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
        xsj = xsj*wg(l)
        call flolag(d,shp,xl,ul,rr,zz,q1,q2,qm,pr,ndm,ndf,ityp)
        mct = mct - 1
        if(mct.lt.0) then
          mct = 50
          if(ityp.eq.1) then
                         write(iow,4001) o,head
            if(ior.lt.0) write(*  ,4001) o,head
          else
                         write(iow,4002) o,head
            if(ior.lt.0) write(*  ,4002) o,head
          end if
        end if
                     write(iow,4003) n,ma,rr,zz,q1,q2,qm,pr
        if(ior.lt.0) write(*  ,4003) n,ma,rr,zz,q1,q2,qm,pr
410   continue
      else if(ityp.ge.5) then
        if(.not.flgsec(1).and..not.flgsec(2).and..not.flgsec(3).and.
     +    .not.flgsec(4)) then
c....     print shear stresses
          mct = mct - 1
          if(mct.lt.0) then
            mct = 50
                         write(iow,4005) o,head
            if(ior.lt.0) write(*  ,4005) o,head
          end if
          do 411 l = 1,lint
            call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
            xsj = xsj*wg(l)
            call sectau(d,xl,ul,shp,ndm,ndf,nel,ityp,x1,x2,t13,t23,tmax,
     +                  xa,ksym)
                         write(iow,4003) n,ma,x1,x2,t13,t23
            if(ior.lt.0) write(*  ,4003) n,ma,x1,x2,t13,t23
411       continue
        else
c....     compute resultant warping and section properties
          call secpro(d,xl,ul,ix,shp,ndm,ndf,xa,flgsec,ityp,ksym)
        end if
      end if
      return 
c.... compute heat capacity (mass) matrix
5     if(ityp.gt.4) return
      ll=2
      if(nel.gt.4) ll = 3
      call pgauss(ll,lint,sg,tg,wg)
      do 502 l = 1,lint
      call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
      xsj = xsj*wg(l)
      if (kat.eq.2) then
        rr = 0.0d0
        do 501 i = 1,nel
          rr = rr+shp(3,i)*xl(1,i)
501     continue
        xsj = xsj*rr
      end if
      j1 = 1
      do 502 j = 1,nel
        p(j1) = p(j1) + d(5)*shp(3,j)*xsj
        j1 = j1 + ndf
502   continue
      return
c.... compute the nodal flow values
8     ll=2
      call pgauss(ll,lint,sg,tg,wg)
      if(ityp.le.4) then
        istv = -4
        call stcnlag(ix,d,xl,ul,shp,strea,strea(1+numnp),ndf,2,nel,
     1               numnp,sg,tg,wg,lint,ityp)
      else if(ityp.ge.5) then
        istv = -5
        if(.not.flgsec(1).and..not.flgsec(2).and..not.flgsec(3).and.
     +    .not.flgsec(4)) 
     +     call stcntors(ix,d,xl,ul,shp,strea,strea(1+numnp),
     +           ndf,ndm,nel,numnp,sg,tg,wg,lint,xa,ksym,ityp)

      end if
      return
c                
c.... calculate errors
9     ll=2
      if(nel.gt.4) ll = 3
      call pgauss(ll,lint,sg,tg,wg)
      call stcn12e(ix,d,xl,ul,shp,strea,strea(1+numnp),ndf,ndm,nel,
     1             numnp,numel,sg,tg,wg,lint,ityp,e_ome,xa,ksym)
      return
c
c.... formats
1000  format(' Input:itype(1=heat,2=flow-horiz,3=flow-vert,4=flow-free',
     1       ',5=warping)'/'       geom.(1=plane,2=axisym)'/3x,'>',$)
1001  format(' Input:K-1, K-2, angle 1-x, rho*h, rho*c '/3x,'>',$)
1002  format(' Input:K-1, K-2, angle 1-x, rho*h, rho*c '/3x,'>',$)
1003  format(' Input:K-1, K-2, angle 1-x, rho*h, rho*c, rho*g, 
     1         surf.elev.'/3x,'>',$)
1004  format(' Input:K-1, K-2, angle 1-x, rho*h, rho*c '/3x,'>',$)
1005  format(' Input:G-1, G-2'/3x,'>',$)
2001  format(3x,'Linear heat conduction element'//
     +  5x,'1-conductivity       ',e12.4/
     +  5x,'2-conductivity       ',e12.4/
     +  5x,'1-x angle deg.       ',e12.4/
     +  5x,'heat source*density  ',e12.4/
     +  5x,'specific heat*density',e12.4/
     +  5x,a12,' analysis')
2002  format(3x,'Horizontal flow in porous media'//
     +  5x,'1-permeability       ',e12.4/
     +  5x,'2-permeability       ',e12.4/
     +  5x,'1-x angle deg.       ',e12.4/
     +  5x,'flow source*density  ',e12.4/
     +  5x,'specific flow*density',e12.4/
     +  5x,a12,' analysis')
2003  format(3x,'Vertical flow in porous media'//
     +  5x,'1-permeability       ',e12.4/
     +  5x,'2-permeability       ',e12.4/
     +  5x,'1-x angle deg.       ',e12.4/
     +  5x,'flow source*density  ',e12.4/
     +  5x,'specific flow*density',e12.4/
     +  5x,'density *g           ',e12.4/
     +  5x,'surface elev.        ',e12.4/
     +  5x ,a12,' analysis')
2004  format(3x,'Simulated free surface porous media'//
     +  5x,'1-permeability       ',e12.4/
     +  5x,'2-permeability       ',e12.4/
     +  5x,'1-x angle deg.       ',e12.4/
     +  5x,'flow source*density  ',e12.4/
     +  5x,'specific flow*density',e12.4/
     +  5x,a12,' analysis')
2005  format(3x,'Warping function for arbitrary cross sections'//
     1  5x,'1-shear modulus',e12.4/
     2  5x,'2-shear modulus',e12.4/
     3  5x,a12,' analysis') 
2006  format(3x,'Stress function for arbitrary cross sections'/
     1  5x,a12,' analysis')
2007  format(3x,'Warping function from Q for arbitrary cross sections'/
     1  5x,'area .............. A    ',e12.4/
     2  5x,'1-moment of inertia I_y  ',e12.4/
     3  5x,'2-moment of inertia I_z  ',e12.4/
     4  5x,'2-moment of inertia I_yz ',e12.4/
     5  5x,'center of gravity   y_s  ',e12.4/
     6  5x,'center of gravity   z_s  ',e12.4/
     7  5x,'1-shear force       Q_y  ',e12.4/
     8  5x,'2-shear force       Q_z  ',e12.4/
     9  5x,'Poisson ratio       nue  ',e12.4/
     +  5x,'center of stress f. y_0  ',e12.4/
     1  5x,'center of stress f. z_0  ',e12.4/
     2  5x,a12,' analysis')
2008  format(3x,'Stress function from Q for arbitrary cross sections'/
     1  5x,'area .............. A    ',e12.4/
     2  5x,'1-moment of inertia I_y  ',e12.4/
     3  5x,'2-moment of inertia I_z  ',e12.4/
     4  5x,'2-moment of inertia I_yz ',e12.4/
     5  5x,'center of gravity   y_s  ',e12.4/
     6  5x,'center of gravity   z_s  ',e12.4/
     7  5x,'1-shear force       Q_y  ',e12.4/
     8  5x,'2-shear force       Q_z  ',e12.4/
     9  5x,'Poisson ratio       nue  ',e12.4/
     +  5x,'center of stress f. y_0  ',e12.4/
     1  5x,'center of stress f. z_0  ',e12.4/
     2  5x,a12,' analysis')
4001  format(a1,20a4//'  H e a t    t r a n s f e r'//
     1 ' elem  mat',4x,'1-coord',4x,'2-coord',4x,
     2 '1-flow',5x,'2-flow',3x,'max flow',1x,'temperature'/)
4002  format(a1,20a4//'  F l o w    i n    p o r o u s    m e d i a'//
     1 ' elem  mat',4x,'1-coord',4x,'2-coord',4x,
     2 '1-flow',5x,'2-flow',3x,'max flow',3x,'pressure'/)
4003  format(2i5,0p2f11.3,1p4e11.3)
4005  format(a1,20a4//
     +        '  W a r p i n g   f o r   c r o s s   s e c t i o n s '//
     1 ' elem  mat',4x,'1-coord',4x,'2-coord',2x,
     2 '13-stress',2x,'23-stress'/)
      end
c
      subroutine flolag(d,shp,xl,ul,rr,zz,q1,q2,qm,pr,ndm,ndf,ityp)
c----------------------------------------------------------------------+
c     calculate flows                                                  |
c----------------------------------------------------------------------+
      USE eldata
      implicit real*8 (a-h,o-z)
      real*8 xl(ndm,*)
      real*8 ul(ndf,*),d(*),shp(3,*)
c.... compute flows at current point
      rr = 0.d0
      zz = 0.d0
      q1 = 0.d0
      q2 = 0.d0
      pr = 0.d0
      do 401 i = 1,nel
        ull = ul(1,i)
        rr = rr + shp(3,i)*xl(1,i)
        zz = zz + shp(3,i)*xl(2,i)
        q1 = q1 + shp(1,i)*ull
        q2 = q2 + shp(2,i)*ull
        if(ityp.eq.4) then
          if (ull.gt.0.0d0) then
            ull = sqrt(ull+ull)
          else if(ull.lt.0.0d0) then
            write(*,2000) i,n
            ull = 0.d0
          end if
        end if
        pr = pr + shp(3,i)*ull
401   continue
      q1 = -d(1)*q1
      q2 = -d(2)*q2
      qm =  sqrt(q1*q1 + q2*q2)
      if(ityp.eq.3) pr = d(6)*(pr + d(7)-zz)
      return
2000  format(' WARNING ** head negative at local node',i2,
     1  ' of element',i4,'.  Set zero.')
      end
c
      subroutine stcnlag(ix,d,xl,ul,shp,dt,st,ndf,ndm,nel,numnp,
     1                  sg,tg,wg,lint,ityp)
c----------------------------------------------------------------------+
c     plot flows and pressure/temperature                              |
c----------------------------------------------------------------------+
      implicit real*8 (a-h,o-z)
      integer*4 ix(*)
      real*8 dt(numnp),st(numnp,*),xl(ndm,*),sg(*),tg(*),wg(*)
      real*8 ul(ndf,*),shp(3,*),d(*)
      do 301 l = 1,lint
        call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
        xsj = xsj*wg(l)
        call flolag(d,shp,xl,ul,rr,zz,q1,q2,qm,pr,ndm,ndf,ityp)
c....   sum flows for error analysis
        call enerel12(d,q1,q2,xsj)
c
        do 300 ii = 1,nel
          ll = abs(ix(ii))
          if(ll.gt.0) then
            pr   = ul(1,ii)
            if(ityp.eq.3) pr = d(6)*(pr + d(7) - xl(2,ii))
            if(ityp.eq.4) pr = sqrt(abs(pr+pr))
            xsji = xsj*shp(3,ii)
            dt(ll)   = dt(ll)   +    xsji
            st(ll,1) = st(ll,1) + q1*xsji  ! flow_x
            st(ll,2) = st(ll,2) + q2*xsji  ! flow_y
            st(ll,3) = st(ll,3) + qm*xsji  ! max flow
            st(ll,4) = st(ll,4) + pr*xsji  ! pressure/temperature
          end if
300     continue
301   continue
      return
      end
c
      subroutine secpro(d,xl,ul,ix,shp,ndm,ndf,xa,flgsec,ityp,ksym)
c-----------------------------------------------------------------------
c     section properties and warping function
c     output on stre,all:
c     1st call:   xa(1) = int[  ] dA      => A
c                 xa(2) = int[ y ]dA      => S_z
c                 xa(3) = int[ z ]dA      => S_y
c                 xa(4) = int[y*y]dA      => I_z
c                 xa(5) = int[z*z]dA      => I_y
c                 xa(6) = int[y*z]dA      => I_yz
c                 xa(7) = int[ w ]dA      => I_w ,for w^ = w - 1/A{int[w]dA}
c                 xa(14)= y_S
c                 xa(15)= z_S
c        (ityp5)  xa(16)= int[t23*y ]dA   (ityp6)  xa(16)= int[w*y]dA 
c        (ityp5)  xa(17)= int[t13*z ]dA   (ityp6)  xa(17)= int[w*z]dA 
c        (ityp5)  xa(18)= int[t23*y^2]dA                                 
c        (ityp5)  xa(19)= int[t13*z^2]dA                                 
c
c     2nd call:   xa(8) = int[w^*y]dA
c                 xa(9) = int[w^*z]dA
c                 xa(12) = y_M
c                 xa(13) = z_M
c     3rd call:   xa(10) = int[w~*w~]dA   => C_M
c                 xa(11) = int[ *** ]dA   => I_T    Paper FG/RS/WW, p.12
c
c     ITYP=7
c     output on stre,all:
c                 xa(1) = int[w,y]dA    w*,y = w,y -1/A*int[w,y]
c                 xa(2) = int[w,z]dA 
c                 xa(4) = int[w  ]dA 
c     2nd call:   xa(3) = int[T_xz*y-T_xy*z]dA
c                 xa(5) = int[T_xy*T_xy]dA 
c                 xa(6) = int[T_xz*T_xz]dA 
c     3rd call:   --
c
c-----------------------------------------------------------------------
      USE cdata
      USE eldata
      USE iofile
      implicit double precision(a-h,o-z)
      logical flgsec(4)
      dimension d(*),xl(ndm,*),ul(ndf,*),ix(*),sg(16),tg(16),wg(16),
     +          shp(3,*),xa(*)
      save xiys,xizs,xiyzs,ym,zm,ys,zs
      save gitww1
      data tol12 /1.e-5/
      pi = datan(1.0d0)*4.d0
      if(flgsec(1)) gitww1=0.d0
      if(flgsec(1)) call pzero(xa,19)
      if(flgsec(1)) flgsec(1) = .false.
      facsym = 1.0d0
      if (ksym.ne.0) facsym = 2.0d0
      if (ksym.eq.3) facsym = 4.0d0
      if (ityp.eq.5.or.ityp.eq.6) then 
c....   loop over gauss points
        ll=2
        if(nel.gt.4) ll = 3
        call pgauss(ll,lint,sg,tg,wg)
        do 100 l = 1,lint
          call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
          xsj = xsj*wg(l)
          w = 0.0d0
          y = 0.0d0
          z = 0.0d0
          do ii = 1,nel
            w = w + shp(3,ii)*ul(1,ii)
            y = y + shp(3,ii)*xl(1,ii)
            z = z + shp(3,ii)*xl(2,ii)
          end do       
          if (flgsec(2)) then
            xa(1) = xa(1) +       xsj
            xa(2) = xa(2) + y   * xsj
            xa(3) = xa(3) + z   * xsj
            xa(4) = xa(4) + y*y * xsj
            xa(5) = xa(5) + z*z * xsj
            xa(6) = xa(6) + y*z * xsj
            xa(7) = xa(7) + w   * xsj
            if(ityp.eq.5) then 
              call sectau(d,xl,ul,shp,ndm,ndf,nel,ityp,x1,x2,t13,t23,
     +                    tmax,xa,ksym)
              xa(16)= xa(16)+ t23*y   * xsj
              xa(17)= xa(17)+ t13*z   * xsj
              xa(18)= xa(18)+ t23*y*y * xsj
              xa(19)= xa(19)+ t13*z*z * xsj
            else if(ityp.eq.6) then 
              xa(16)= xa(16)+ w*y * xsj
              xa(17)= xa(17)+ w*z * xsj
            end if
          end if
          if(.not. flgsec(2) .and. flgsec(3)) then 
            wq = w - xa(7)/xa(1)
            if (ksym.ne.0) wq = w              ! in case of symmetry
            xa(8) = xa(8) + wq * y * xsj
            xa(9) = xa(9) + wq * z * xsj
          end if
          if(.not. flgsec(3) .and. flgsec(4)) then 
            wq = w - xa(7)/xa(1)
            if (ksym.ne.0) wq = w              ! in case of symmetry
            xws   = wq + ym*(z-zs) - zm*(y-ys)
            xa(10)= xa(10) + xws*xws *xsj
            wy = 0.0d0
            wz = 0.0d0
c
            do ii = 1,nel
              wi = xa(7)/xa(1)
              if (ksym.ne.0) wi = 0.0d0
              ww = (ul(1,ii) - wi) + ym*xl(2,ii) - zm*xl(1,ii)
cww           correct but same result for:
cww           ww = (ul(1,ii)-wi)
cww     +             + ym*(xl(2,ii)-xa(15)) - zm*(xl(1,ii)-xa(14))
              wy = wy + shp(1,ii)*ww
              wz = wz + shp(2,ii)*ww
            end do   
            yym = y-ym
            zzm = z-zm
cww g=1     xa(11) = xa(11) + (yym*yym + zzm*zzm + wz*yym - wy*zzm)*xsj
            xa(11) = xa(11) + d(2)*yym*(yym + wz)*xsj
     1                      + d(1)*zzm*(zzm - wy)*xsj
c...test ww alternative calculation of I_T
cww            wy = 0.0d0
cww            wz = 0.0d0
cww            do ii = 1,nel
cww              ww = ul(1,ii)
cww              wy = wy + shp(1,ii)*ww
cww              wz = wz + shp(2,ii)*ww
cww            end do   
cww            gitww1 = gitww1 + d(2)*y*(y+wz)*xsj + d(1)*z*(z-wy)*xsj
C
          end if
100     continue
c....   print results of first call
        if(.not.flgsec(1) .and.flgsec(2) .and. n.eq.numel) then
          flgsec(2) = .false.
          ys = xa(2) / xa(1)
          zs = xa(3) / xa(1)
          if (ksym.eq.1 .or. ksym.eq.3) zs = 0.0d0   ! symmetry on 1-axis
          if (ksym.eq.2 .or. ksym.eq.3) ys = 0.0d0   ! symmetry on 2-axis
          xa(14)= ys
          xa(15)= zs
          xiys  = (xa(5) - zs*zs*xa(1)) * facsym
          xizs  = (xa(4) - ys*ys*xa(1)) * facsym
          xiyzs =  xa(6) - ys*zs*xa(1)
          if (ksym.ne.0)  xiyzs = 0.0d0
          xad  = 0.5d0*(xiys-xizs)
          if (dabs(xiyzs).lt.tol12 .and. dabs(xad).lt.tol12) then
            a = 0.0d0
          else
            a = datan2(-xiyzs,xad) 
          end if
          a1 =  0.5d0*(xiys+xizs)
          a2 =  0.5d0*(xiys-xizs)
          xieta  = a1 + a2*dcos(a) - xiyzs*dsin(a)
          xizeta = a1 - a2*dcos(a) + xiyzs*dsin(a)
          write(iow,4000)
          write(*  ,4000)
          write(iow,4001) xa(1)*facsym,xa(5),xa(4),xa(6),xiys,xizs,xiyzs
     +                   ,xieta,xizeta,a*90.d0/pi,ys,zs
          write(*  ,4001) xa(1)*facsym,xa(5),xa(4),xa(6),xiys,xizs,xiyzs
     +                   ,xieta,xizeta,a*90.d0/pi,ys,zs
          if (ityp.eq.6) then 
            ystre = xa(16)/xa(7)
            zstre = xa(17)/xa(7)
            if(ksym.eq.1.or.ksym.eq.3)zstre=0.0d0! symmetry on 1-axis
            if(ksym.eq.2.or.ksym.eq.3)ystre=0.0d0! symmetry on 2-axis
            write(iow,4004) 2.d0*xa(7)*facsym,ystre,zstre
            write(*  ,4004) 2.d0*xa(7)*facsym,ystre,zstre
          end if
          return
        end if
c....   print results of second call
        if(.not.flgsec(2) .and.flgsec(3) .and. n.eq.numel) then 
          flgsec(3) = .false.
          if (ityp.eq.6) return
          det =  xizs*xiys - xiyzs*xiyzs
          ym  = ( (xa(8)*xiyzs - xa(9)*xizs ) / det ) * facsym
          zm  = ( (xa(8)*xiys  - xa(9)*xiyzs) / det ) * facsym
          if (ksym.eq.1 .or. ksym.eq.3) zm = 0.0d0   ! symmetry on 1-axis
          if (ksym.eq.2 .or. ksym.eq.3) ym = 0.0d0   ! symmetry on 2-axis
          xa(12) = ym
          xa(13) = zm
          write(iow,4002) ym,zm
          write(*  ,4002) ym,zm
          return
        end if
c....   print results of third call
        if(.not.flgsec(3) .and.flgsec(4) .and. n.eq.numel) then
          flgsec(4) = .false.
          if (ityp.eq.6) return
          ystre = xa(18)/xa(16)/2.d0
          zstre = xa(19)/xa(17)/2.d0
          if(ksym.eq.1.or.ksym.eq.3)zstre=0.0d0! symmetry on 1-axis
          if(ksym.eq.2.or.ksym.eq.3)ystre=0.0d0! symmetry on 2-axis
          write(iow,4003) xa(11)*facsym, xa(10)*facsym, ystre,zstre
          write(*  ,4003) xa(11)*facsym, xa(10)*facsym, ystre,zstre
cww          write(*,*) '  GI_T-ww1',gitww1*facsym
        end if
      else if(ityp.eq.7) then 
        area = d(7)
c....   loop over gauss points
        ll=2
        if(nel.gt.4) ll = 3
        call pgauss(ll,lint,sg,tg,wg)
        do 101 l = 1,lint
          call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
          xsj = xsj*wg(l)
          y   = 0.0d0
          z   = 0.0d0
          wq  = 0.0d0
          wky = 0.0d0
          wkz = 0.0d0
          do ii = 1,nel
            y   = y   + shp(3,ii)*xl(1,ii)
            z   = z   + shp(3,ii)*xl(2,ii)
            wq  = wq  + shp(3,ii)*ul(1,ii)
            wky = wky + shp(1,ii)*ul(1,ii)
            wkz = wkz + shp(2,ii)*ul(1,ii)
          end do       
          if (flgsec(2)) then
c.....      for nue .ne.0
            yy   = y - d(16)      ! coordinates   yq = y - y0
            zz   = z - d(17)
            det = d(8)*d(9)-d(10)*d(10)
            a1  = (d(13)*d(8)-d(14)*d(10))/det      
            a2  = (d(14)*d(9)-d(13)*d(10))/det      
            dt13= -d(15)/(1.d0+d(15))*0.5d0*a1*zz*zz
            dt23=  d(15)/(1.d0+d(15))*0.5d0*a2*yy*yy
            xa(1) = xa(1) + (wky-dt13) * xsj / area
            xa(2) = xa(2) + (wkz+dt23) * xsj / area
            xa(4) = xa(4) +  wq        * xsj / area
          end if
          if(.not. flgsec(2) .and. flgsec(3)) then 
             call sectau(d,xl,ul,shp,ndm,ndf,nel,ityp,x1,x2,t13,t23,
     +                   tmax,xa,ksym)
             xa(3) = xa(3) + (t23*y - t13*z) * xsj
             xa(5) = xa(5) + t13* t13 * xsj
             xa(6) = xa(6) + t23* t23 * xsj                        
          end if
          if(.not. flgsec(3) .and. flgsec(4)) then 
          end if
101     continue
c....   print results of first call
        if(.not.flgsec(1) .and.flgsec(2) .and. n.eq.numel) then
          flgsec(2) = .false.
          return
        end if
c....   print results of second call
        if(.not.flgsec(2) .and.flgsec(3) .and. n.eq.numel) then 
          flgsec(3) = .false.
          ym = 0.d0
          zm = 0.d0
          if (dabs(d(14)).gt.0.d0) ym =  xa(3)/d(14) * facsym
          if (dabs(d(13)).gt.0.d0) zm = -xa(3)/d(13) * facsym
          if (ksym.eq.1 .or. ksym.eq.3) zm = 0.0d0   ! symmetry on 1-axis
          if (ksym.eq.2 .or. ksym.eq.3) ym = 0.0d0   ! symmetry on 2-axis
          write(iow,4002) ym,zm
          write(*  ,4002) ym,zm
c....     shear correction factors
          asy = 0.d0
          if(dabs(d(13)).gt.0.d0) 
     +    asy = d(13)*d(13)/( (xa(5)+xa(6))*facsym*area)
          asz = 0.d0
          if(dabs(d(14)).gt.0.d0) 
     +    asz = d(14)*d(14)/( (xa(5)+xa(6))*facsym*area)
          write(iow,4005) asy,asz
          write(*  ,4005) asy,asz
          return
        end if
c....   print results of third call
        if(.not.flgsec(3) .and.flgsec(4) .and. n.eq.numel) then
          flgsec(4) = .false.
        end if
      else if(ityp.eq.8) then 
cww          flgsec(2) = .false.
cww          flgsec(3) = .false.
cww          flgsec(4) = .false.
cww          xa(1) = 1.d0
C
        area = d(7)
c....   loop over gauss points
        ll=2
        if(nel.gt.4) ll = 3
        call pgauss(ll,lint,sg,tg,wg)
        do 102 l = 1,lint
          call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
          xsj = xsj*wg(l)
          y   = 0.0d0
          z   = 0.0d0
          wq  = 0.0d0
          wky = 0.0d0
          wkz = 0.0d0
          do ii = 1,nel
            y   = y   + shp(3,ii)*xl(1,ii)
            z   = z   + shp(3,ii)*xl(2,ii)
            wq  = wq  + shp(3,ii)*ul(1,ii)
            wky = wky + shp(1,ii)*ul(1,ii)
            wkz = wkz + shp(2,ii)*ul(1,ii)
          end do       
          if (flgsec(2)) then
            xa(1) = xa(1) + wky * xsj / area
            xa(2) = xa(2) + wkz * xsj / area
            xa(4) = xa(4) + wq  * xsj / area
          end if
          if(.not. flgsec(2) .and. flgsec(3)) then 
             call sectau(d,xl,ul,shp,ndm,ndf,nel,ityp,x1,x2,t13,t23,
     +                   tmax,xa,ksym)
             xa(3) = xa(3) + (t23*y - t13*z) * xsj
          end if
          if(.not. flgsec(3) .and. flgsec(4)) then 
          end if
102     continue
c....   print results of first call
        if(.not.flgsec(1) .and.flgsec(2) .and. n.eq.numel) then
          flgsec(2) = .false.
          return
        end if
c....   print results of second call
        if(.not.flgsec(2) .and.flgsec(3) .and. n.eq.numel) then 
          flgsec(3) = .false.
          ym = 0.d0
          zm = 0.d0
          if (dabs(d(14)).gt.0.d0) ym =  xa(3)/d(14) * facsym
          if (dabs(d(13)).gt.0.d0) zm = -xa(3)/d(13) * facsym
          if (ksym.eq.1 .or. ksym.eq.3) zm = 0.0d0   ! symmetry on 1-axis
          if (ksym.eq.2 .or. ksym.eq.3) ym = 0.0d0   ! symmetry on 2-axis
          write(iow,4002) ym,zm
          write(*  ,4002) ym,zm
          return
        end if
c....   print results of third call
        if(.not.flgsec(3) .and.flgsec(4) .and. n.eq.numel) then
          flgsec(4) = .false.
        end if
      end if
      return
4000  format(3x,'S E C T I O N   P R O P E R T I E S',/)
4001  format(3x,'Area of Cross-Section   :',e14.6,/,
     1       3x,'Moment of Inertia       :',
     1       8x,' I_11',11x,' I_22',11x,' I_12',9x,'alpha(ø)',/,
     1       3x,'Referenz Basis          :',3(e16.8)/,
     1       3x,'Center of Gravity  S    :',3(e16.8)/,
     1       3x,'Principal Axis  in S    :',2(e16.8),16x,e16.8,/,
     1       3x,'Center of Gravity 1_S   :',e16.8/,
     1       3x,'                  2_S   :',e16.8/)
4002  format(3x,'Center of Shear   1_M   :',e16.8/,
     1       3x,'                  2_M   :',e16.8/)
4003  format(3x,'St.Venant Stiffness GI_T:',e16.8/,
     1       3x,'Warping Modulus   C_M   :',e16.8//,
     2       3x,'Center of Stress  1_SF  :',e16.8/,
     3       3x,'Function          2_SF  :',e16.8)
4004  format(3x,'St.Venant Modulus I_T   :',e16.8/,
     1       3x,'Center of Stress  1_SF  :',e16.8/,
     2       3x,'Function          2_SF  :',e16.8)
4005  format(3x,'Shear corr.factor asy   :',e16.8/,
     1       3x,'                  asz   :',e16.8/)
4006  format(3x,'eq.23 2nd term:',e16.8)
      end                                       
c
      subroutine stcntors(ix,d,xl,ul,shp,dt,st,ndf,ndm,nel,numnp,
     1                    sg,tg,wg,lint,xa,ksym,ityp)
c----------------------------------------------------------------------+
c     plot warping functions  and shear stresses                       |
c----------------------------------------------------------------------+
      implicit real*8 (a-h,o-z)
      integer*4 ix(*)
      real*8 dt(numnp),st(numnp,*),xl(ndm,*),sg(*),tg(*),wg(*),xa(*)
      real*8 ul(ndf,*),shp(3,*),d(*)
      do 301 l = 1,lint
        call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
        xsj = xsj*wg(l)
        ys  = d(11)
        zs  = d(12)
        call sectau(d,xl,ul,shp,ndm,ndf,nel,ityp,x1,x2,t13,t23,tmax,
     +              xa,ksym)
c....   sum shear stresses for error analysis
        call enerel12(d,t13,t23,xsj)
        do 300 ii = 1,nel
          ll = abs(ix(ii))
          if(ll.gt.0) then
            w        = ul(1,ii)
            xsji     = xsj*shp(3,ii)
            dt(ll)   = dt(ll)   +      xsji
            st(ll,1) = st(ll,1) + t13 *xsji  ! shear stress T_13
            st(ll,2) = st(ll,2) + t23 *xsji  ! shear stress T_23
            st(ll,3) = st(ll,3) + tmax*xsji  ! shear stress T_max
            if (ityp.eq.5) then
              wq  = ul(1,ii)-xa(7)/xa(1)
              if (ksym.ne.0) wq = ul(1,ii)
              wt  = wq+xa(12)*(xl(2,ii)-xa(15))-xa(13)*(xl(1,ii)-xa(14))
              st(ll,4) = st(ll,4) + wq *xsji  ! unit warping  w_q
              st(ll,5) = st(ll,5) + wt *xsji  ! main warping  w_t
            else if (ityp.eq.6) then
              st(ll,4) = st(ll,4) + w  *xsji  ! stress function       
            else if (ityp.eq.7) then           ! unit warping from shear
              facsym1 = 1.0d0
              facsym2 = 1.0d0
              if (ksym.eq.1) then             ! symmetry on 1-axis
                facsym1 = 0.0d0
                facsym2 = 2.0d0
              else if (ksym.eq.2) then         ! symmetry on 2-axis
                facsym1 = 2.0d0
                facsym2 = 0.0d0
              else if (ksym.eq.3) then         ! symmetry on 1/2-axis
                facsym1 = 4.0d0*d(13)
                facsym2 = 4.0d0*d(14)
              end if
              if (dabs(d(14)).gt.0.d0) then
                if (ksym.eq.1) then            ! symmetry on 1-axis
                  facsym1 = 2.0d0
                  facsym2 = 0.0d0
                else if (ksym.eq.2) then        ! symmetry on 2-axis
                  facsym1 = 0.0d0
                  facsym2 = 2.0d0
                end if
              end if
              st(ll,4) = st(ll,4) 
     +                 + (w - xa(4)*facsym2
     +                      - xa(1)*(xl(1,ii)-ys)*facsym1 
     +                      - xa(2)*(xl(2,ii)-zs)*facsym2)*xsji
            end if
          end if
300     continue
301   continue
      return
      end
c
      subroutine sectau(d,xl,ul,shp,ndm,ndf,nel,ityp,x1,x2,t13,t23,tmax,
     +                  xa,ksym)
c----------------------------------------------------------------------+
c     shear stresses from warping function                             |
c----------------------------------------------------------------------+
      implicit double precision(a-h,o-z)
      dimension xl(ndm,*),ul(ndf,*),shp(3,*),d(*),xa(*)
c...  coordinates and derivatives of warping function
      x1  = 0.d0
      x2  = 0.d0
      wk1 = 0.d0
      wk2 = 0.d0
      do i = 1,nel
        x1 = x1   + shp(3,i)*xl(1,i)
        x2 = x2   + shp(3,i)*xl(2,i)
        wk1 = wk1 + shp(1,i)*ul(1,i)
        wk2 = wk2 + shp(2,i)*ul(1,i)
      end do
c...  shear stresses ( assume G = 1, theta = 1)
      if (ityp.eq.5) then
cww        gmod  = 1.d0
        gmod1  = d(1)
        gmod2  = d(2)
        theta = 1.d0
        t13   = gmod1*theta*(-x2+wk1)
        t23   = gmod2*theta*( x1+wk2)
      else if (ityp.eq.6) then
        t13   =  wk2
        t23   = -wk1
      else if (ityp.eq.7 .or. ityp.eq.8) then
C
C>>>>> fuer 8 pruefen
C
        facsym1 = 1.0d0
        facsym2 = 1.0d0
        if (ksym.eq.1) then            ! symmetry on 1-axis
          facsym1 = 0.0d0
          facsym2 = 2.0d0
        else if (ksym.eq.2) then        ! symmetry on 2-axis
          facsym1 = 2.0d0
          facsym2 = 0.0d0
        else if (ksym.eq.3) then        ! symmetry on 1/2-axis
          facsym1 = 4.0d0*d(13)
          facsym2 = 4.0d0*d(14)
        end if
        gmod = 1.d0
        if (ityp.eq.7) then
          if (dabs(d(14)).gt.0.d0) then
            if (ksym.eq.1) then            ! symmetry on 1-axis
              facsym1 = 2.0d0
              facsym2 = 0.0d0
            else if (ksym.eq.2) then        ! symmetry on 2-axis
              facsym1 = 0.0d0
              facsym2 = 2.0d0
            end if
          end if
          t13  = d(13)/d(7) + gmod*( wk1-xa(1)*facsym1 )
          t23  = d(14)/d(7) + gmod*( wk2-xa(2)*facsym2 )
c.....    add for nue .ne.0
          yy   = x1 - d(16)      ! coordinates   yq = y - y0
          zz   = x2 - d(17)
          det = d(8)*d(9)-d(10)*d(10)
          a1  = (d(13)*d(8)-d(14)*d(10))/det      
          a2  = (d(14)*d(9)-d(13)*d(10))/det      
          dt13= -d(15)/(1.d0+d(15))*0.5d0*a1*zz*zz
          dt23=  d(15)/(1.d0+d(15))*0.5d0*a2*yy*yy
          t13  = t13  -  dt13
          t23  = t23  +  dt23
        else if (ityp.eq.8) then
           yy  = x1 - d(11)      ! coordinates   yq = y - ys
           zz  = x2 - d(12)
           det = d(8)*d(9)-d(10)*d(10)
           ay  = (d(13)*d(8)-d(14)*d(10))/2.d0/det      
           az  = (d(14)*d(9)-d(13)*d(10))/2.d0/det      
           t13  =  wk2  ! for dtau  
           t23  = -wk1  !           
cww           t13  =  wk2 - ay*yy*yy ! for tau    ?*facsym1  
cww           t23  = -wk1 - az*zz*zz !            ?*facsym2  
        end if
      end if
      tmax  = dsqrt(t13*t13 + t23*t23)
      return
      end
c
      subroutine enerel12(d,q1,q2,xsj)
c----------------------------------------------------------------------
c.....sum energy for error analysis
c----------------------------------------------------------------------
      USE errin1
      implicit double precision (a-h,o-z)
      dimension d(*)
      u_om(1) = u_om(1) + (q1*q1/d(1) + q2*q2/d(2))*xsj
      u_om(2) = u_om(2) + (q1*q1      + q2*q2     )*xsj
      return
      end
c
      subroutine stcn12e(ix,d,xl,ul,shp,dt,st,ndf,ndm,nel,numnp,numel,
     1                    sg,tg,wg,lint,ityp,e_ome,xa,ksym)
c----------------------------------------------------------------------+
c     errors from energy of flow/stress differences                    |
c----------------------------------------------------------------------+
      USE errin1
      USE errin2
      implicit real*8 (a-h,o-z)
      integer*4 ix(*)
      real*8 dt(numnp),st(numnp,*),xl(ndm,*),sg(*),tg(*),wg(*)
      real*8 ul(ndf,*),shp(3,*),d(*),xa(*)
      real*8 e_ome(*),e_ome12(numerr)
c.... sum element energy
c.... intial values for element errsors
      e_ome12 = 0.0d0
c
      do 301 l = 1,lint   ! Gauss-Points
c.....  Gauss Point flows from mat. law
        call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
        xsj = xsj*wg(l)
        if(ityp.le.4) then
          call flolag(d,shp,xl,ul,rr,zz,q1,q2,qm,pr,ndm,ndf,ityp)
        else if(ityp.ge.5) then
          call sectau(d,xl,ul,shp,ndm,ndf,nel,ityp,x1,x2,q1,q2,tmax,
     +                xa,ksym)
        end if
c.....  Gauss Point flows from nodal interpolation
        q1n = 0.0d0
        q2n = 0.0d0 
        do 300 ii = 1,nel
          ll = abs(ix(ii))
          if(ll.gt.0) then
            xsji = shp(3,ii)
            q1n = q1n + st(ll,1) *xsji  ! flow_x
            q2n = q2n + st(ll,2) *xsji  ! flow_y
          end if
300     continue
c.....  Diff. flows  at gauss point
        dq1 = q1n - q1
        dq2 = q2n - q2
c.....  element errors
        e_ome12(1)= e_ome12(1) + (dq1*dq1/d(1) + dq2*dq2/d(2))*xsj
        e_ome12(2)= e_ome12(2) + (dq1*dq1      + dq2*dq2     )*xsj
301   continue
c.... plot/print errors
      call elmterr(ix,xl,ndm,numel,e_ome12,e_ome)
      return
      end
c
