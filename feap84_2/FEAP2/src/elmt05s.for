      subroutine elmt05(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c--------------------------------------------------------------------
c.... 3-9 plane/axisymmetric linear elastic element routine 
c--------------------------------------------------------------------
c     Input: E, nu, rho, ityp,  [pts/stiff, pts/stre]
c
c     d-field
c     d( 1) = e/(1-nu^2) or e(1-nu)/((1+nu)(1-nu)] = E_11
c     d( 2) = nu*E_11
c     d( 3) = e/[2*(1-nu)] = E_33 = G
c     d( 4) = rho = gamma/g [F/L**3]/ [S/T**2]
c     d( 5) = l (int.points for K)
c     d( 6) = k (int.points for Sigma)
c     d( 7) = Y0 for error analysis (not documented) 
c     d( 8) =
c     d( 9) = Base temperature T_o or Temp difference T_o-T_1
c     d(10) = e*alpha_t/(1-nu) or e*alpha_t/(1-2nu)
c     d(11) = f_x  [F/L**3] 
c     d(12) = f_y  [F/L**3]  
c     d(13) = alpha_t
c     d(14) = thick h
c     d(15) = ityp 0,1=ESZ,2=EVZ,3=AXIS
c     d(16) = e
c     d(17) = nu
c     d(18) = -nu/e
c--------------------------------------------------------------------
c     q(01) = fx  [F/L**2] !!                                       
c     q(02) = fy  [F/L**2] !!
c--------------------------------------------------------------------
c     error analysis added
c     elmt no. used for plot, see inord, ipord
c     term 2*pi added for axisymm.
c     
c--------------------------------------------------------------------
      USE bdata
      USE cdata
      USE eldata
      USE errin1
      USE errin2
      USE errin3
      USE iofile
      USE pdata6
      USE pdata10
      USE qload
      USE strnam
      implicit double precision (a-h,o-z)
      dimension   xl(ndm,*),tl(*),sig(6),ix(*),
     1       d(*),ul(ndf,*),s(nst,*),p(*),eps(4),sigr(4),shp(3,9),
     1       sg(16),tg(16),wg(16),di(7),ipord3(4),ipord6(7),vl(ndf*nel)
      dimension h1(*),h2(*),h3(*)
      character wd(3)*12
      data wd/'Plane Stress','Plane Strain','Axisymmetric'/
      data ipord3 /1,2,3,1/
      data ipord6 /1,4,2,5,3,6,1/
      ielno = 5
      pi = 4.0d0*datan(1.0d0)
c.... go to correct array processor
      go to(1,2,3,4,5,4,2,8,9,2,2,2,2,14,2,2,2,2,2,2,2,3), isw
      return
c.... input material properties
1     if(ior.lt.0) write(*,5000)
      call dinput(di,7)
      e    = di(1)
      xnu  = di(2)
      d(4) = di(3)
      ityp = di(4)
      l    = di(5)
      k    = di(6)
      d(7) = di(7) ! y0
      if(ior.lt.0) write(*,5001)
      call dinput(di,5)
      d(14)= di(1)  ! thickness
      d(11)= di(2)  ! f_x
      d(12)= di(3)  ! f_y
      alp  = di(4)  
      t0   = di(5)
c.... set material parameter type and flags
      ityp = max(1,min(ityp,3))
      j = min(ityp,2)
      d(1) = e*(1.+(1-j)*xnu)/(1.+xnu)/(1.-j*xnu)
      d(2) = xnu*d(1)/(1.+(1-j)*xnu)
      d(3) = e/2./(1.+xnu)
      if(d(14).le.0.0d0 .or. ityp.eq.2) d(14) = 1.0 ! thickness: ityp=1->h, ityp=2 always 1 
      if(ityp.eq.3) d(14)=0.d0
      d(15) = ityp
      d(16) = e
      d(17) = xnu
      d(18) = -xnu/e
cww   l = min(4,max(1,l))
cww   k = min(4,max(1,k))
c.... set gauss points from input or direct
      if(l.eq.0) then        
        if(nen.le.4) then
          l=2
        else
          l=3
        end if
      end if
      if(k.eq.0) then        
        if(nen.le.4) then
          k=1
        else
          k=2
        end if
      end if
      d(5) = l    ! Gausspoints for K 
      d(6) = k    ! Gausspoints for Stre
      d(9) = t0   ! Temp-diff
      d(10)= e*alp/(1.-j*xnu)
      d(13) = alp ! alpha_t 
                   write(iow,2000) wd(ityp),d(16),d(17),d(4),l,k,d(14),
     1                                      d(11),d(12),alp,t0
      if(ior.lt.0) write(  *,2000) wd(ityp),d(16),d(17),d(4),l,k,d(14),
     1                                      d(11),d(12),alp,t0
c.... define node numbering for plot mesh routine, see pltord
      if(nen.eq.3) then
        inord(ielno) = 4
        do ii = 1,4
          ipord(ii,ielno) = ipord3(ii)
        end do
      else if(nen.eq.6) then
        inord(ielno) = 7
        do ii = 1,7
          ipord(ii,ielno) = ipord6(ii)
        end do
      end if
c.... description of stresses  
      strsus( 1) =  '  Stress  S_xx '
      strsus( 2) =  '  Stress  S_xy '
      strsus( 3) =  '  Stress  S_yy '
      strsus( 4) =  '  Stress  S_zz '
      if(ityp.eq.1) strsus( 4) =  '  Stress  S_v  '
      if(ityp.eq.3) strsus( 4) =  '  Stress  S_pp ' 
      strsus( 5) =  '  Stress  S_1  '
      strsus( 6) =  '  Stress  S_2  '
      strsus( 7) =  '  ANGLE phi_1  '
      strsus( 8) =  '  Eps     E_xx '
      strsus( 9) =  '  Eps     E_xy '
      strsus(10) =  '  Eps     E_yy '
      strsus(11) =  '               '
      if(ityp.eq.1) strsus(11) =  '  Eps     E_zz '
      if(ityp.eq.3) strsus(11) =  '  Eps     E_pp ' 
      do is =12,25
        strsus(is) = '               '
      end do
2     return
c.... stiffness/residual computation
3     l    = d(5)
      ityp = d(15)
      call pgauss(l,lint,sg,tg,wg)
      d1 = d(10)*d(14)
c.... compute integrals of shape functions
      do 340 l = 1,lint
        call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
        xsj = xsj*wg(l)
c....   compute temperature field: with input via TEMP/BTEM: d(9)=T_0, without d(9)=-dt=T_o-T_1
        ta = -d(9)
        xx = 0. 
        do 300 i = 1,nel
          ta = ta + shp(3,i)*tl(  i)
          xx = xx + shp(3,i)*xl(1,i)
300     continue
        ta = d1*ta
        if(ityp.le.2) then
          dv  = xsj
          xsj = 0.0
        else
          xsj = xsj*2.0d0*pi ! ww
          dv  = xsj*xx
        end if
c....   loop over rows
        j1 = 1
        do 330 j = 1,nel
          w11 = shp(1,j)*dv
          w12 = shp(2,j)*dv
          w22 = shp(3,j)*xsj
c....     compute the gravity and thermal loads
          call qload05(qxg,qyg,d,aqloa,numel,n,mqloa,propq,isw)

          p(j1  ) = p(j1  ) + qxg*shp(3,j)*dv +(w11+w22)*ta
          p(j1+1) = p(j1+1) + qyg*shp(3,j)*dv + w12*ta
          if(isw.eq.22) goto 331
          
c....     loop over columns (symmetry noted)
          k1 = j1
c....     plane problem (fast computation form)
          if(ityp.le.2) then
            do 310 k = j,nel
              s(j1  ,k1  ) = s(j1  ,k1  ) + w11*shp(1,k)
              s(j1  ,k1+1) = s(j1  ,k1+1) + w11*shp(2,k)
              s(j1+1,k1  ) = s(j1+1,k1  ) + w12*shp(1,k)
              s(j1+1,k1+1) = s(j1+1,k1+1) + w12*shp(2,k)
              k1 = k1 + ndf
310         continue
          else
c....       axisymmetric problem
            a11 = d(1)* w11 + d(2)*w22
            a21 = d(2)* w11 + d(1)*w22
            a31 = d(2)*(w11+w22)
            a41 = d(3)* w12
            a12 = d(2)* w12
            a32 = d(1)* w12
            a42 = d(3)* w11
            do 320 k = j,nel
              w11 = shp(1,k)
              w12 = shp(2,k)
              w22 = shp(3,k)/xx
              s(j1  ,k1  ) = s(j1  ,k1  ) +  w11*a11 + w22*a21 + w12*a41
              s(j1+1,k1  ) = s(j1+1,k1  ) + (w11 + w22)*a12 + w12*a42
              s(j1  ,k1+1) = s(j1  ,k1+1) +  w12*a31 + w11*a41
              s(j1+1,k1+1) = s(j1+1,k1+1) +  w12*a32 + w11*a42
              k1 = k1 + ndf
320         continue
          end if
331     j1 = j1 + ndf
330     continue
340   continue
      if(isw.eq.22) return
c.... assemble plane stiffness matrix from integrals and material props.
      if(ityp.le.2) then
        d1 =  d(1)*d(14)
        d2 =  d(2)*d(14)
        d3 =  d(3)*d(14)
        nsl = nel*ndf
        do 350 j = 1,nsl,ndf
          do 350 k = j,nsl,ndf
            w11 = s(j,k)
            w12 = s(j,k+1)
            w21 = s(j+1,k)
            w22 = s(j+1,k+1)
            s(j  ,k  ) = d1*w11 + d3*w22
            s(j  ,k+1) = d2*w12 + d3*w21
            s(j+1,k  ) = d2*w21 + d3*w12
            s(j+1,k+1) = d1*w22 + d3*w11
350       continue
      end if
c.... make stiffness symmetric
      do 360 j = 1,nst
        do 360 k = j,nst
          s(k,j) = s(j,k)
360   continue
c.... compute a residual
      call pmove(ul,vl,ndf*nel)
      do 370 j = 1,nst
        do 370 k = 1,nst
          p(j) = p(j) - s(j,k)*vl(k)
370   continue
c.... calculate eigenvalues
cww   call elemev(s,nel,ndf,nst)
c
      return
c.... output element variables  isw 4,6
4     l    = d(5)
      k    = d(6)
      ityp = d(15)
      if(isw.eq.4) l = k
      call pgauss(l,lint,sg,tg,wg)
c.... compute element stresses, strains, and forces
      do 440 l = 1,lint
c....   compute element shape functions
        call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
c....   compute stresses, strains and coordinates
        call stre05(d,xl,ul,tl,shp,eps,sigr,sigv,xx,yy,ndm,ndf,nel,ityp)
        if(isw.eq.4) then
          sig(1) = sigr(1)
          sig(2) = sigr(2)
          sig(3) = sigr(3)
          sig(4) = sigr(4)
          call pstres(sig,sig(5),sig(6),ang)
c....     output stresses and strains
          mct = mct - 2
          if(mct.le.0) then
            if(ityp.eq.1) then ! with sigv
                           write(iow,2003) o,head
              if(ior.lt.0) write(*  ,2003) o,head
            else 
                           write(iow,2001) o,head
              if(ior.lt.0) write(*  ,2001) o,head
            end if
            mct = 50
          end if
          if(ityp.eq.1) sig(4) = sigv  
                       write(iow,2002) n,xx,sig,ma,yy,eps,ang
          if(ior.lt.0) write(*  ,2002) n,xx,sig,ma,yy,eps,ang
        else if(isw.eq.6) then
c....     compute internal forces
          if(ityp.le.2) then
            dv      = xsj*wg(l)*d(14)
            sigr(4) = -d(11)*dv
          else
            xsj     = xsj*wg(l)
            xsj     = xsj*2.0d0*pi !ww
            dv      = xsj*xx 
            sigr(4) = sigr(4)*xsj - d(11)*dv
          end if
          j1 = 1
          do 430 j = 1,nel
            p(j1)   = p(j1)   - (shp(1,j)*sigr(1)+shp(2,j)*sigr(2))*dv
     1                        -  shp(3,j)*sigr(4)
            p(j1+1) = p(j1+1) - (shp(1,j)*sigr(2)+shp(2,j)*sigr(3))*dv
     1                        + d(12)*shp(3,j)*dv
            j1 = j1 + ndf
430       continue
        end if
440   continue
      return
c.... compute lumped mass matrix
5     l    = d(5)
      ityp = d(15)
      call pgauss(l,lint,sg,tg,wg)
      do 520 l = 1,lint
c....   compute shape functions
        call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.true.) 
        dv = wg(l)*xsj*d(4)*d(14)
        if(ityp.gt.2) then
          xx = 0.
          do 500 i = 1,nel
            xx = xx + shp(3,i)*xl(1,i)
500       continue
          dv = dv*xx
          dv = dv*2.0d0*pi ! ww
        end if
c....   for each node j compute db = rho*shape*dv*thickness
        j1 = 1
        do 510 j = 1,nel
          p(j1  ) = p(j1) + shp(3,j)*dv
          p(j1+1) = p(j1)
          j1 = j1 + ndf
510     continue
520   continue
      return
c.... compute the nodal stress values 
8     istv = 11
      if(iplma(ma).eq.0)  return ! only if MATN
cww   l    = d(5)
      l    = d(6)
      ityp = d(15)
      if(nel.eq.8) then
        l=3
        call pnewcot(l,lint,sg,tg,wg) ! l=3 Int.Points=8nodes+Center
      else
        call pgauss(l,lint,sg,tg,wg)
      end if

      call stcn05(ix,d,xl,ul,tl,shp,strea,strea(1+numnp),ndf,ndm,nel,
     1   numnp,numel,sg,tg,wg,sigr,eps,lint,ityp)
      return
c.... compute the stress errors
cww9  l    = d(5)
9     l    = d(6)
      ityp = d(15)
      call pgauss(l,lint,sg,tg,wg)
      call ster05(ix,d,xl,ul,tl,shp,strea,strea(1+numnp),ndf,ndm,
     1            numnp,numel,sg,tg,wg,sigr,eps,lint,ityp,e_ome)
      return
c.... calculate stresses sig(i) at center of element
14    l    = 1
      k    = 1
      ityp = d(15)
      call pgauss(l,lint,sg,tg,wg)
      call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
      call stre05(d,xl,ul,tl,shp,eps,sigr,sigv,xx,yy,ndm,ndf,nel,ityp)
      sig(1) = sigr(1)
      sig(2) = sigr(2)
      sig(3) = sigr(3)
      sig(4) = sigr(4)
      if(ityp.eq.1) sig(4) = sigv  
      call pstres(sig,sig(5),sig(6),ang)
c
      if(nfp.gt.11.or.nfp.lt.1) return
      if(nfp.ge.1.and.nfp.le.6) strp=sig(nfp)
      if(nfp.eq.7)              strp=ang
      if(nfp.ge.8.and.nfp.le.11) strp=eps(nfp-7)
c
      if(flfp) then
c....   calculate extreme values
        xmaxf = max(xmaxf,strp)
        xminf = min(xminf,strp)
      else
c....   plot stresses 
        call pppcolf(strp) 
        call plot9(iel,ix,xl,ndm,nel,1)
      end if
      return
c
c.... formats for input-output
2000  format(/5x,a12,' Linear Elastic Element'//
     1 10x,'Modulus',e18.5/10x,'Poisson ratio',f8.5/10x,
     2 'Density',e18.5/10x,'Gauss pts/dir',i3/10x,'Stress pts',i6/
     3 10x,'Thickness',e16.5/10x,'1-gravity',e16.5/10x,'2-gravity',
     4 e16.5/10x,'Alpha',e20.5/10x,
     5 'Base temp T_0 or Temp diff. T_0-T_1',e16.5/)
2001  format(a1,20a4//5x,'Element Stresses'//' elmt  1-coord',
     1    2x,'11-stress',2x,'12-stress',2x,'22-stress',2x,'33-stress',
     2    3x,'1-stress',3x,'2-stress'/' matl  2-coord',2x,'11-strain',
     3    2x,'12-strain',2x,'22-strain',2x,'33-strain',6x,'angle'/
     4    39(' -'))
2003  format(a1,20a4//5x,'Element Stresses'//' elmt  1-coord',           ! with sigv
     1    2x,'11-stress',2x,'12-stress',2x,'22-stress',2x,' v-stress',
     2    3x,'1-stress',3x,'2-stress'/' matl  2-coord',2x,'11-strain',
     3    2x,'12-strain',2x,'22-strain',2x,'33-strain',6x,'angle'/
     4    39(' -'))
2002  format(i5,0p1f9.3,1p6e11.3/i5,0p1f9.3,1p4e11.3,0p1f11.2/)
cww5000  format(' Input: E, nu, rho, pts/stiff, pts/stre',
5000  format(' Input: E, nu, rho',
     1 ', type(1=stress,2=strain,3=axism)',/)
5001  format(' Input: Thickness, 1-body force, 1-body force, alpha,'
     1       ,' Temp-base T_0 or Temp-diff T_0-T_1'/)
      end
c
      subroutine stre05(d,xl,ul,tl,shp,eps,sig,sigv,xx,yy,ndm,ndf,nel,
     +                  ityp)
c----------------------------------------------------------------------
c.... stress calculation                                              |
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension   xl(ndm,*),tl(*),d(*),ul(ndf,*),shp(3,*),eps(4),sig(4)
c.... compute strains and coordinates
      call pzero(eps,4)
      xx = 0.0
      yy = 0.0
      ta = -d(9)
      do 100 j = 1,nel
        xx = xx + shp(3,j)*xl(1,j)
        yy = yy + shp(3,j)*xl(2,j)
        ta = ta + shp(3,j)*tl(  j) ! Temp from tl-field
        eps(1) = eps(1) + shp(1,j)*ul(1,j)
        eps(2) = eps(2) + shp(1,j)*ul(2,j) + shp(2,j)*ul(1,j)
        eps(3) = eps(3) + shp(2,j)*ul(2,j)
        if(ityp.eq.3) eps(4) = eps(4) + shp(3,j)*ul(1,j)
100   continue
      tda = ta*d(10)
c.... compute stresses
c.... stress - circumferential 
      if(ityp.eq.3) then
c....   strain - circumferential 
        if(xx.ne.0.0d0) then
          eps(4) = eps(4)/xx 
        else
          eps(4) = eps(1)
        end if			
c....   stress - circumferential 
        sig(4) = d(1)*eps(4) + d(2)*(eps(1) + eps(3)) - tda
      else if(ityp.eq.2) then
c....   stress - EVZ 
        sig(4) = d(2)*(eps(1) + eps(3)) - tda
      else
c....   stress - ESZ 
        sig(4) = 0.0
      end if
c.... stresses
      sig(1) = d(1)*eps(1) + d(2)*(eps(3) + eps(4)) - tda 
      sig(2) = d(3)*eps(2)
      sig(3) = d(1)*eps(3) + d(2)*(eps(1) + eps(4)) - tda
      if(ityp.eq.1) eps(4) = d(18)*(sig(1) + sig(3)) + ta*d(13)
c     sigma_v
      sigv= 0.d0  
      if(ityp.eq.1) sigv = sqrt(
     + sig(1)*sig(1)+sig(3)*sig(3)-sig(1)*sig(3)+3.d0*sig(2)*sig(2))
c
      return
      end
c
      subroutine stcn05(ix,d,xl,ul,tl,shp,dt,st,ndf,ndm,nel,numnp,numel,
     1                sg,tg,wg,sig,eps,lint,ityp)
c ----------------------------------------------------------------------
c.... stress projection  and energy                                    |
c ----------------------------------------------------------------------
      USE errin1
      USE errin2
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),xl(ndm,*),tl(*),ix(*),
     1          eps(4),sig(4),ul(ndf,*),shp(3,9),d(*),sg(9),tg(9),wg(9)
      gr = (1.+d(17))/d(16)
      do 130 l = 1,lint
        call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
        call stre05(d,xl,ul,tl,shp,eps,sig,sigv,xx,yy,ndm,ndf,nel,ityp)

c....   additional terms element energy
        s4=sig(4)
        if(ityp.eq.1) s4=0.d0
        u_om(1)=u_om(1)
     +         +(d(18)*(sig(1)+sig(3)+s4)**2+gr*sig(2)**2)*xsj
c....   element energy
        ntyp = 4
        if(ityp.eq.1) ntyp = 3
        do 110 i = 1,ntyp
          u_om(1) = u_om(1) + gr*sig(i)**2*xsj
          u_om(2) = u_om(2) +    sig(i)**2*xsj
110     continue
        y0=d(7)
        u_om(3) = u_om(3) + y0*y0
c....   stress projection
        do 120 ii = 1,nel
          xsji = xsj*shp(3,ii)*wg(l)
          if(nel.eq.8) xsji = xsj*shp(3,ii) ! no weight,most on mid-side!    
          ll = abs(ix(ii))
          if(ll.gt.0) then
            dt(ll) = dt(ll) + xsji
            st(ll, 1) = st(ll, 1) + sig(1)*xsji
            st(ll, 2) = st(ll, 2) + sig(2)*xsji
            st(ll, 3) = st(ll, 3) + sig(3)*xsji
            st(ll, 4) = st(ll, 4) + sig(4)*xsji
            if(ityp.eq.1) st(ll, 4) = st(ll, 4) + sigv*xsji
            st(ll, 8) = st(ll, 8) + eps(1)*xsji
            st(ll, 9) = st(ll, 9) + eps(2)*xsji
            st(ll,10) = st(ll,10) + eps(3)*xsji
            if(ityp.eq.1) st(ll,11) = st(ll,11) + eps(4)*xsji
            if(ityp.eq.3) st(ll,11) = st(ll,11) + eps(4)*xsji
          end if
120     continue
130   continue
      return
      end
c
      subroutine ster05(ix,d,xl,ul,tl,shp,dt,st,ndf,ndm,
     1   numnp,numel,sg,tg,wg,sig,eps,lint,ityp,e_ome)
c ----------------------------------------------------------------------
c.... error calculation                                                |
c ----------------------------------------------------------------------
      USE eldata
      USE errin1
      USE errin2
      USE iofile
      USE pdata2
      implicit double precision (a-h,o-z)
      dimension   dt(numnp),st(numnp,*),xl(ndm,*),tl(*),
     1 ix(*),shp(3,4),sig(6),sigp(4),dsig(4),d(*),eps(4),ul(ndf,*),
     2 sg(16),tg(16),wg(16),e_ome(numel,*),e_ome05(numerr)
c.... intial values for element errors
      e_ome05 = 0.0d0
      gr = (1.+d(17))/d(16)
      do 200 ii = 1,lint
c....   stresses at Gauss points
        call shape(sg(ii),tg(ii),xl,shp,xsj,ndm,nel,ix,.false.)
        xsj = xsj*wg(ii)
        call stre05(d,xl,ul,tl,shp,eps,sig,sigv,xx,yy,ndm,ndf,nel,ityp)
c....   stresses at nodes
        call pzero(sigp,4)
        do 110 i = 1,nel
          ll = iabs(ix(i))
          if(ll.ne.0) then
            do 105 j = 1,4
              sigp(j) = sigp(j) + shp(3,i)*st(ll,j)
105         continue
          end if
110     continue
c.....  stress differences Sigma_Node-Sigma_GP  
        do i = 1,4
          dsig(i) = sigp(i)- sig(i)
        end do 
c.....  element energy from differences 1,2
        ntyp=4
        if(ityp.eq.1) ntyp=3
        do i = 1,ntyp
          e_ome05(1)= e_ome05(1) + gr*(dsig(i)**2)*xsj
          e_ome05(2)= e_ome05(2) +    (dsig(i)**2)*xsj
        end do 
c....   additional term
        s4=dsig(4)
        if(ityp.eq.1) s4=0.d0
        e_ome05(1)= e_ome05(1) 
     +         +(d(18)*(dsig(1)+dsig(3)+s4)**2+gr*dsig(2)**2)*xsj

c.....  element 'energy' = Sigma_V  3
c       e_ome05(3)= e_ome05(3) + sig(3)*sig(3) !######## sigv
        e_ome05(3)= e_ome05(3) + sigv*sigv 

200   continue
c.....plot/print/add errors
      call elmterr(ix,xl,ndm,numel,e_ome05,e_ome)
      return
      end
c
      subroutine qload05(qxg,qyg,d,q,numel,n,mqloa,propq,isw)
c----------------------------------------------------------
c.... define gravity loads
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*)
      if(isw.eq.22) then
        qxg = 0.d0 
        qyg = 0.d0 
        if(mqloa.ne.1) then 
          qxg = q(n,1)*propq
          qyg = q(n,2)*propq
        end if 
      else
        qxg = d(11) 
        qyg = d(12) 
      end if
      return
      end
