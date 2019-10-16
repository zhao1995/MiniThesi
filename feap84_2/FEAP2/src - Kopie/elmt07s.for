      subroutine elmt07(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-------------------------------------------------------+
c     quadrilateral linear DKQ plate element - 4 - node |
c     modified for all ndf and nel                      |
c-------------------------------------------------------+
c     open 
c     Betonbemessung ist noch nicht richtig!!!
c     Einbau nur vorläufig für m_x+_m_xy
c
c     Einbau Querkräfte unter isw 4,8
c     Kragarm Einzellast ok
c     Kragarm unter Gleichast ok
c     Quadratplatte ok aber deutlich zu klein
c
c-------------------------------------------------------+
c                                                       |
c     d(1)  = Eh**3/(1-nu**2)/12                        |
c     d(2)  = nu                                        |
c     d(3)  = (1-nu)/2*Eh**3/(1-nu**2)/12               |
c     d(4)  =                                           |
c     d(5)  =                                           |
c     d(6)  = q                                         |
c     d(7)  = rho                                       |
c     d(8)  = rho*h                                     |
c     d(9)  = h**2/12                                   |
c     d(10) = c                                         |
c     d(11) = f_cd                                      |
c     d(12) = f_yd                                      |
c     d(13) = d_x                                       |
c     d(14) = d_y                                       |
c     d(15) = h                                         |
c     d(16) = -a_t*dt/h                                 |
c     d(17) =                                           |
c     d(18) =                                           |
c-------------------------------------------------------+
c     q(01) = q                                         |
c     q(02) = dt                                        |
c-------------------------------------------------------+
      USE bdata
      USE cdata
      USE eldata
      USE errin2
      USE errin3
      USE iofile
      USE mdat2
      USE pdata7
      USE pdata10
      USE pltran
      USE prisdat
      USE qload
      USE strnam
      implicit double precision (a-h,o-z)
      dimension xl(ndm,*),tl(*),ix(*),di(10),sl(12,12),sigp(7),
     1        d(*),ul(*),s(nst,*),p(nst),shp(3,8),yl(3,4),
     1        sg(4),tg(4),wg(4),sig(5),shpw(3,4),as(4),ptemp(12),
     1        sigq(4,3),signo(4,3)
      dimension h1(*),h2(*),h3(*)
      character comp*25,text*15
      common /elcom071/  b(4),c(4),aa(4),bb(4),cc(4),dd(4),ee(4)
!$OMP THREADPRIVATE (/elcom071/)  
      common /elcom072/  bm(3,12)
!$OMP THREADPRIVATE (/elcom072/)  
c.... go to correct processor
      go to (1,2,3,4,5,3,7,8,9,2,2,2,2,14,2,2,2,2,2,2,2,22), isw
      return
c.... input material properties
1     if(ior.lt.0) write(*,5000)
      call dinput(di,10)
c
      emod  = di(1)
      xnu   = di(2)
      thick = di(3)
      d(6)  = di(4)
      d(7)  = di(5)
      d(10) = di(6) 
      d(11) = di(7)
      d(12) = di(8)
      d(13) = di(9)
      d(14) = di(10)
      d(15) = thick
       
      d(1) = di(1)/(1.-di(2)*di(2))*thick**3/12.0
      d(2) = di(2)*d(1)
      d(3) = 0.5*(d(1)-d(2))

      d(8) = di(5)*thick
      d(9) = thick*thick/12.

      if(ior.lt.0) write(*,5001)
      call dinput(di,2)

      at  = di(1)
      dt  = di(2)
      d(16)=-at*dt/thick

      if(d(10).ge.0.d0)then
         comp='(tension/compression)'
      else
         comp='(  only  compression)'
      end if   
c
                   write(iow,2000) emod,xnu,thick,d(6),d(7),d(10),comp
      if(ior.lt.0) write(*  ,2000) emod,xnu,thick,d(6),d(7),d(10),comp
      if(d(11).gt.0.0d0) then
                     write(iow,2001) d(11),d(12),d(13),d(14)
        if(ior.lt.0) write(*  ,2001) d(11),d(12),d(13),d(14)
      end if
      if(at.gt.0.0d0) write(iow,2006) at,dt
      if(at.gt.0.0d0.and.ior.lt.0) write(*,2006) at,dt
  
      if(ndm.eq.3.and.ndf.eq.3) ipla = 1
c.... description of stresses  
      strsus( 1) = '  MOMENT m_11  '
      strsus( 2) = '  MOMENT m_12  '
      strsus( 3) = '  MOMENT m_22  '
      strsus( 4) = 'FOUND. PRESSURE'
      strsus( 5) = '  MOMENT m_1   '
      strsus( 6) = '  MOMENT m_2   '
      strsus( 7) = '  ANGLE Phi_1  '
      strsus( 8) = '  S-FORCE q_1  '
      strsus( 9) = '  S-FORCE q_2  '                
      strsus(10) = 'STEEL as_x(Bot)'
      strsus(11) = 'STEEL as_x(Top)'                
      strsus(12) = 'STEEL as_y(Bot)'                
      strsus(13) = 'STEEL as_y(Top)'                
      do is =14,25
        strsus(is) = '               '
      end do
c...  names for principal moments
      nptyp = 3 
c
2     return
c.... compute the element tangent array
3     cb = d(10)
      call jacqudpl(xl,ndm)
      q  = d(6)
      d6 = q/4.0
      l  = 2
      call pgauss(l,lint,sg,tg,wg)
      call pzero(sl,144)
      temp = (d(1)+d(2))*d(16)
      do 310 l = 1,lint
        call qushp8pl(sg(l),tg(l),shp,xl,ndm,xsj)
        xsj = xsj*wg(l)
        call dktqbmpl(shp(1,5),shp,ndf)
c....   compute weighted jacobian and d-matrix constants
        d1 = d(1)*xsj
        d2 = d(2)*xsj
        d3 = d(3)*xsj
        d6x= d6  *xsj
c....   compute the element load vector
        p(      1) = d6x*(1.-sg(l))*(1.-tg(l)) + p(      1)
        p(  ndf+1) = d6x*(1.+sg(l))*(1.-tg(l)) + p(  ndf+1)
        p(2*ndf+1) = d6x*(1.+sg(l))*(1.+tg(l)) + p(2*ndf+1)
        p(3*ndf+1) = d6x*(1.-sg(l))*(1.+tg(l)) + p(3*ndf+1)
c....   compute element load vector due to temperature
        do it = 1,12
          ptemp(it) = (bm(1,it)+bm(2,it)) * temp * xsj
        end do
c....   copy into load vector
        kt = 0
        do it = 1,4
          jt=it-1
          p(jt*ndf+1) = ptemp(kt+1) + p(jt*ndf+1)
          p(jt*ndf+2) = ptemp(kt+2) + p(jt*ndf+2)
          p(jt*ndf+3) = ptemp(kt+3) + p(jt*ndf+3)
          kt = kt+3
        end do

c....   compute contribution to element stiffness for this point
        do 300 i = 1,12
          dn1 = d1*bm(1,i) + d2*bm(2,i)
          dn2 = d2*bm(1,i) + d1*bm(2,i)
          dn3 = d3*bm(3,i)
          do 300 j = i,12
            sl(i,j) = sl(i,j) + dn1*bm(1,j) + dn2*bm(2,j) + dn3*bm(3,j)
300     continue
c....   elastic foundation
        if(dabs(cb).gt.0.0d0.and.isw.eq.3) then
          call shapef(sg(l),tg(l),xl,shpw,xsjw,ndm,.false.)
c.....    calculate pressure from actual displacements (not with time!)
          press = 0.d0
          pressw= 0.d0
          do i = 1,4
            pw = ul(nen*ndf+ndf*(i-1)+1)
            pressw = pressw +pw*pw
            press =  press+dabs(cb)*shpw(3,i)*ul(nen*ndf+ndf*(i-1)+1)
c....       geht das bei ndf =6??
          end do     
          if(pressw.eq.0.d0) press = -1.d0 ! only for first step, to have springs 
c....     add term always/only for negative values
          if(cb.gt.0.d0.or.(cb.lt.0.d0.and.press.lt.0.0d0)) then
            efo1 = dabs(cb)*xsjw*wg(l)
            ii = 1
            do 320 i = 1,4
              elfo =  efo1*shpw(3,i)
              jj = ii
              do 321 j = i,4
                s(ii,jj) = s(ii,jj) + elfo*shpw(3,j)
321           jj = jj + ndf
320         ii = ii + ndf
          end if
        end if  
310   continue
c.... copy sl into stiffness matrix (blocs 3*3)
      i1 = 0
      i2 = 0
      do i = 1,4
        j1 = i1
        j2 = i2
        do j = i,4 
          do ii =1,3
            do jj = 1,3
              s(i1+ii,j1+jj) = s(i1+ii,j1+jj) + sl(i2+ii,j2+jj)
            end do
          end do
         j1 = j1 + ndf
         j2 = j2 + 3
        end do
        i1 = i1 + ndf
        i2 = i2 + 3
      end do
c.... make stiffness symmetric
      do 330 i = 2,nst
        do 330 j = 1,i-1
          s(i,j) = s(j,i)
330   continue
c.... compute the element residual vector: R = P - K*v
      do 350 j = 1,nst
        uu = ul(j)
        if(uu.ne.0.0d0) then
          do 340 i = 1,nst
            p(i) = p(i) - s(i,j)*uu
340       continue
        end if
350   continue
c.... calculate eigenvalues
cww   call elemev(s,nel,ndf,nst)
c
      return
c.... compute and output the element variables
4     continue
c.... calculate nodal values of moments for shear forces
      l=2
      call pgauss(l,lint,sg,tg,wg)
      call jacqudpl(xl,ndm)
      do l = 1,lint
        call strepl(d,xl,ul,shp,xsj,ndm,ndf,sg(l),tg(l),sig)
        do i=1,3
          sigq(l,i)=sig(i)
        end do
      end do
      call nodval(sigq,signo,3)
c      write(*,*) 'GP-values'
c      do i = 1,4 
c                     write(iow,2007) (sigq(i,k),k=1,1)
c        if(ior.lt.0) write(*  ,2007) (sigq(i,k),k=1,1)
c      end do
c      write(*,*) 'nodal-values'
c
c      do i = 1,4 
c                     write(iow,2007) (signo(i,k),k=1,1)
c        if(ior.lt.0) write(*  ,2007) (signo(i,k),k=1,1)
c      end do
c2007  format(6x,3(1x,e13.5))

c.... calculate moments and shear forces at center
      l = 1
      call pgauss(l,lint,sg,tg,wg)
      call jacqudpl(xl,ndm)
c.... moments      
      call strepl(d,xl,ul,shp,xsj,ndm,ndf,sg,tg,sig)
      xx = 0.25*(xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4))
      yy = 0.25*(xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4))
c.... shear forces from equilibrium
      do l = 1,lint
        call shapef(sg(l),tg(l),xl,shpw,xsjw,ndm,.false.)
        sig(4) = 0.d0
        sig(5) = 0.d0
        do i = 1,4
          sig(4) = sig(4) + shpw(1,i)*signo(i,1)+ shpw(2,i)*signo(i,3)  
          sig(5) = sig(5) + shpw(2,i)*signo(i,2)+ shpw(1,i)*signo(i,3)  
        end do
      end do
c...  pressure sigma = c*u
      press = 0.d0
      cb = d(10)
      if(dabs(cb).gt.0.0d0) then
        press = 0.d0
        do l = 1,lint
          call shapef(sg(l),tg(l),xl,shpw,xsjw,ndm,.false.)
          do i = 1,4
            press =  press+dabs(cb)*shpw(3,i)*ul(ndf*(i-1)+1)
          end do     
        end do
        if(cb.lt.0.d0.and.press.ge.0.d0) press = 0.d0 
      end if
c...  design for concrete
      text = ' '
      if(d(11).gt.0.0d0) call design07(d,sig,as,text)
c
      mct = mct - 1
      if(mct.le.0) then
        mct = 50
                     write(iow,2002) o,head
        if(ior.lt.0) write(*  ,2002) o,head
        if(d(11).gt.0.0d0) then
                       write(iow,2004) 
          if(ior.lt.0) write(*  ,2004)
        end if
      end if
                   write(iow,2003) n,ma,xx,yy,sig,press
      if(ior.lt.0) write(*  ,2003) n,ma,xx,yy,sig,press
      if(d(11).gt.0.0d0) then
                     write(iow,2005) text,as(1),as(2),as(3),as(4)
        if(ior.lt.0) write(*  ,2005) text,as(1),as(2),as(3),as(4)
      end if
c.... save stresses for time history plots  ww 27.12.2006????
      do i = 1,3
        bm(i,1) = sig(i)
      end do

      return
c.... compute consistant and lumped mass matrix, CMAS added WW
5     l = 2
      call pgauss(l,lint,sg,tg,wg)
      do 510 l = 1,lint
c....   compute shape functions
        call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.true.)
        xsj = wg(l)*xsj*d(8)
c....   for each node j compute db = rho*shape*xsj
        i1 = 1
        do 500 i = 1,nel
          pm = shp(3,i)*xsj
          p(i1  ) = p(i1)   + pm
          p(i1+1) = p(i1+1) + pm*d(9)
          p(i1+2) = p(i1+2) + pm*d(9)   
          j1 = i1
          do 501 j = i,nel
            sm = pm*shp(3,j)
            s(i1,j1)     = s(i1,j1)     + sm
            s(i1+1,j1+1) = s(i1+1,j1+1) + sm*d(9)
            s(i1+2,j1+2) = s(i1+2,j1+2) + sm*d(9)
  501     j1 = j1 + ndf
  500   i1 = i1 + ndf
  510 continue
c.... make mass symmetric
      do 530 i = 2,nst
        do 530 j = 1,i-1
          s(i,j) = s(j,i)
530   continue
      return
c.... compute the surface loads: const load q= d(1)! in direc. 3
7     call jacqudpl(xl,ndm)
      d6 = d(1)/4.0
      l = 2
      call pgauss(l,lint,sg,tg,wg)
      do 710 l = 1,lint
        call qushp8pl(sg(l),tg(l),shp,xl,ndm,xsj)
        xsj = xsj*wg(l)
        call dktqbmpl(shp(1,5),shp,ndf)
        d6x= d6  *xsj
c.... compute the element load vector
        p(      1) = d6x*(1.-sg(l))*(1.-tg(l)) + p(      1)
        p(  ndf+1) = d6x*(1.+sg(l))*(1.-tg(l)) + p(  ndf+1)
        p(2*ndf+1) = d6x*(1.+sg(l))*(1.+tg(l)) + p(2*ndf+1)
        p(3*ndf+1) = d6x*(1.-sg(l))*(1.+tg(l)) + p(3*ndf+1)
710   continue
      return
c 
8     istv = 13
      if(iplma(ma).eq.0)  return ! only if MATN
      cb = d(10)
c.... stress contours via stress projection from Gauss-points
c     robust for skew meshes
      l = 2
      call pgauss(l,lint,sg,tg,wg)
      call stcnpl1(ix,d,xl,ul,shp,strea,strea(1+numnp),ndm,ndf,numnp,
     1               sg,tg,cb,as)
c.... stress contours via stress calculation at nodes
c     better for regular meshes
c      call stcnpl2(ix,d,xl,ul,shp,strea,strea(1+numnp),ndm,ndf,numnp,
c     +             cb,as)
      return         
c.... error calculation
9     l = 2
      call pgauss(l,lint,sg,tg,wg)
c.... sum internal energy of stress differences for error analysis
      call ster07(ix,d,xl,ul,tl,shp,strea,strea(1+numnp),
     1     e_ome,ndf,ndm,numnp,numel,sg,tg,wg,lint)
      return         
c.... calculate stresses sig(i) at center of element
c.... for qx,qy not implemented
14    l = 1
      call pgauss(l,lint,sg,tg,wg)
      call jacqudpl(xl,ndm)
      call strepl(d,xl,ul,shp,xsj,ndm,ndf,sg,tg,sig)
      sigp(1)=sig(1)
      sigp(2)=sig(3)
      sigp(3)=sig(2)
      call pstres(sigp,sigp(5),sigp(6),sigp(7))
c...  pressure sigma = c*u
      press = 0.d0
      cb = d(10)
      if(dabs(cb).gt.0.0d0) then
         press = 0.d0
         do l = 1,lint
            call shapef(sg(l),tg(l),xl,shpw,xsjw,ndm,.false.)
            do i = 1,4
              press =  press+dabs(cb)*shpw(3,i)*ul(ndf*(i-1)+1)
            end do     
        end do
        if(cb.lt.0.d0.and.press.ge.0.d0) press = 0.d0 
      sigp(4)=press
c...  concrete values
      if(d(11).gt.0.0d0) call design07(d,sig,as,text)
      end if
c
      if(nfp.gt.11.or.nfp.lt.1) return
      if(nfp.ge.1.and.nfp.le.7)  strp=sigp(nfp)
      if(nfp.ge.8.and.nfp.le.11) strp=as(nfp-7)
c
      if(flfp) then
c....   calculate extreme values
        xmaxf = max(xmaxf,strp)
        xminf = min(xminf,strp)
      else
c....   plot stresses 
c....   color      
        call pppcolf(strp) 
c....   coordinates 3d
        call pzero(yl,12)
        do i=1,ndm
          do j=1,4
            yl(i,j)=xl(i,j) 
          end do
        end do
c.....  transform ccordinates for rot        
        call plxtrn(yl,tra,vr,3,4)
c...... plot element
        call plot9(iel,ix,yl,ndm,nel,1)
      end if
      return
c
c.... element load vector
22    call jacqudpl(xl,ndm)
      l = 2
      call pgauss(l,lint,sg,tg,wg)
      call qload07(d6,temp,aqloa,d,numel,n,mqloa,propq)
      d6 = d6/4.0
      do l = 1,lint
        call qushp8pl(sg(l),tg(l),shp,xl,ndm,xsj)
        xsj = xsj*wg(l)
        call dktqbmpl(shp(1,5),shp,ndf)
c....   compute weighted jacobian and d-matrix constants
        d1 = d(1)*xsj
        d2 = d(2)*xsj
        d3 = d(3)*xsj
        d6x= d6  *xsj
c....   compute the element load vector
        p(      1) = d6x*(1.-sg(l))*(1.-tg(l)) + p(      1)
        p(  ndf+1) = d6x*(1.+sg(l))*(1.-tg(l)) + p(  ndf+1)
        p(2*ndf+1) = d6x*(1.+sg(l))*(1.+tg(l)) + p(2*ndf+1)
        p(3*ndf+1) = d6x*(1.-sg(l))*(1.+tg(l)) + p(3*ndf+1)
c....   compute element load vector due to temperature
        do it = 1,12
          ptemp(it) = (bm(1,it)+bm(2,it)) * temp * xsj
        end do
c....   copy into load vector
        kt = 0
        do it = 1,4
          jt=it-1
          p(jt*ndf+1) = ptemp(kt+1) + p(jt*ndf+1)
          p(jt*ndf+2) = ptemp(kt+2) + p(jt*ndf+2)
          p(jt*ndf+3) = ptemp(kt+3) + p(jt*ndf+3)
          kt = kt+3
        end do
      end do
      return
c.... formats
2000  format(5x,'DKQ Plate Bending Element: Material constants'//
     + 10x,'E - modulus ............................. =',e15.5/
     + 10x,'nu  poisson ratio ....................... =',e15.5/
     + 10x,'h   thickness ........................... =',e15.5/
     + 10x,'q   loading ............................. =',e15.5/
     + 10x,'rho density ............................. =',e15.5/
     + 10x,'c   el.found. ........................... =',e15.5,a)
2001  format(10x,'Values for Concrete, Design due to DIN 1045-1 '/
     + 10x,'f_cd (C30/37: 17.0 10^3 kN/m**2)......... =',e15.5/
     + 10x,'f_yd (Bst 500/550 S/M 435 10^3 kN/m**2).. =',e15.5/
     + 10x,'Position Steel d_x ...................... =',e15.5/
     + 10x,'Position Steel d_y ...................... =',e15.5)
2002  format(a1,20a4/10x,
     1'DKQ Plate Moments, Shear forces and Foundation Pressure'/
     2 ' el',1x,'mat',2x,'1-coord',3x,'2-coord',
     3 4x,'moment_11',4x,'moment_22',4x,'moment_12',
     4 4x,'force q_1',4x,'force q_2',4x,'pressure ')
2003  format(i3,i3,2f10.3,6(1x,e12.5))
2004  format(26x,
     1 4x,'as_x(bot)',4x,'as_x(top)',4x,'as_y(bot)',4x,'as_y(top)')
2005  format(6x,a15,5x,4(1x,e12.5))
2006  format(
     1 10x,'temp. coefficient         =',e15.5/
     2 10x,'temp. difference(top-bot) =',e15.5)
5000  format(' Input: E, nu, thick, uniform load, density, el.found.'/
     +       ' Input: beta_R,beta_S,h_x,h_y'/)
5001  format(' Input: alpha_t,delta_t'/)
      end
c----------------------------------------------------------------------
c
      subroutine jacqudpl(xl,ndm)
      implicit double precision (a-h,o-z)
      dimension xl(ndm,*)
      common/elcom071/b(4),c(4),aa(4),bb(4),cc(4),dd(4),ee(4)
!$OMP THREADPRIVATE (/elcom071/)  
      do 100 i = 1,4
        k = mod(i,4) + 1
        b(i) = xl(2,k) - xl(2,i)
        c(i) = xl(1,i) - xl(1,k)
        sql = b(i)*b(i) + c(i)*c(i)
        aa(i) = 1.5*c(i)/sql
        bb(i) = 0.75*b(i)*c(i)/sql
        cc(i) = (0.25*c(i)*c(i) - 0.5*b(i)*b(i))/sql
        dd(i) =-1.5*b(i)/sql
        ee(i) = (0.25*b(i)*b(i) - 0.5*c(i)*c(i))/sql
100   continue
      return
      end
c----------------------------------------------------------------------
      subroutine dktqbmpl(shm,shn,ndf)
      implicit double precision (a-h,o-z)
      dimension shm(3,4),shn(3,4)
      common/elcom071/b(4),c(4),aa(4),bb(4),cc(4),dd(4),ee(4)
!$OMP THREADPRIVATE (/elcom071/)  
      common/elcom072/bm(3,12)
!$OMP THREADPRIVATE (/elcom072/)  
      i1 = 1
      do 110 i = 1,4
        j = mod(i+2,4) + 1
        i2 = i1 + 1
        i3 = i2 + 1
        bm(1,i1) = aa(i)*shm(1,i) - aa(j)*shm(1,j)
        bm(1,i2) = bb(i)*shm(1,i) + bb(j)*shm(1,j)
        bm(1,i3) = cc(i)*shm(1,i) + cc(j)*shm(1,j) - shn(1,i)
        bm(2,i1) = dd(i)*shm(2,i) - dd(j)*shm(2,j)
        bm(2,i2) =-ee(i)*shm(2,i) - ee(j)*shm(2,j) + shn(2,i)
        bm(2,i3) =-bb(i)*shm(2,i) - bb(j)*shm(2,j)
        bm(3,i1) = aa(i)*shm(2,i) - aa(j)*shm(2,j)
     1           + dd(i)*shm(1,i) - dd(j)*shm(1,j)
        bm(3,i2) =-ee(i)*shm(1,i) - ee(j)*shm(1,j) + shn(1,i) - bm(2,i3)
        bm(3,i3) = cc(i)*shm(2,i) + cc(j)*shm(2,j) - shn(2,i) - bm(1,i2)
        i1 = i1 + 3
110   continue
      return
      end
c----------------------------------------------------------------------
      subroutine qushp8pl(s,t,shp,xl,ndm,xsj)
c.... shape function routine for 8-node serendipity elements
      implicit double precision (a-h,o-z)
      dimension   xl(ndm,*),si(4),ti(4),shp(3,8)
      data si/-1.,1.,1.,-1./,ti/-1.,-1.,1.,1./
      xs = 0.0
      xt = 0.0
      ys = 0.0
      yt = 0.0
      do 100 i = 1,4
        ss = si(i)*s
        tt = ti(i)*t
        sn = si(i)*(1.+ tt)
        tn = ti(i)*(1.+ ss)
        xs = xs + sn*xl(1,i)
        xt = xt + tn*xl(1,i)
        ys = ys + sn*xl(2,i)
        yt = yt + tn*xl(2,i)
        shp(1,i) = 0.25*sn*(ss+ss+tt)
        shp(2,i) = 0.25*tn*(ss+tt+tt)
        shp(3,i) = 0.25*(1.+ss)*(1.+tt)*(-1.+ss+tt)
100   continue
      xsj = (xs*yt-ys*xt)/4.0
      xs = xs/xsj
      xt = xt/xsj
      ys = ys/xsj
      yt = yt/xsj
      xsj = xsj/4.0
      do 110 i = 5,7,2
        ss = si(i-4)*s
        tt = ti(i-4)*t
        shp(1,i) = -s*(1.+tt)
        shp(2,i) = 0.5*ti(i-4)*(1.-s*s)
        shp(3,i) = 0.5*(1.-s*s)*(1.+tt)
        shp(1,i+1) = -0.5*si(i-4)*(1.-t*t)
        shp(2,i+1) = -t*(1.-ss)
        shp(3,i+1) = 0.5*(1.-ss)*(1.-t*t)
110   continue
      do 120 i = 1,8
        sn       = yt*shp(1,i) - ys*shp(2,i)
        shp(2,i) = xs*shp(2,i) - xt*shp(1,i)
        shp(1,i) = sn
120   continue
      return
      end
c
      subroutine strepl(d,xl,ul,shp,xsj,ndm,ndf,sg,tg,sig)
c----------------------------------------------------------------------
c     purpose: Calculate straains and stresses at integration point     
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xl(ndm,*),d(*),ul(*),eps(3),shp(3,8),sig(3),sg(4),tg(4)
      common /elcom072/ bm(3,12)
!$OMP THREADPRIVATE (/elcom072/)  
c.... strains
      call qushp8pl(sg,tg,shp,xl,ndm,xsj)
      call dktqbmpl(shp(1,5),shp,ndf)
      do 420 i = 1,3
        eps(i) = 0.0
cww     do 410 j = 1,12
cww       eps(i) = eps(i) + bm(i,j)*ul(j)
cww410  continue
          do jj = 1,4
          do kk = 1,3
            i1 = (jj-1)*3  + kk
            j1 = (jj-1)*ndf+ kk
              eps(i)=eps(i)+bm(i,i1)*ul(j1)
          end do
        end do
420   continue
c.... strains due to temperature loads
      eps(1) = eps(1) - d(16) 
      eps(2) = eps(2) - d(16) 
c.... stresses
      sig(1) = d(1)*eps(1) + d(2)*eps(2)
      sig(2) = d(2)*eps(1) + d(1)*eps(2)
      sig(3) = d(3)*eps(3)
      return
      end
c
      subroutine stcnpl1(ix,d,xl,ul,shp,dt,st,ndm,ndf,numnp,sg,tg,cb,as)
c----------------------------------------------------------------------
c.... Purpose: stress projection and energy                             
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ix(4)
      dimension st(numnp,*),dt(numnp),xl(ndm,*),ss(4),tt(4),as(4),
     1          d(*),ul(*),shp(3,8),sg(4),tg(4),sig(5),cmat(3,3),csig(3)
     2         ,sigq(4,3),signo(4,3),shpw(3,4) 
      character*15 text
      data ss/-1.0,1.0,1.0,-1.0/,tt/-1.0,-1.0,1.0,1.0/
      call pzero(as,4)
c.... set elasticity matrix
      call pzero(cmat,9)
      cmat(1,1) = d(1)
      cmat(1,2) = d(2)
      cmat(2,2) = d(1)
      cmat(3,3) = d(3)
      call msym (cmat,3,3,2)         ! C    = C^T
      call invert(cmat,3,3)          ! C    = C^(-1)

      call jacqudpl(xl,ndm)
c.... calculate nodal values of moments for shear forces
      do l = 1,4
        call strepl(d,xl,ul,shp,xsj,ndm,ndf,sg(l),tg(l),sig)
        do i=1,3
          sigq(l,i)=sig(i)
        end do
      end do
      call nodval(sigq,signo,3)
c
c.... calculate stresses
      do 410 l = 1,4
c....   moments
        call strepl(d,xl,ul,shp,xsj,ndm,ndf,sg(l),tg(l),sig)
c....   shear forces
        call shapef(sg(l),tg(l),xl,shpw,xsjw,ndm,.false.)
        sig(4) = 0.d0
        sig(5) = 0.d0
        do i = 1,4
          sig(4) = sig(4) + shpw(1,i)*signo(i,1)+ shpw(2,i)*signo(i,3)  
          sig(5) = sig(5) + shpw(2,i)*signo(i,2)+ shpw(1,i)*signo(i,3)  
        end do
c
c....   compute element energy
        call enerel07(sig,cmat,csig,xsj,3)
c
c...    design for concrete
        if(d(11).gt.0.0d0) then
          call design07(d,sig,as,text)
        end if
          do 400 ll = 1,4
            ii = abs(ix(ll))
            if(ii.eq.0) go to 400
            xsji = xsj*(1.+ss(ll)*sg(l))*(1.+tt(ll)*tg(l))
            sz = dabs(sig(1)+sig(2)+sig(3))
            if(sz.lt.1e-5) goto 400
            dt(ii) = dt(ii) + xsji
            st(ii, 1) = st(ii,1) + sig(1)*xsji
            st(ii, 2) = st(ii,2) + sig(3)*xsji
            st(ii, 3) = st(ii,3) + sig(2)*xsji
            press = dabs(cb)*ul(ndf*(l-1)+1) 
            if(cb.lt.0.d0.and.press.gt.0.d0) press = 0.d0
            st(ii, 4) = st(ii, 4) + press *xsji
            st(ii, 8) = st(ii, 8) + sig(4)*xsji
            st(ii, 9) = st(ii, 9) + sig(5)*xsji
            st(ii,10) = st(ii,10) + as (1)*xsji
            st(ii,11) = st(ii,11) + as (2)*xsji
            st(ii,12) = st(ii,12) + as (3)*xsji
            st(ii,13) = st(ii,13) + as (4)*xsji
400       continue
410   continue
      return
      end
c
      subroutine stcnpl2(ix,d,xl,ul,shp,dt,st,ndm,ndf,numnp,cb,as)
c----------------------------------------------------------------------
c.... stress projection and energy                                    |
c.... for qx,qy not implemented
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ix(4)
      dimension st(numnp,*),dt(numnp),xl(ndm,*),ss(4),tt(4),
     1          d(*),ul(*),shp(3,8),sig(3),as(4),cmat(3,3),csig(3)
      character*15 text
      data ss/-1.0,1.0,1.0,-1.0/,tt/-1.0,-1.0,1.0,1.0/
      call pzero(as,4)
c.... set elasticity matrix
      call pzero(cmat,9)
      cmat(1,1) = d(1)
      cmat(1,2) = d(2)
      cmat(2,2) = d(1)
      cmat(3,3) = d(3)
      call msym (cmat,3,3,2)         ! C    = C^T
      call invert(cmat,3,3)          ! C    = C^(-1)

      call jacqudpl(xl,ndm)
      do 400 l = 1,4
          call strepl(d,xl,ul,shp,xsj,ndm,ndf,ss(l),tt(l),sig)
c....   compute element energy
        call enerel07(sig,cmat,csig,xsj,3)
c...    design for concrete
        if(d(11).gt.0.0d0) then
          call design07(d,sig,as,text)
        end if

          ii = ix(l)
          if(ii.le.0) go to 400
          xsji = xsj*(1.+ss(l)*ss(l))*(1.+tt(l)*tt(l))
          dt(ii)   = dt(ii) + xsji
          st(ii, 1) = st(ii, 1) + sig(1)*xsji
          st(ii, 2) = st(ii, 2) + sig(3)*xsji
          st(ii, 3) = st(ii, 3) + sig(2)*xsji
          press = dabs(cb)*ul(ndf*(l-1)+1) 
          if(cb.lt.0.d0.and.press.gt.0.d0) press = 0.d0
          st(ii, 4) = st(ii, 4) + press *xsji
          st(ii, 8) = st(ii, 8) + as (1)*xsji
          st(ii, 9) = st(ii, 9) + as (2)*xsji
          st(ii,10) = st(ii,10) + as (3)*xsji
          st(ii,11) = st(ii,11) + as (4)*xsji
400     continue
      return
      end
c
      subroutine pressure(press,cb,ul,ndf,nel)
c--------------------------------------------------------------------+
c     calculate pressure due to elastic foundation: sigma = c*u      |
c--------------------------------------------------------------------+
      implicit double precision (a-h,o-z)
      dimension press(4),ul(ndf,*)
      do i = 1,nel
        press(i) = dabs(cb)*ul(1,i)
      end do
      return
      end
c
      subroutine design07(d,sig,as,text)
c--------------------------------------------------------------------+
c     concrete design  (from general design diagram) DIN 1045-1      |
c     mue_Ed = m / (b*d**2*f_cd)                                     |
c     b      = 1                                                     |
c     f_cd   = alpha *f_ck/gamma_c                                   |
c     zeta   = 0.5 + sqrt(0.25 - 0.5*mue_Ed)       (Approximation)   |
c     z      = zeta*d                                                |  
c     ohne Druckbewehrung                                            |
c     as1    = m_Ed /(z*sig_s1d)                                     |
c                                                                    |
c     mit Druckbewehrung                                             |
c     m_Ed_lim   = mue_Ed_lim*(b*d**2*f_cd)                          |
c     Delta m_Ed = m_ed - m_Ed_lim                                   |
c     as1    = m_Ed_lim /(z*sig_s1d)+[Delta m_Ed /((d-d_2)*sig_s1d)] |
c     as2    =                       [Delta m_Ed /((d-d_2)*sig_s2d)] |
c                                                                    |
c     f_yd   = f_yk/gamma_s                                          |
c     Sigma_sd = f_yd                    für mue_Ed < mue_Ed_lim     |
c     Sigma_sd = f_yd(4.906-10.5 mue_Ed) für mue_Ed > mue_Ed_lim     |
c     mue_Ed_lim = 0.372                                             | 
c                                                                    |
c     load and material are given with respect to ultimate load
c     e.g. gamma_Q =1.35
c          gamma_s =1.15 -> f_yd
c          gamma_c =1.5  -> f_cd
c--------------------------------------------------------------------+
      implicit double precision (a-h,o-z)
      dimension d(*),sig(3),as(4)
      character*15 text
      gmued  = 0.372d0
      f_cd   = d(11)
      f_yd   = d(12)
      dx     = d(13)
      dy     = d(14)
      h      = d(15)
      asxadd = 0.d0
      asyadd = 0.d0
      call pzero(as,4)
      text = ' '
c.... design
c.....m_x, m_xy added not correct, see EC2  
      if(sig(1) .ge. 0.0d0) smx=      sig(1)  + dabs(sig(3))  
      if(sig(1) .lt. 0.0d0) smx= dabs(sig(1)) + dabs(sig(3)) 
      xmued = smx /(dx*dx*f_cd)
      val   = 0.25d0 - 0.5d0*xmued
      if(val.lt.0.d0) then ! provisorial 
        write(*,*) 'zeta  not available '  
        val=0.d0
      end if
      zeta  = 0.5d0  + dsqrt(val)
      sigma_sd = f_yd 
      if(xmued .gt. gmued) then
        sigma_sd = f_yd*(4.906-10.5*xmued)
        text = 'not economic!  ' 
        if(sigma_sd.lt.0.d0) then ! provisorial
          write(*,*)  'Sigma_sd < 0, load too high '  
          sigma_sd=1.e-3
        end if  
        dsmx  =  (xmued - gmued)*dx*dx*f_cd
        smx   =  smx  - dsmx
c....   approx. distance for steel d-d_2  
        dmd2x  = dx - (h-dx)
        asxadd = dsmx/(sigma_sd*dmd2x)
      end if
      asx = smx / (zeta*dx*sigma_sd)
      if(sig(1) .ge. 0.0d0) then
        as(1) = asx + asxadd
        as(2) =       asxadd
      elseif(sig(1) .lt. 0.0d0) then
        as(1) =       asxadd
        as(2) = asx + asxadd
      end if
c.....m_y m_xy added not correct, see EC2
      if(sig(2) .ge. 0.0d0) smy=      sig(2)  + dabs(sig(3))
      if(sig(2) .lt. 0.0d0) smy= dabs(sig(2)) + dabs(sig(3))
      ymued = smy /(dy*dy*f_cd)
      val   = 0.25d0 - 0.5d0*ymued
      if(val.lt.0.d0) then ! provisorial 
        write(*,*) 'zeta  not available '  
        val=0.d0
      end if
      zeta  = 0.5d0  + dsqrt(val)
      sigma_sd = f_yd 
      if(ymued .gt. gmued) then
        sigma_sd = f_yd*(4.906-10.5*ymued)
        text = 'not economic!  ' 
        if(sigma_sd.lt.0.d0) then ! provisorial
          write(*,*)  'Sigma_sd < 0, load too high '  
          sigma_sd=1.e-3
        end if  
        dsmy  =  (ymued - gmued)*dy*dy*f_cd
        smy   =  smy  - dsmy
c....   approx. distance for steel d-d_2  
        dmd2y  = dy - (h-dy)
        asyadd = dsmy/(sigma_sd*dmd2y)
      end if
      asy = smy / (zeta*dy*sigma_sd)
      if(sig(2) .ge. 0.0d0) then
        as(3) = asy + asyadd
        as(4) =       asyadd
      elseif(sig(2) .lt. 0.0d0) then
        as(3) =       asyadd
        as(4) = asy + asyadd
      end if
      return
      end
c
      subroutine designalt(d,sig,as)
c--------------------------------------------------------------------+
c     concrete design  (from general design diagram) DIN 1045 alt    |
c     ms    = m / (b*h**2*beta_R)                                    |
c     k_z   = 0.5 + sqrt(0.25 - 0.5*gamma*ms)     (Approximation)    |
c     asz   = gamma * mb /(k_z*h*beta_S)+[gamma * me /((h*)*beta_S)] |
c     asd   =                            [gamma * me /((h*)*beta_S)] |
c     b     = 1                                                      |
c     gamma = 1.75                                                   |
c--------------------------------------------------------------------+
      implicit double precision (a-h,o-z)
      dimension d(*),sig(3),as(4)
      smg    = 0.337d0
      betar  = d(11)
      betas  = d(12)
      thickx = d(13)
      thicky = d(14)
      thick  = d(15)
      gamma  = 1.75d0
      dsmxg  = 0.d0
      dsmyg  = 0.d0
      asxz   = 0.d0
      asyz   = 0.d0
      asxd   = 0.d0
      asyd   = 0.d0
      call pzero(as,4)
c.... design
c.....m_x
      if(sig(1) .ge. 0.0d0) smx=      sig(1)  * gamma
      if(sig(1) .lt. 0.0d0) smx= dabs(sig(1)) * gamma
      smxg  = smx /(thickx*thickx*betar)
      if(smxg .gt. smg) then
        dsmxg = (smxg - smg)
        dsmx  =  dsmxg*thickx*thickx*betar
        smxg  =  smg
        smx   =  smx  - dsmx
c....   approx. distance for steel h-h'  
        dhx   = thickx - (thick-thickx)
c....   tension
        asxz  = asxz + dsmx/(betas*dhx)
c....   compression
        asxd  = asxd + dsmx/(betas*dhx)
      end if
      val   = 0.25d0 - 0.5d0*smxg
      zk    = 0.5d0 + dsqrt(val)
      asxz  = asxz  + smx / (zk*thickx*betas)
      if(sig(1) .ge. 0.0d0) then
        as(1) = as(1) + asxz
        as(2) = as(2) + asxd
      elseif(sig(1) .lt. 0.0d0) then
        as(2) = as(2) + asxz
        as(1) = as(1) + asxd
      end if
c.....m_y
      if(sig(2) .ge. 0.0d0) smy=      sig(2)  * gamma
      if(sig(2) .lt. 0.0d0) smy= dabs(sig(2)) * gamma
      smyg  = smy /(thicky*thicky*betar)
      if(smyg .gt. smg) then
        dsmyg = (smyg - smg)
        dsmy  =  dsmyg*thicky*thicky*betar
        smyg  =  smg
        smy   =  smy  - dsmy
c....   approx. distance for steel h-h'  
        dhy   = thicky - (thick-thicky)
c....   tension
        asyz  = asyz + dsmy/(betas*dhy)
c....   compression
        asyd  = asyd + dsmy/(betas*dhy)
      end if
      val   = 0.25d0 - 0.5d0*smyg
      zk    = 0.5d0 + dsqrt(val)
      asyz  = asyz  + smy / (zk*thicky*betas)
      if(sig(2) .ge. 0.0d0) then
        as(3) = as(3) + asyz
        as(4) = as(4) + asyd
      elseif(sig(2) .lt. 0.0d0) then
        as(4) = as(4) + asyz
        as(3) = as(3) + asyd
      end if
      return
      end
c
      subroutine ster07(ix,d,xl,ul,tl,shp,dt,st,e_ome,
     1      ndf,ndm,numnp,numel,sg,tg,wg,lint)
c ----------------------------------------------------------------------
c.... error calculation                                                |
c ----------------------------------------------------------------------
      USE eldata
      USE errin1
      USE iofile
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),xl(ndm,*),tl(*),ix(*),shp(3,4),
     1   sig(3),sigp(4,3), d(*),ul(ndf,*),sg(*),tg(*),wg(*),sigg(4,3)
      dimension xsje(4),cmat(3,3),e_ome(*)
      call pzero(sigp,4*3)
c.... elasticity matrix (sym!)
      call pzero(cmat,9)
      cmat(1,1) = d(1)
      cmat(1,2) = d(2)
      cmat(2,2) = d(1)
      cmat(3,3) = d(3)
      call jacqudpl(xl,ndm)
      do ii = 1,lint
c....   stresses from material law
        call strepl(d,xl,ul,shp,xsj,ndm,ndf,sg(ii),tg(ii),sig)
        do i = 1,3 
          sigg(ii,i) = sig(i)
        end do
c
c....   stresses from nodal points
        call shape(sg(ii),tg(ii),xl,shp,xsj,ndm,nel,ix,.false.)
        xsje(ii) = xsj*wg(ii)
        do i = 1,nel
            ll = iabs(ix(i))
          if(ll.ne.0) then
            sigp(ii,1) = sigp(ii,1) + shp(3,i)*st(ll,1)
            sigp(ii,2) = sigp(ii,2) + shp(3,i)*st(ll,3)  ! s_x,s_y,s_xy
            sigp(ii,3) = sigp(ii,3) + shp(3,i)*st(ll,2)
            end if
        end do
      end do      
      numsig = 3
      call erro07(ix,xl,ndm,numel,lint,xsje,numsig,
     1            sigg,sigp,cmat,e_ome)    
      end
c 
      subroutine erro07(ix,xl,ndm,numel,lint,xsje,numsig,
     1                  sigg,sigp,cmat,e_ome)
c ----------------------------------------------------------------------
c.... calculate element stress differences
c ----------------------------------------------------------------------
      USE errin1
      USE errin2
      implicit double precision (a-h,o-z)
      dimension sigg(lint,numsig),sigp(lint,numsig),xsje(lint),
     1     dsig(3),csig(3),cmat(numsig,*),e_ome07(numerr),e_ome(numel,*)

      call msym  (cmat,numsig,numsig,2)        ! C   = C^T
      call invert(cmat,numsig,numsig)          ! C   = C^(-1)
c.... intial values for element errors
      e_ome07 = 0.0d0
      do i = 1,lint
c.... stress differences
        do j = 1,numsig
          dsig(j) = sigp(i,j) - sigg(i,j)
        end do
        call mvmul (cmat,dsig,numsig,numsig,csig) ! csig = C^(-1)*dsig
c....   element errors
        do j = 1,numsig
          e_ome07(1) = e_ome07(1) +  dsig(j)*csig(j)*xsje(i)
          e_ome07(2) = e_ome07(2) + (dsig(j)**2)    *xsje(i)
        end do
      end do
c.....plot/print/add errors
      call elmterr(ix,xl,ndm,numel,e_ome07,e_ome)
      return
      end
c
      subroutine enerel07(sig,cmat,csig,xsj,numsig)
c ----------------------------------------------------------------------
c.... compute element energy
c ----------------------------------------------------------------------
      USE errin1
      implicit double precision (a-h,o-z)
      dimension sig(*),cmat(numsig,*),csig(*)
      call mvmul (cmat,sig,numsig,numsig,csig) ! csig = C^(-1)*sig
      do i = 1,numsig
        u_om(1) = u_om(1) + sig(i)*csig(i)*xsj 
        u_om(2) = u_om(2) + sig(i)**2     *xsj 
      end do
      end
c      
      subroutine qload07(d6,temp,q,d,numel,n,mqloa,propq)
c----------------------------------------------------------
c.... calculate loads from macro qloa
c
c     q(01) = q 
c     q(02) = dT
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*)
      d6   = 0.d0 
      temp = 0.d0
      if(mqloa.ne.1) then
        d6   =             q(n,1)*propq 
        temp = (d(1)+d(2))*q(n,2)*propq 
      end if 
      return
      end
c      
      subroutine nodval(xgp,xno,n)
c-----------------------------------------------------------------------
c     purpose: calculate nodal values of a stress resultant from
c              values at Gauss-points
c
c     input:   xgp(4,n) - values at gauss-points
c              n        - number of values 
c
c     output:  xno(4,n) - values at nodes
c
c     valid for: 4 node, 4 Gauss-points
c
c     theory: m_GP(xi_p,eta_q) = sum N_k(xi_p,eta_q)*m_k at 4 Gauss points
c             Aik*m_k = m_GP 
c             m_k = Aik^-1*m_GP
c
c     see also 
c     Hinton E, Campbell JS. Local and global smoothing of discontinuous 
c     finite element functions  using a least squares method. 
c     IJNME 1974; 8:461– 480.
c     and
c     Kapuria S., and Kulkarni, S. D. 
c     An improved discrete Kirchhoff quadrilateral element based
c     on third-order zigzag theory for static analysis of composite
c     and sandwich plates, IJNME 2007; 69:1948–1981(eq. 43) 
c
c     ww ibs uka 12/06
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ai(4,4),xgp(4,3),xno(4,3) 

c.....coefficient matrix 
      cw = dsqrt(3.d0)*0.5d0 
      c1 =  1.d0 + cw
      c2 = -0.5d0
      c3 =  1.d0 - cw  
  
      ai(1,1)=c1
      ai(2,2)=c1
      ai(3,3)=c1
      ai(4,4)=c1

      ai(1,2)=c2  
      ai(2,3)=c2  
      ai(3,4)=c2  
      ai(1,4)=c2  

      ai(2,1)=c2  
      ai(3,2)=c2  
      ai(4,3)=c2  
      ai(4,1)=c2  

      ai(1,3)=c3  
      ai(2,4)=c3  
      ai(3,1)=c3  
      ai(4,2)=c3  
 
c.....calculate nodal values
      do j = 1,n 
        do i = 1,4
          xno(i,j) = 0.d0
          do k = 1,4 
            xno(i,j) = xno(i,j)+ ai(i,k)*xgp(k,j)
          end do
        end do
      end do
      return
      end 