      subroutine mate3d15(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c
c     Purpose: calculate S and C for 
c               a linear elastic transversal isotropic material law 
c
c     Inputs:
c         h1(nh) - history array h1
c         d(md)  - local d-array
c         Eps    - strains  
c         isw    - solution option from element 
c
c     Input material parameter:                
c         E_11   - Young's modulus
c         E_22   - Young's modulus
c         v_12   - Poisson's ratio
c         G_12   - Shear   modulus  G_12 = G_13
c         G_23   - Shear   modulus
c         phi    - Number of layers
c
c     Outputs:
c         md = 30 - number of used data for control of d-array //mdx in main subroutine
c         nh = 0 - number of history parameter at Gauss-Point
c         sig    - stresses
c         cmat   - tangent modulus
c
c     Allocation of d-array: (from d(21) from elmt45 subroutine)
c         d(1) = E_11 - Young's modulus
c         d(2) = E_22 - Young's modulus
c         d(3) = v_12 - Poisson's ratio
c         d(4) = G_12 - Shear   modulus  G_12 = G_13
c         d(5) = G_23 - Shear   modulus
c
c
c     Input data 4. line - data for failure models d(7)-d(2)      
c     read failure values        d(7)-d(13); Cuntze: d(7)-d(14)
c     set          values        d(15)-d(22)
c     d(7)= iftyp (failure criterion)
c            0 = no damage
c            1 = Tsai-Wu criteria
c            2 = Goyal criteria
c            3 = Puck criteria
c            4 = Cuntze criteria 
c
c-----------------------------------------------------------------------
c
c     Input data 5. line - data for degradation models d(23)-d(29)      
c     read degradation values    d(23)-d(27)
c         set          values    d(28)-d(29)
c     d(23)= idtyp (degradation type: 0 or damage criterion)
c            0 = no degradation
c            1 = Tsai-Wu degradation
c            2 = Goyal degradation
c            3 = Chang&Lessard degradation
c            4 = Puck degradation
c            5 = linear degradation
c-----------------------------------------------------------------------
c
c     Comments WW
c     getestet nur mit ELMT45 und einem 1 Element in Dickenrichtung 
c     der Schale sowie 1-n Schichten im Element 
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension eps(6),sig(6),cmat(6,6),d(*),h1(*),h2(*),hdam(6)
     '          ,plout(10)
c
      if(isw.eq.1) then
c....   input data , write d
        call mati3d15(d,md,nh)
c
      elseif (isw.eq.3.or.isw.eq.4.or.isw.eq.8) then
c
        call mod3d15(cmat,cappa,d,eps,idam,fi,isw,h1,h2,hdam,plout)
c
c....   stresses sig(i)=cmat(i,j)*eps(j) (mvmul = matrix*vector multiplication)
        call mvmul(cmat,eps,6,6,sig)
c      
      end if 
c        
      return
      end
c
      subroutine mati3d15(d,md,nh)
c-----------------------------------------------------------------------
c
c     mate3d15 Linear elastic transversal isotropic with damage
c     input material parameter 
c     and some additional expanded input stuff by M. Krawiec
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*)
c
      if(ior.lt.0) write(*,1001)
1001  format(' Input: E_1,E_2,nu_12,G_12,G_23,')
c
      call dinput(d,5)
c
      v23   = 0.5d0*d(2)/d(5)-1.d0     
c       
                  write(iow,1002) (d(i),i=1,5),v23
      if(ior.lt.0)write(*  ,1002) (d(i),i=1,5),v23
c
c     Input failure parameters
      call paraset15(d,1)
c
c     Input degradation parameters
      call paraset15(d,2)
c
c     update md for ply values
      md    = 30
      d(30) = md
c
      iftyp = d (7)
      idtyp = d (23)   
c
c.... set values for h-array & check
      nh = 0
c
      if(iftyp.eq.0)      then  ! dummy
        nh=0  
      else if(iftyp.eq.1) then
        if(idtyp.eq.1) nh=1   ! AERNNOVA criteria + degradation
      else if(iftyp.eq.2) then
c... Extended HASHIN criteria + degradation AERNNOVA
        if(idtyp.eq.1) nh=3   
c... Extended HASHIN criteria + degradation edited CRC-ACS  model
        if(idtyp.eq.2) nh=3  
      else if(iftyp.eq.3) then              ! Puck criteria
        if(idtyp.eq.1.or.idtyp.eq.2)then
          write(*,1009) ; stop
        endif
        if(idtyp.eq.0) nh=6   ! no damage
        if(idtyp.eq.3) nh=6   ! CHANG&LESSARD degradation
        if(idtyp.eq.4) nh=9   ! PUCK degradation
      else if(iftyp.eq.4) then
        if(idtyp.eq.1.or.idtyp.eq.2)then
          write(*,1009) ; stop
        endif
        if(idtyp.eq.0) nh=5    ! no damage
        if(idtyp.eq.3) nh=5    ! CHANG&LESSARD degradation
        if(idtyp.eq.5) nh=8   ! linear degradation
      else
        write(*,1008) ; stop
      end if

                  write(iow,1003) nh
      if(ior.lt.0)write(*  ,1003) nh

c      
c....formats
1002  format(5x,'Linear elastic transversal isotropic material data',/,
     +  5x,'elastic modulus E_1 .......',g12.5,/,
     +  5x,'elastic modulus E_2 .......',g12.5,/,
     +  5x,'poissons ratio v_12 .......',g12.5,/,
     +  5x,'shear modulus  G_12 .......',g12.5,/,
     +  5x,'shear modulus  G_23 .......',g12.5,/,
     +  5x,'poissons ratio v_23 .......',g12.5,/)
1003  format(5x,'length nh of h1,h2 ........',i12)
1008  format(5x,'Invalid failure criteria !')
1009  format(5x,'Invalid degradation and failure model combination !')
c
      return
      end 
c
c
      subroutine paraset15(d,isdam)
c----------------------------------------------------------------------
c
c       Purpose:
c       read in the d field
c
c     isdam =1:  read failure     values    d(7)-d(14)
c                         iftyp 4 cuntze    d(7)-d(15)         
c                set              values    d(15)-d(22)
c                         iftyp 4 cuntze    d(16)-d(17)
c                    
c     isdam =2:  read degradation values    d(23)-d(28)
c                set              values    d(29)-d(30)
c
c            
c              Former Values of elmt45
c
c     isw =1:  read failure     values    d(21)-d(28)
c              set              values    d(29)-d(36)
c     isw =2:  read degradation values    d(37)-d(44)
c              set              values    d(45)-d(46)
c
c
c....     data description
c  7      iftyp 1 or 2
c  8      X_t   Bruchspannung bei Zug in Faserrichtung          sig+11
c  9      X_c   Bruchspannung bei Druck in Faserrichtung        sig-11 
c  10     Y_t   Bruchspannung bei Zug quer zur Faserrichtung    sig+22 
c  11     Y_c   Bruchspannung bei Druck quer zur Faserrichtung  sig-22 
c  12     S_12  Bruchschubspannung 12                           tau 12               
c  13     S_13  Bruchschubspannung 13                           tau 13    
c  14     S_23  Bruchschubspannung 23                           tau 23
c
c  7      iftyp 3
c  8      Rt11  Zugfestigkeit in Faserrichtung                  sig+11
c  9      Rc11  Druckfestigkeit in Faserrichtung                sig-11 
c  10     Rt22  Zugfestigkeit quer zur Faserrichtung            sig+22 
c  11     Rc22  Druckfestigkeit quer zur Faserrichtung          sig-22 
c  12     R21   Festigkeit quer / längs                         tau 12
c  13     pt21  Neigungsparameter                       
c  14     pc21  Neigungsparameter                      
c  15     Ra22  Widerstand der Bruchfläche            
c  16     p22   Neigungsparameter                          
c  17     tauc  max. Schubspannung                         
c
c  7      iftyp 4
c  8      Rt11  Zugfestigkeit in Faserrichtung
c  9      Rc11  Druckfestigkeit in Faserrichtung
c  10     Rt22  Zugfestigkeit quer zur Faserrichtung
c  11     Rc22  Druckfestigkeit quer zur Faserrichtung
c  12     R21   Schubfestigkeit     
c  13     b21   Parameter IFF2            sinnvoll: 0.05<b21<0.15 // 0.1
c  14     b22   Parameter IFF3            sinnvoll: 1.00<b22<1.15 // 1.1
c  15     rop   Interaktionsparameter     sinnvoll: rop im Bereich von 3.0 
c  16     tmin  min. tau für IFF2 (Zugbereich)
c  17     tmax  max. tau für IFF2 (Druckbereich)
c
c----------------------------------------------------------------------
c
      USE iofile
        implicit double precision(a-h,o-z)   
        dimension d(*)
c
       mdx=d(30)
      if(isdam.eq.1) then     
c
        call pzero(d(7),16)
        call dinput(d(7),9)
c
c....   read failure data       
        iftyp=d(7)

c....  check for iftyp and set the failure values
      if(iftyp.eq.1) then
c......   calculate parameters for Tsai-Wu criterion
          fxy  =0.5d0
          d(15)=1.d0/d(8)-1.d0/d(9)   ! q1 = 1/X_t - 1/X_c
          d(16)=1.d0/(d(8)*d(9))      ! a1 = 1/(X_t * X_c)
          d(17)=1.d0/d(10)-1.d0/d(11) ! q2 = 1/Y_t - 1/Y_c
          d(18)=1.d0/(d(10)*d(11))    ! a2 = 1/(Y_t * Y_c)
          d(19)=2.d0*fxy/dsqrt(d(8)*d(9)*d(10)*d(11))
          d(20)=1.d0/(d(12)*d(12))    ! a4 = 1/(S_12 * S_12)
          d(21)=1.d0/(d(13)*d(13))    ! a5 = 1/(S_13 * S_13)
          d(22)=1.d0/(d(14)*d(14))    ! a6 = 1/(S_23 * S_23)
          write(iow,2001) (d(k),k=7,22)       
      else if(iftyp.eq.2) then
c......   calculate parameters for extended Hashin criterion
          d(15) = 1.d0/d(8)/d(8)      ! qxt
          d(16) = 1.d0/d(9)/d(9)      ! qxc
          d(17) = 1.d0/d(10)/d(10)    ! qyt
          d(18) = 1.d0/d(11)/d(11)    ! qyc
          d(19) = 1.d0/d(12)/d(12)    ! qs12
          d(20) = 1.d0/d(13)/d(13)    ! qs13
          d(21) = 1.d0/d(14)/d(14)    ! qs23
          write(iow,2002) (d(k),k=7,21)
        else if(iftyp.eq.3) then
c......   calculate parameters for Puck criterion
          d(15) = d(12)/2.d0/d(14)*(dsqrt(1.d0+2.d0*d(14)*d(11)/
     +            d(12))-1.d0)                      ! ra22
          d(16) = d(14)*d(15)/d(12)              ! p22
          d(17) = d(12)*dsqrt(1.d0+2.d0*d(16))      ! tauc
          write(iow,2003) (d(k),k=7,17)
        else if(iftyp.eq.4) then
c......   calculate parameters for Cuntze criterion
          tmin = d(12)
          tmax = d(12)
          do i=1,30
            tmin = (d(12)**3.d0-2*d(13)*d(10)*tmin**2.d0)**(1.d0/3.d0)
            tmax = (d(12)**3.d0+2*d(13)*d(11)*tmax**2.d0)**(1.d0/3.d0)
          enddo
          d(16) = tmin
          d(17) = tmax
          write(iow,2004) (d(k),k=7,17)
        else
          write(iow,600) iftyp
          stop
      endif
c
c.... READ DEGRADATION VALUES
       else if(isdam.eq.2) then
c  
        call pzero(d(23),8)
        call dinput(d(23),6)
        idtyp=d(23)
c
        iftyp=d(7)
        if(iftyp.eq.0) then
          if(idtyp.eq.1) then
            write(iow,3001) (d(k),k=23,23)
          else
            write(iow,601) idtyp
            stop
          endif
        else if(iftyp.eq.1) then
          if(idtyp.eq.1) then 
            write(iow,3001) (d(k),k=23,23)
          else
            write(iow,601) idtyp
            stop
          endif
        else if(iftyp.eq.2) then
          if(idtyp.eq.2) then
            write(iow,3002) (d(k),k=23,23)
          else
            write(iow,601) idtyp
            stop
          endif
        else if(iftyp.eq.3) then
          if(idtyp.eq.3.or.idtyp.eq.4) then 
            write(iow,3003) (d(k),k=23,24)
          else
            write(iow,601) idtyp
            stop
          endif
        else if(iftyp.eq.4) then
          if(idtyp.eq.3.or.idtyp.eq.5) then
            write(iow,3004) (d(k),k=23,24)
          else
            write(iow,601) idtyp
            stop
          endif
        end if
       endif
c
c.... formats
2001  format(
     + 5x,'TSAI-WU failure criterion (values without sign) ',/,
     + 5x,'Versagenskriterium                                ',f12.1,/,
     + 5x,'X_t   Bruchspannung Zug in Faserrichtung          ',e12.5,/,
     + 5x,'X_c   Bruchspannung Druck in Faserrichtung        ',e12.5,/,
     + 5x,'Y_t   Bruchspannung Zug quer zur Faserrichtung    ',e12.5,/,
     + 5x,'Y_c   Bruchspannung Druck quer zur Faserrichtung  ',e12.5,/,
     + 5x,'S_12  Bruchspannung Schub 12                      ',e12.5,/,
     + 5x,'S_13  Bruchspannung Schub 13                      ',e12.5,/,
     + 5x,'S_23  Bruchspannung Schub 23                      ',e12.5,/,
     + 5x,'q1  = 1/X_t - 1/X_c                               ',e12.5,/,
     + 5x,'a1  = 1/(X_t * X_c)                               ',e12.5,/,
     + 5x,'q2  = 1/Y_t - 1/Y_c                               ',e12.5,/,
     + 5x,'a2  = 1/(Y_t * Y_c)                               ',e12.5,/,
     + 5x,'a12 = 2F_xy/dsqrt(X_t*X_c*Y_t*Y_c)                ',e12.5,/,
     + 5x,'a4  = 1/(S_12 * S_12)                             ',e12.5,/,
     + 5x,'a5  = 1/(S_13 * S_13)                             ',e12.5,/,
     + 5x,'a6  = 1/(S_23 * S_23)                             ',e12.5)
2002  format(
     + 5x,'Extended HASHIN failure criterion (values without sign)',/,
     + 5x,'Versagenskriterium                               ',f12.1,/,
     + 5x,'X_t   Bruchspannung Zug in Faserrichtung         ',e12.5,/,
     + 5x,'X_c   Bruchspannung Druck in Faserrichtung       ',e12.5,/,
     + 5x,'Y_t   Bruchspannung Zug quer zur Faserrichtung   ',e12.5,/,
     + 5x,'Y_c   Bruchspannung Druck quer zur Faserrichtung ',e12.5,/,
     + 5x,'S_12  Bruchschubspannung 12                      ',e12.5,/,
     + 5x,'S_13  Bruchschubspannung 13                      ',e12.5,/,
     + 5x,'S_23  Bruchschubspannung 23                      ',e12.5,/,
     + 5x,'qxt = 1/(X_t  *  X_t)                            ',e12.5,/,
     + 5x,'qxc = 1/(X_c  *  X_c)                            ',e12.5,/,
     + 5x,'qyt = 1/(Y_t  *  Y_t)                            ',e12.5,/,
     + 5x,'qyc = 1/(Y_c  *  Y_c)                            ',e12.5,/,
     + 5x,'qs12= 1/(S_12 * S_12)                            ',e12.5,/,
     + 5x,'qs13= 1/(S_13 * S_13)                            ',e12.5,/,
     + 5x,'qs23= 1/(S_23 * S_23)                            ',e12.5)
2003  format(
     + 5x,'Puck criterion (value without dimension)',/,
     + 5x,'Versagenkriterium                                 ',f12.1,/,
     + 5x,'Zugfestigkeit in Faserrichtung                    ',e12.5,/,
     + 5x,'Druckfestigkeit in Faserrichtung                  ',e12.5,/,
     + 5x,'Zugfestigkeit quer zur Faserrichtung              ',e12.5,/,
     + 5x,'Druckfestigkeit quer zur Faserrichtung            ',e12.5,/,
     + 5x,'Festigkeit quer / längs                           ',e12.5,/,
     + 5x,'Neigungsparameter                                 ',e12.5,/,
     + 5x,'Neigungsparameter                                 ',e12.5,/,
     + 5x,'Widerstand der Bruchfläche                        ',e12.5,/,
     + 5x,'Neigungsparameter                                 ',e12.5,/,
     + 5x,'max. Schubspannung                                ',e12.5)
2004  format(
     + 5x,'Cuntze criterion (value without dimension)',/,
     + 5x,'Versagenkriterium                                 ',f12.1,/,
     + 5x,'Zugfestigkeit in Faserrichtung                    ',e12.5,/,
     + 5x,'Druckfestigkeit in Faserrichtung                  ',e12.5,/,
     + 5x,'Zugfestigkeit quer zur Faserrichtung              ',e12.5,/,
     + 5x,'Druckfestigkeit quer zur Faserrichtung            ',e12.5,/,
     + 5x,'Festigkeit quer / längs                           ',e12.5,/,
     + 5x,'Parameter IFF2                                    ',e12.5,/,
     + 5x,'Parameter IFF3                                    ',e12.5,/,
     + 5x,'Interaktionsparameter                             ',e12.5,/,
     + 5x,'taumin (IFF2 - Zug)                               ',e12.5,/,
     + 5x,'taumax (IFF2 - Druck)                             ',e12.5) 
3001  format(
     + 5x,'Degradation model                                ',f12.1,/)
3002  format(
     + 5x,'Degradation model                                ',f12.1,/)
3003  format(
     + 5x,'Degradation model                                ',f12.1,/,
     + 5x,'Degradation parameter by Puck (1) / Knops (2)     ',f12.1)
3004  format(
     + 5x,'Degradation model                                ',f12.1)
 600  format(1x,'--- !! UNKNOWN FAILURE TYPE !! ---',/,1x,
     +'your failure type input:',i5)
 601  format(1x,'--- !! INCOMPATIBLE OR UNKNOWN DEGRADATION TYPE !! ---'
     +,/,1x,'your degradation type input:',i5)
c
      return
      end  
c
c
      subroutine mod3d15(cmat,cappa,d,eps,idam,fi,isw,h1,h2,hdam,plout)
c----------------------------------------------------------------------
c     Elasticity matrix modification for update according to stresses
c     with failure criterion and degradation model 
c.... if failure criterion was fullfilled earlier: E1=E2, E2=0.1*E2
c.... hdam  =>  thicknesses of damaged layers (for M, FM, M+FM, or F)
c----------------------------------------------------------------------
      USE iofile
      implicit double precision(a-h,o-z)
      dimension cmat(6,6),h1(*),h2(*),eps(6),d(*),hdam(6),plout(10)
c
c.... layer data            
      pi = 4.0d0 * datan(1.0d0)
c.... material parameter input
      e1  = d(1)
      e2  = d(2)
      v12 = d(3)
      g12 = d(4)  
      g23 = d(5)  
c
      v21   = v12*e2/e1     
      v23   = 0.5d0*d(2)/d(5)-1.d0     
      cappa  = (5.d0/6.d0)
c
c.... damage parameter input
      iftyp  = d(7)
      idtyp  = d(23)
      idam   = 0
      fi     = 0.d0
c
c...    data for actual layer
        phig = d(6)
        phi = phig / 180.0d0 * pi
        c  = dcos(phi)
        c2 = c * c
        c3 = c * c * c
        c4 = c**4
        s  = dsin(phi)
        s2 = s * s
        s3 = s * s * s
        s4 = s**4
c
c....   local elasticity matrix of layer including damage, modification
        if(iftyp.eq.0) then
          call elast3d0(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,g12,g23,
     +     cappa)          
        else if(iftyp.eq.1) then
          call elast3d1(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,
     +     g12,g23,cappa,eps,d,h1,h2,iftyp,c2,s2,c,s,idam)
        else if(iftyp.eq.2) then
          call elast3d2(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,
     +     g12,g23,cappa,eps,d,h1,h2,idtyp,c2,s2,c,s,idam)
        else if(iftyp.eq.3) then
          call elast3d3(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,g12,
     *         g23,cappa,eps,d,h1,h2,idtyp,c2,s2,c,s)
        else if(iftyp.eq.4) then
          call elast3d4(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,g12,
     *         g23,cappa,eps,d,h1,h2,idtyp,c2,s2,c,s)
        end if
c
        if(iftyp.eq.0.or.iftyp.eq.1.or.iftyp.eq.2.or.iftyp.eq.3.or.
     +     iftyp.eq.4) then 
c
c....   calculate global material parameters (x-y)
        cmat(1,1) = c4*h11l + 2*c2*s2*(h12l+2*h33l) + s4*h22l
        cmat(1,2) = c2*s2*(h11l+h22l-4*h33l) + (c4+s4)*h12l
        cmat(1,4) = c3*s*(h11l-h12l-2*h33l) - s3*c*(h22l-h12l-2*h33l)
        cmat(2,2) = s4*h11l + 2*c2*s2*(h12l+2*h33l) + c4*h22l
        cmat(2,4) = s3*c*(h11l-h12l-2*h33l) - c3*s*(h22l-h12l-2*h33l)
        cmat(4,4) = c2*s2*(h11l+h22l-2*h12l-2*h33l) + (c4+s4)*h33l
        q11g = c2*q11l + s2*q22l
        q12g = c*s*(q11l-q22l)
        q22g = s2*q11l + c2*q22l
c
c.... local elasticity matrix, coefficient
       cc  = 1.0d0/(1.0d0+v23)/(1.0d0-v23-2.0d0*v12*v21)
       c22 = cc*e2*(1.d0-v12*v21)
       c12 = cc*e2*v12*(1.0d0+v23)
       c13 = c12
       c23 = cc*e2*(v23+v12*v21) 
c
c....   Normal parts
        cmat(3,3) = cmat(2,2)             ! Volumetrisches Schalenelement
        cmat(5,5) = q11g
        cmat(6,6) = q22g
c
c....   tangential parts
        cmat(1,3) = c2*c13 + s2*c23 ! Volumetrisches Schalenelement
        cmat(2,3) = s2*c13 + c2*c23 ! Volumetrisches Schalenelement
        cmat(3,4) = c*s*(c13-c23)   ! Volumetrisches Schalenelement
        cmat(5,6) = q12g
c
c....   3 parts alternative
        cmat(3,3) = e2
        cmat(1,3) = 0
        cmat(2,3) = 0
        cmat(3,4) = 0
c
c.... symmetry
       cmat(2,1) = cmat(1,2)
       cmat(3,1) = cmat(1,3)
       cmat(3,2) = cmat(2,3)
       cmat(4,1) = cmat(1,4)
       cmat(4,2) = cmat(2,4)
       cmat(4,3) = cmat(3,4)
       cmat(6,5) = cmat(5,6)
c
       if(isw.eq.4.or.isw.eq.8) then
         call hdam4515(d,h2,hdam)
         do i=1,6
           plout(i) = hdam(i)
         enddo
       endif
      endif
c
      return
      end
c
c---------------------------------------------------------------------
c-------------subroutines elasts from down here-----------------------
c---------------------------------------------------------------------
c
c
      subroutine elast3d0(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,g12,
     + g23,cappa)
c----------------------------------------------------------------------
c     Calculate local elasticity matrix including damage
c     Degradation type = 1   =>   AERNNOVA model
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
c.... undamaged
      e1q  = e1
      e2q  = e2
      f    = 1.d0/(1.d0-v12*v12*e2q/e1q)
      h11l = f*e1q
      h12l = f*v12*e2q
      h22l = f*e2q
      h33l = g12
      q11l = g12*cappa
      q22l = g23*cappa
c      
      return
      end
c
c
      subroutine elast3d1(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,
     +  g12,g23,cappa,eps,d,h1,h2,idtyp,c2,s2,c,s,idam)
c----------------------------------------------------------------------
c     Calculate local elasticity matrix including damage
c     Degradation type = 1   =>   AERNNOVA model
c     If failure then
c     E1 => E2   ;   E2 => 0.1*E2
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension epsg(5),epsl(5),eps(6),d(*),h1(*),h2(*)
c
       idtyp=d(23)
c
c.... undamaged
      e1q  = e1
      e2q  = e2
      f    = 1.d0/(1.d0-v12*v12*e2q/e1q)
      h11l = f*e1q
      h12l = f*v12*e2q
      h22l = f*e2q
      h33l = g12
      q11l = g12*cappa
      q22l = g23*cappa
c      
c.... tolerance for undamaged state (fdam<=1+tol)
      tol  = 1.d-08 + 1.d0
c
      idv1=h1(1) ! failure criterion: 0=undam. 1=dam
      idv2=0
c        
      if(idv1.eq.1) then ! damaged
        if(idtyp.eq.1) then   ! Degradation model 1
           e1q  = e2
            e2q  = 0.1d0*e2
            f    = 1.d0/(1.d0-v12*v12*e2q/e1q)
            g12q = 0.1d0*g12
            g23q = 0.1d0*g23
c
            h11l = f*e1q
            h12l = f*v12*e2q
            h22l = f*e2q
            h33l = g12q
            q11l = g12q*cappa
            q22l = g23q*cappa
        end if
      else ! check if now damaged 
c
c....   global strains at position zs
        epsg(1) = eps(1)
        epsg(2) = eps(2)
        epsg(3) = eps(4)
        epsg(4) = eps(5)
        epsg(5) = eps(6)
c
c....   local  strains at position zs:  eps_l = t_e * eps_g
        epsl(1) =      c2*epsg(1) +     s2*epsg(2) +   s *c *epsg(3)
        epsl(2) =      s2*epsg(1) +     c2*epsg(2) -   s *c *epsg(3)
        epsl(3) = -2*s*c *epsg(1) + 2*s*c *epsg(2) + (c2-s2)*epsg(3)
        epsl(4) =      c *epsg(4) +     s *epsg(5) 
        epsl(5) =     -s *epsg(4) +     c *epsg(5) 
c
c....   local stresses
        s11 =h11l*epsl(1) + h12l*epsl(2) 
        s22 =h12l*epsl(1) + h22l*epsl(2) 
        s12 =h33l*epsl(3)                 
        s13 =q11l*epsl(4)                 
        s23 =q22l*epsl(5)
c        
        if(idtyp.eq.1) then
c....     Tsai-Wu-criterion        
c         q1*s11+a1*s11*s11+a12*s11*s22+q2*s22+a2*s22*s22
c         +a4*s12*s12+a5*s13*s13+a6*s23*s23<1
          q1 =d(15)
          a1 =d(16)
          q2 =d(17)
          a2 =d(18)
          a12=d(19)
          a4 =d(20)
          a5 =d(21)
          a6 =d(22)
c
          fdam = q1*s11+a1*s11*s11+a12*s11*s22+q2*s22+a2*s22*s22
     1         + a4*s12*s12 + a5*s13*s13+ a6*s23*s23
c
          if(fdam.gt.tol) then
            if(idtyp.eq.1) then   ! Degradation model 1
c....         modify material data  
              e1q  = e2
              e2q  = 0.1d0*e2
              f    = 1.d0/(1.d0-v12*v12*e2q/e1q)
              g12q = 0.1d0*g12
              g23q = 0.1d0*g23
            end if
c....       modify local elasticity matrix of layer
            h11l = f*e1q
            h12l = f*v12*e2q
            h22l = f*e2q
            h33l = g12q
            q11l = g12q*cappa
            q22l = g23q*cappa
c
c....       set failure index to 1
            idv2=1
          end if
        end if 
      end if
c
      h2(1)= max(idv1,idv2)  ! no undamaging possible
      ih2=h2(1)
      idam=max(idam,ih2)
c       
      return
      end
c
c
c
            subroutine elast3d2(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,
     +  g12,g23,cappa,eps,d,h1,h2,idtyp,c2,s2,c,s,idam)
c----------------------------------------------------------------------
c     Calculate local elasticity matrix including damage
c     Degradation type = 1   =>   AERNNOVA model
c                      = 2   =>   edited CRC-ACS  model
c
c     idv1     =>  internal input history parameter
c     idv2     =>  internal output history parameter
c     idam     =>  damage mode identifier
c     hdam(i)  =>  height of damaged layers (for M, FM, M+FM, F)
c                  i=1 => Matrix cracking (M)
c                  i=2 => Fiber-Matrix Shearing (FM)
c                  i=3 => Matrix cracking + Fiber-matrix shearing (M+FM)
c                  i=4 => Fiber Fracture (F)
c
c |IDAM| FAILURE MODE                |ABBREV|   PARAMETERS REDUCED    |
c |  0 | no failure                  |  --  |            -            |
c |  1 | matrix                      |   M  |     C12,C22             |
c |  2 | fiber-matrix shear          |  FM  |             C55,C66     |
c |  3 | matrix + fiber-matrix shear | MFM  |     C12,C22,C55,C66     |
c |  4 | fiber breakage              |  F   | C11,C12,C22,C55,C66     |
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension d(*),eps(6),h1(3),h2(3),epsg(5),epsl(5)
c
      idv1 = h1(1)
      idvt = h1(2)
      idvc = h1(3)
      idv2 = 0
      if(idtyp.eq.2) then
        fact = 0.01d0
        factf = 0.01d0
      end if
c
c.... tolerance for undamaged state (fdam<=1+tol)
      tol  = 1.d-08 + 1.d0
c
c.... set local elasticity matrix
      f = 1.d0/(1.d0-v12*v12*e2/e1)
        if(idv1.eq.0) then             ! no damage
        h11l = f*e1
        h12l = f*v12*e2
        h22l = f*e2
        h33l = g12
        q11l = g12*cappa
        q22l = g23*cappa
      else if(idv1.eq.1) then        ! M
        h11l = f*e1
        h12l = fact*f*v12*e2
        h22l = fact*f*e2
        h33l = g12              
        q11l = g12*cappa
        q22l = g23*cappa        
      else if(idv1.eq.2) then        ! FM
        h11l = f*e1
        h12l = f*v12*e2
        h22l = f*e2
        h33l = fact*g12
        q11l = fact*g12*cappa
        q22l = g23*cappa
      else if(idv1.eq.3) then        ! MFM
        h11l = f*e1
        h12l = fact*f*v12*e2
        h22l = fact*f*e2
        h33l = fact*g12
        q11l = fact*g12*cappa
        q22l = g23*cappa
      else if(idv1.eq.4) then        ! F
        h11l = factf*f*e1
        h12l = factf*f*v12*e2
        h22l = factf*f*e2
        h33l = factf*g12
        q11l = factf*g12*cappa
        q22l = factf*g23*cappa
      end if
c
c.... check for damage if not yet broken (failure = fiber breakage = F)
      if(idtyp.eq.2) then
c
        if(idv1.lt.4) then
c.....  global strains at position zs
        epsg(1) = eps(1)
        epsg(2) = eps(2)
        epsg(3) = eps(4)
        epsg(4) = eps(5)
        epsg(5) = eps(6)
c.....  local  strains at position zs:  eps_l = t_e * eps_g
        epsl(1) =      c2*epsg(1) +     s2*epsg(2) +   s *c *epsg(3)
        epsl(2) =      s2*epsg(1) +     c2*epsg(2) -   s *c *epsg(3)
        epsl(3) = -2*s*c *epsg(1) + 2*s*c *epsg(2) + (c2-s2)*epsg(3)
        epsl(4) =      c *epsg(4) +     s *epsg(5)
        epsl(5) =     -s *epsg(4) +     c *epsg(5)
c.....  local stresses
        s11 =h11l*epsl(1) + h12l*epsl(2)
        s22 =h12l*epsl(1) + h22l*epsl(2)
        s12 =h33l*epsl(3)
        s13 =q11l*epsl(4)
        s23 =q22l*epsl(5)
c
c.....    Extended HASHIN criteria
c.....    parameters for extended Hashin criterion
          yc   = d(11)
          qxt  = d(15)
          qxc  = d(16)
          qyt  = d(17)
          qs12 = d(19)
          qs13 = d(20)
          qs23 = d(21)
c.....    check for transverse matrix cracking
          if(idv1.le.2) then
            if(s22.ge.0.d0) then
              fdam = s22*s22*qyt
     +             + s12*s12*qs12 + s13*s13*qs13 + s23*s23*qs23
              if(fdam.gt.tol) idvt = 1
            else if(s22.lt.0.d0) then
              fdam = s22/yc*(yc*yc*qs23/4.d0-1.d0) + s22*s22*qs23/4.d0
     +             + s12*s12*qs12 + s13*s13*qs13 + s23*s23*qs23
              if(fdam.gt.tol) idvc = 1
            end if
            if(fdam.gt.tol) then
              if(idv1.eq.0)  idv2 = 1
              if(idv1.eq.2)  idv2 = 3
            end if
          end if
c       
c......   check for fiber-matrix shearing
          if(idv1.eq.0.or.idv1.eq.1) then
            if(s11.ge.0.d0) fdam = s12*s12*qs12+s13*s13*qs13
            if(s11.lt.0.d0) fdam = s12*s12*qs12+s13*s13*qs13+s11*s11*qxc
            if(fdam.gt.tol) then
              if(idv1.eq.0)  idv2 = 2
              if(idv1.eq.1)  idv2 = 3
            end if
          end if
c       
c......   check for fiber breakage
          if(s11.ge.0.d0)  fdam = s11*s11*qxt
          if(s11.lt.0.d0)  fdam = s11*s11*qxc
          if(fdam.gt.tol)  idv2 = 4
c
        endif
      end if
c
      idam=max(idv1,idv2)      
c.... modify local elasticity matrix of layer
      if(idam.eq.1) then
        h11l = f*e1
        h12l = fact*f*v12*e2
        h22l = fact*f*e2
        h33l = g12            
        q11l = g12*cappa
        q22l = g23*cappa       
      else if(idam.eq.2) then
        h11l = f*e1
        h12l = f*v12*e2
        h22l = f*e2
        h33l = fact*g12
        q11l = fact*g12*cappa
        q22l = g23*cappa
      else if(idam.eq.3) then
        h11l = f*e1
        h12l = fact*f*v12*e2
        h22l = fact*f*e2
        h33l = fact*g12
        q11l = fact*g12*cappa
        q22l = g23*cappa
      else if(idam.eq.4) then
        h11l = factf*f*e1
        h12l = factf*f*v12*e2
        h22l = factf*f*e2
        h33l = factf*g12
        q11l = factf*g12*cappa
        q22l = factf*g23*cappa
      end if
c
c.... assign history variable
      h2(1) = idam
      h2(2) = idvt
      h2(3) = idvc
c       
      return
      end
c
c
      subroutine elast3d3(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,g12,
     *         g23,cappa,eps,d,h1,h2,idtyp,c2,s2,c,s)
c----------------------------------------------------------------------
c     Calculate local elasticity matrix including damage
c     In-plane damage detection type = 3 => Puck criterion
c     Degradation type = 3               => Chang & Lessard degradation
c     Degradation type = 4               => Puck degradation model
c
c     | IDAM | FAILURE MODE                | PARAMETERS REDUCED    |
c     |   0  | no failure                  |           -           |
c     |   1  | matrix cracking             |        C44,C66        |
c     |   2  | fibre fracture              |C11,C12,C22,C44,C55,C66|
c     
c     hdam(i)  =>  thickness of damaged layers / total thickness
c                  i=1 => Matrix cracking: General
c                  i=2 => Matrix cracking: Mode A
c                  i=3 => Matrix cracking: Mode B
c                  i=4 => Matrix cracking: Mode C
c                  i=5 => Fiber Fracture
c                  i=6 => Fracture angle phi 
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      real*8    lim,phic
      dimension d(*),epsg(5),epsl(5),eps(6),h1(*),h2(*)
c
c.... read history variables
      idv1 = h1(1)
      idva = h1(2)
      idvb = h1(3)
      idvc = h1(4)
      emin = h1(5)
      phic = h1(6)
      if(idtyp.eq.4)then
        etae = h1(7) ; if(etae.eq.0.d0) etae=1.d0
        etag = h1(8) ; if(etag.eq.0.d0) etag=1.d0
        etaneg = h1(9) ; if(etaneg.eq.0.d0) etaneg=1.d0
      endif
c
c.... set constant values
      tol0  = (d(24)+d(25))/150000.d0
      tol1  = 1.d0+1.d-06
      fact  = 0.4d0
      factf = 0.01d0
      etaae = 1.d0
      etaag = 1.d0
      etab  = 1.d0
      etac  = 1.d0
c
c.... set values initially to zero
      idv2  = 0
      eiffa = 0.d0
      eiffb = 0.d0
      eiffc = 0.d0
      fdam  = 0.d0
c
c.... global strains
      epsg(1) = eps(1)
      epsg(2) = eps(2)
      epsg(3) = eps(4)
      epsg(4) = eps(5)
      epsg(5) = eps(6)
c
c.... local strains at position zs: eps_1 = t_e * eps_g
      epsl(1) =      c2*epsg(1) +     s2*epsg(2) +   s *c *epsg(3)
      epsl(2) =      s2*epsg(1) +     c2*epsg(2) -   s *c *epsg(3)
      epsl(3) = -2*s*c *epsg(1) + 2*s*c *epsg(2) + (c2-s2)*epsg(3)
      epsl(4) =      c *epsg(4) +     s *epsg(5)
      epsl(5) =     -s *epsg(4) +     c *epsg(5)
c
c.... emin
      if(epsl(2).lt.emin) emin=epsl(2)
c
c.... set local elasticity matrix without damage
      f = 1.d0/(1.d0-v12**2*e2/e1)
      h11l = f*e1
      h12l = f*e2*v12
      h22l = f*e2
      h33l = g12
      q11l = g12*cappa
      q22l = g23*cappa
c.... local stresses without damage
      s22q = h12l*epsl(1) + h22l*epsl(2)
      s12q = h33l*epsl(3)
c
c.... set local elasticity matrix
      if(idtyp.eq.3)then
        etae=fact ; etag=fact ; etaneg=fact
      endif
      if(idv1.eq.0)then
        h11l = f*e1
        h12l = f*e2*v12
        h22l = f*e2
        h33l = g12
        q11l = g12*cappa
        q22l = g23*cappa
      elseif(idv1.eq.1)then
        if(idva.eq.1.and.epsl(2).gt.emin)then
          h11l = f*e1
          h12l = etae*f*e2*v12
          h22l = etae*f*e2
          h33l = etag*g12
          q11l = g12*cappa
          q22l = etag*g23*cappa
        elseif(idva.eq.1.and.epsl(2).eq.emin)then
          h11l = f*e1
          h12l = f*e2*v12
          h22l = f*e2
          h33l = etag*g12
          q11l = g12*cappa
          q22l = etag*g23*cappa
        elseif(idvb.eq.1.or.idvc.eq.1)then
          h11l = f*e1
          h12l = f*e2*v12
          h22l = f*e2
          h33l = etaneg*g12
          q11l = g12*cappa
          q22l = etaneg*g23*cappa
        endif
      elseif(idv1.eq.2)then
        h11l = factf*f*e1
        h12l = factf*f*e2*v12
        h22l = factf*f*e2
        h33l = factf*g12
        q11l = factf*g12*cappa
        q22l = factf*g23*cappa
      endif
c.... local stresses
      s11 = h11l*epsl(1) + h12l*epsl(2)
      s22 = h12l*epsl(1) + h22l*epsl(2)
      s12 = h33l*epsl(3)
c
c.... read material strength parameters
      rt11 = d(8)
      rc11 = d(9)
      rt22 = d(10)
      rc22 = d(11)
      r21  = d(12)
      pt21 = d(13)
      pc21 = d(14)
      ra22 = d(15)
      p22  = d(16)
      t21c = d(17)
c
c.... check for fibre fracture
      if(s11.ge.0.d0) rs11=abs(s11/rt11)
      if(s11.lt.0.d0) rs11=abs(s11/rc11)
      fw=1.d0
      if(rs11.ge.tol1) then
        idv2=2 ; idva=0 ; idvb=0 ; idvc=0
      endif
c
c.... check for matrix cracking
      if(idv1.ne.2.and.idv2.ne.2)then
        if(s22.ge.tol0.and.idvb.eq.0.and.idvc.eq.0)then
          fdam=1.d0/fw*(dsqrt((abs(s12)/r21)**2+(1.d0-pt21*rt22/r21)
     +         **2*(abs(s22)/rt22)**2)+pt21*abs(s22)/r21)
          if(idva.eq.0.and.idvb.eq.0.and.idvc.eq.0.and.fdam.gt.tol1)then
            idv2=1 ; idva=1
          endif
          if(idva.eq.1)then
            eiffa=1.d0/fw*(dsqrt((abs(s12q)/r21)**2+(1.d0-pt21*rt22/
     +            r21)**2*(abs(s22q)/rt22)**2)+pt21*abs(s22q)/r21)
            if(eiffa.gt.tol1) call puckpara045(1,eiffa,etaae,etaag)
          endif
        elseif(s22.le.-tol0.and.idva.eq.0)then
          srat=abs(s12/s22)
          lim=abs(t21c/ra22)
          if(srat.gt.lim.and.idvc.eq.0)then
            fdam=1.d0/fw*(1.d0/r21*(dsqrt(abs(s12)**2+(pc21*abs(s22))
     +           **2)+pc21*abs(s22)))
            if(idva.eq.0.and.idvb.eq.0.and.idvc.eq.0.and.fdam.gt.tol1)
     +       then
              idv2=1 ; idvb=1
            endif
            if(idvb.eq.1)then
              eiffb=1.d0/fw*(1.d0/r21*(dsqrt(abs(s12q)**2+(pc21*
     +              abs(s22q))**2)+pc21*abs(s22q)))
            endif
            if(eiffb.gt.tol1)then
              call puckpara045(2,eiffb,dummy,etab)
              if(s12q.eq.0.d0) rho=acos(0.d0)
              if(s12q.ne.0.d0) rho=atan(abs(s22q/s12q))
              etab=etab*(cos(rho))**2+(sin(rho))**2
           endif
          elseif(srat.le.lim.and.idvb.eq.0)then
            fdam=1.d0/fw*(((abs(s12)/(2.d0*(1.d0+p22)*r21))**2+
     +              (abs(s22)/rc22)**2)*rc22/abs(s22))
            if(idva.eq.0.and.idvb.eq.0.and.idvc.eq.0.and.fdam.gt.tol1)
     +       then
              idv2=1 ; idvc=1
            endif
            if(idvc.eq.1)then
              eiffc=1.d0/fw*(((abs(s12q)/(2.d0*(1.d0+p22)*r21))**2+
     +              (abs(s22q)/rc22)**2)*rc22/abs(s22q))
            endif
            if(eiffc.gt.tol1)then  
              call puckpara045(3,eiffc,dummy,etac)
              if(s12q.eq.0.d0) rho=acos(0.d0)
              if(s12q.ne.0.d0) rho=atan(abs(s22q/s12q))
              etac=etac*(cos(rho))**2+(sin(rho))**2
              if(phic.eq.0.d0) phic=acos(dsqrt(abs(fw*ra22/s22)))
     +                              *180.d0/2.d0/acos(0.d0)
            endif
          else
            continue
          endif
        else
          fdam=1.d0/fw*(abs(s12)/r21)
          if(idva.eq.0.and.idvb.eq.0.and.idvc.eq.0.and.fdam.gt.tol1)then
              idv2=1 ; idva=0 ; idvb=1 ; idvc=0
          endif
          if(idvb.eq.1)then
            eiffb=1.d0/fw*(abs(s12q)/r21)
          endif
          if(eiffb.gt.tol1)then
            call puckpara045(2,eiffb,dummy,etab)
            if(s12q.eq.0.d0) rho=acos(0.d0)
            if(s12q.ne.0.d0) rho=atan(abs(s22q/s12q))
            etab=etab*(cos(rho))**2+(sin(rho))**2
          endif
        endif
      endif
c
c.... reset local elasticity matrix
      idam=max(idv1,idv2)
c
      if(idtyp.eq.3)then
        etae=fact ; etag=fact ; etaneg=fact
      elseif(idtyp.eq.4)then
        etae=min(etae,etaae) ; etag=min(etag,etaag)
        etaneg=min(etaneg,etab,etac)
      endif
      if(idam.eq.0)then
        h11l = f*e1
        h12l = f*e2*v12
        h22l = f*e2
        h33l = g12
        q11l = g12*cappa
        q22l = g23*cappa
      elseif(idam.eq.1)then
        if(idva.eq.1.and.epsl(2).gt.emin)then
          h11l = f*e1
          h12l = etae*f*e2*v12
          h22l = etae*f*e2
          h33l = etag*g12
          q11l = g12*cappa
          q22l = etag*g23*cappa
        elseif(idva.eq.1.and.epsl(2).eq.emin)then
          h11l = f*e1
          h12l = f*e2*v12
          h22l = f*e2
          h33l = etag*g12
          q11l = g12*cappa
          q22l = etag*g23*cappa
        elseif(idvb.eq.1.or.idvc.eq.1)then
          h11l = f*e1
          h12l = f*e2*v12
          h22l = f*e2
          h33l = etaneg*g12
          q11l = g12*cappa
          q22l = etaneg*g23*cappa
        endif
      elseif(idam.eq.2)then
        h11l = factf*f*e1
        h12l = factf*f*e2*v12
        h22l = factf*f*e2
        h33l = factf*g12
        q11l = factf*g12*cappa
        q22l = factf*g23*cappa
      endif
c
c.... write history values
      h2(1) = idam
      h2(2) = idva
      h2(3) = idvb
      h2(4) = idvc
      h2(5) = emin
      h2(6) = phic
      if(idtyp.eq.4)then
        h2(7) = etae
        h2(8) = etag
        h2(9) = etaneg
      endif
c
      return
      end
c
c
      subroutine elast3d4(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,g12,
     *         g23,cappa,eps,d,h1,h2,idtyp,c2,s2,c,s)
c----------------------------------------------------------------------
c     Calculate local elasticity matrix including damage
c     In-plane damage detection type = 4 => Cuntze criterion
c     Degradation type = 0               => No degradation
c     Degradation type = 3               => Chang & Lessard degradation
c     Degradation type = 5               => Puck degradation model
c
c     | IDAM | FAILURE MODE                  | PARAMETERS REDUCED    |
c     |   0  | no failure                    |           -           |
c     |   1  | matrix cracking (compression) |        C44,C66        |
c     |   1  | matrix cracking (tension)     |    C12,C22,C44,C66    |
c     |   2  | fibre fracture                |C11,C12,C22,C44,C55,C66|
c     
c     hdam(i)  =>  thickness of damaged layers / total thickness
c                  i=1 => Matrix cracking: General
c                  i=2 => Matrix cracking: IFF1
c                  i=3 => Matrix cracking: IFF2a (tension)
c                  i=4 => Matrix cracking: IFF2b (compression)
c                  i=5 => Matrix cracking: IFF3
c                  i=6 => Fiber Fracture
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension d(*),epsg(5),epsl(5),h1(9),h2(9),eps(6)
c
c     read history variables
      idam = h1(1)
      iff1 = h1(2)
      iff2 = h1(3)
      iff3 = h1(4)
      emin = h1(5)
      if(idtyp.eq.5)then
        etae = h1(6) ; if(etae.eq.0.d0) etae=1d0
        etag = h1(7) ; if(etag.eq.0.d0) etag=1d0
        etaneg = h1(8) ; if(etaneg.eq.0.d0) etaneg=1d0
      endif
c
c.... set constant values
      tol0  = 1.d-06
      fact  = 0.01d0
      fred1e = 1.d0
      fred2ae = 1.d0
      fred1g = 1.d0
      fred2ag = 1.d0
      fred2bg = 1.d0
      fred3g = 1.d0
c
c.... set values initially to zero
      eff1 = 0.d0
      eff2 = 0.d0
      eff3 = 0.d0
      eff4 = 0.d0
      eff5 = 0.d0
      effr = 0.d0
      uff1 = 0.d0
      uff2 = 0.d0
      uff3 = 0.d0
      uff4 = 0.d0
      uff5 = 0.d0
      uffr = 0.d0
c
c.... global strains
        epsg(1) = eps(1)
        epsg(2) = eps(2)
        epsg(3) = eps(4)
        epsg(4) = eps(5)
        epsg(5) = eps(6)
c
c.... local strains at position zs: eps_1 = t_e * eps_g
      epsl(1) =      c2*epsg(1) +     s2*epsg(2) +   s *c *epsg(3)     
      epsl(2) =      s2*epsg(1) +     c2*epsg(2) -   s *c *epsg(3) 
      epsl(3) = -2*s*c *epsg(1) + 2*s*c *epsg(2) + (c2-s2)*epsg(3) 
      epsl(4) =      c *epsg(4) +     s *epsg(5)                   
      epsl(5) =     -s *epsg(4) +     c *epsg(5)                 
c  
c.... emin
      if(epsl(2).lt.emin) emin=epsl(2)
c
c.... set local elasticity matrix without damage
      f = 1.d0/(1.d0-v12**2*e2/e1)
      h11l = f*e1
      h12l = f*e2*v12
      h22l = f*e2
      h33l = g12
      q11l = g12*cappa
      q22l = g23*cappa
c
c.... local stresses without damage
      s11q = h11l*epsl(1) + h12l*epsl(2)
      s22q = h12l*epsl(1) + h22l*epsl(2)
      s12q = h33l*epsl(3)
      s13q = q11l*epsl(4)
      s23q = q22l*epsl(5)
c
c.... set invariants without damage
      binv2 = s22q
      binv3 = s12q**2.d0 + s13q**2.d0
      binv4 = s22q**2.d0 + 4.d0*s23q**2.d0
      binv5 = s22q*(s13q**2.d0-s12q**2.d0) - 4.d0*s12q*s13q*s23q
c
c.... set local elasticity matrix
      if(idtyp.eq.0) then
        h11l = e1*f
        h12l = e2*v12*f
        h22l = e2*f
        h33l = g12
        q11l = cappa*g12
        q22l = cappa*g23
      elseif(idtyp.eq.3.or.idtyp.eq.5) then
        if(idtyp.eq.3) then
        etae=fact ; etae=fact ; etaneg=fact
        endif
        if(idam.eq.0) then
          h11l = e1*f
          h12l = e2*v12*f
          h22l = e2*f
          h33l = g12
          q11l = cappa*g12
          q22l = cappa*g23
        elseif(idam.eq.1) then
          if((iff1.eq.1.or.iff2.eq.2).and.epsl(2).gt.emin) then     !Matrix-Zug
            h11l = e1*f
            h12l = etae*e2*v12*f      !etae*e2*v12*f
            h22l = etae*e2*f          !etae*e2*f   
            h33l = etag*g12           !etag*g12   
            q11l = cappa*g12
            q22l = etag*cappa*g23     !etag*cappa*g23 
          elseif((iff1.eq.1.or.iff2.eq.2).and.epsl(2).eq.emin) then
            h11l = e1*f
            h12l = e2*v12*f
            h22l = e2*f
            h33l = etaneg*g12
            q11l = cappa*g12
            q22l = etaneg*cappa*g23
          else                                                      !Matrix-Druck
            h11l = e1*f
            h12l = e2*v12*f
            h22l = e2*f
            h33l = etaneg*g12
            q11l = cappa*g12
            q22l = etaneg*cappa*g23
          endif
        elseif(idam.eq.2) then
          h11l = fact*e1*f
          h12l = fact*e2*v12*f
          h22l = fact*e2*f
          h33l = fact*g12
          q11l = fact*cappa*g12
          q22l = fact*cappa*g23
        else 
          stop
        endif
      else
        stop
      endif
c
c.... set local stresses
      s11l = h11l*epsl(1) + h12l*epsl(2)
      s22l = h12l*epsl(1) + h22l*epsl(2)
      s12l = h33l*epsl(3)
      s13l = q11l*epsl(4)
      s23l = q22l*epsl(5)
c
c.... set invariants
      ainv1 = s11l
      ainv2 = s22l
      ainv3 = s12l**2.d0 + s13l**2.d0
      ainv4 = s22l**2.d0 + 4.d0*s23l**2.d0
      ainv5 = s22l*(s13l**2.d0-s12l**2.d0) - 4.d0*s12l*s13l*s23l
c
c.... read materialparameters
      rt11 = d(8)
      rc11 = d(9)
      rt22 = d(10)
      rc22 = d(11)
      r12  = d(12)
      b12  = d(13)
      b22  = d(14)
      rop  = d(15)
      tmin = d(16)
      tmax = d(17)
c
c.... check for fibre fracture
      if(idam.ne.2) then
        if(ainv1.gt.0.d0) eff1= ainv1/rt11
        if(ainv1.lt.0.d0) eff2=-ainv1/rc11
        uff1 = eff1
        uff2 = eff2
        if(abs(eff1).gt.1.d0.or.abs(eff2).gt.1.d0) idam=2
      endif
c
c.... check for matrix fracture
      if(idam.ne.2) then
        if(s22l.ge.tol0) then
          if(abs(s12l/s22l).le.abs(tmin/rt22)) then
            eff3 = (ainv2+sqrt(ainv4))/2.d0/rt22                     !dominantes Versagen
            eff4 = (ainv3**(3.d0/2.d0) + b12*(ainv2*ainv3-ainv5))    !benötigt für interaktion
     *             **(1.d0/3.d0)/r12
            if(eff3.lt.0.d0) eff3 = 0.d0
            if(eff4.lt.0.d0) eff4 = 0.d0
            effr = (eff1**rop + eff2**rop + eff3**rop + 
     +              eff4**rop + eff5**rop)**(1/rop)
            if(eff3.ge.1.d0.or.effr.ge.1.d0) iff1 = 1
            if(eff3.ge.1.d0.or.effr.ge.1.d0) idam = 1
            if(idtyp.eq.5) then
              uff3 = (binv2+sqrt(binv4))/2.d0/rt22                     
              uff4 = (binv3**(3.d0/2.d0) + b12*(binv2*binv3-binv5))    
     *               **(1.d0/3.d0)/r12
              if(uff3.lt.0.d0) uff3 = 0.d0
              if(uff4.lt.0.d0) uff4 = 0.d0
              uffr = (uff1**rop + uff2**rop + uff3**rop + 
     +                uff4**rop + uff5**rop)**(1/rop)
c.... fred1
              if(uffr.ge.1.d0) then
                call cuntzepara045(1,uffr,fred1e,fred1g,fred3g)
              endif
            endif
          elseif(abs(s12l/s22l).gt.abs(tmin/rt22)) then
            eff4 = (ainv3**(3.d0/2.d0) + b12*(ainv2*ainv3-ainv5))    !dominantes Versagen
     *             **(1.d0/3.d0)/r12
            eff3 = (ainv2+sqrt(ainv4))/2.d0/rt22                     !benötigt für interaktion
            if(eff3.lt.0.d0) eff3 = 0.d0
            if(eff4.lt.0.d0) eff4 = 0.d0
            effr = (eff1**rop + eff2**rop + eff3**rop + 
     +              eff4**rop + eff5**rop)**(1/rop)
            if(eff4.ge.1.d0.or.effr.ge.1.d0) iff2 = 2
            if(eff4.ge.1.d0.or.effr.ge.1.d0) idam = 1
            if(idtyp.eq.5) then
              uff4 = (binv3**(3.d0/2.d0) + b12*(binv2*binv3-binv5))    !dominantes Versagen
     *               **(1.d0/3.d0)/r12
              uff3 = (binv2+sqrt(binv4))/2.d0/rt22                     !benötigt für interaktion
              if(uff3.lt.0.d0) uff3 = 0.d0
              if(uff4.lt.0.d0) uff4 = 0.d0
              uffr = (uff1**rop + uff2**rop + uff3**rop + 
     +                uff4**rop + uff5**rop)**(1/rop)
c.... fred2a 
              if(uffr.ge.1.d0) then
                call cuntzepara045(1,uffr,fred2ae,fred2ag,fred3g)
              endif
            endif
          endif
        elseif(s22l.lt.tol0) then
          if(abs(s12l/s22l).le.abs(tmax/rc22)) then
            eff5 = ((b22-1.d0)*ainv2+b22*sqrt(ainv4))/rc22             !dominantes Versagen
            if(abs(s12l/s22l).gt.(2*b12*1.01)) then
              eff4 = (ainv3**(3.d0/2.d0) + b12*(ainv2*ainv3-ainv5))    !benötigt für Interaktion
     *               **(1.d0/3.d0)/r12
            endif
            if(eff5.lt.0.d0) eff5 = 0.d0
            if(eff4.lt.0.d0) eff4 = 0.d0
            effr = (eff1**rop + eff2**rop + eff3**rop + 
     +              eff4**rop + eff5**rop)**(1/rop)
            if(eff5.ge.1.d0.or.effr.ge.1.d0) iff3 = 1
            if(eff5.ge.1.d0.or.effr.ge.1.d0) idam = 1
            if(idtyp.eq.5) then
              uff5 = ((b22-1.d0)*binv2+b22*sqrt(binv4))/rc22            
              if(abs(s12q/s22q).gt.(2*b12*1.01)) then
                uff4 = (binv3**(3.d0/2.d0) + b12*(binv2*binv3-binv5))   
     *                 **(1.d0/3.d0)/r12
              endif
              if(uff5.lt.0.d0) eff5 = 0.d0
              if(uff4.lt.0.d0) eff4 = 0.d0
              uffr = (uff1**rop + uff2**rop + uff3**rop + 
     +                uff4**rop + uff5**rop)**(1/rop)
c.... fred3
              if(uffr.ge.1.d0) then 
                call cuntzepara045(2,uffr,dummy,fred3g,dummy)
              endif
            endif
          elseif(abs(s12l/s22l).gt.abs(tmax/rc22)) then
            eff4 = (ainv3**(3.d0/2.d0) + b12*(ainv2*ainv3-ainv5))      !dominantes Versage
     *             **(1.d0/3.d0)/r12
            eff5 = ((b22-1.d0)*ainv2+b22*sqrt(ainv4))/rc22             !benötigt für Interaktion
            if(eff4.lt.0.d0) eff4 = 0.d0
            if(eff5.lt.0.d0) eff5 = 0.d0
            effr = (eff1**rop + eff2**rop + eff3**rop + 
     +              eff4**rop + eff5**rop)**(1/rop)
            if((eff4.ge.1.d0.or.effr.ge.1.d0).and.iff2.ne.2) iff2 = 1
            if((eff4.ge.1.d0.or.effr.ge.1.d0).and.iff2.ne.2) idam = 1
            if(idtyp.eq.5) then
              if(abs(s12q/s22q).gt.(2*b12*1.01)) then
                uff4 = (binv3**(3.d0/2.d0) + b12*(binv2*binv3-binv5))   
     *                 **(1.d0/3.d0)/r12
              endif
              uff5 = ((b22-1.d0)*binv2+b22*sqrt(binv4))/rc22         
              if(uff4.lt.0.d0) uff4 = 0.d0
              if(uff5.lt.0.d0) uff5 = 0.d0
              uffr = (uff1**rop + uff2**rop + uff3**rop + 
     +                uff4**rop + uff5**rop)**(1/rop)
c.... fred2b
              if(uffr.ge.1.d0) then
                call cuntzepara045(2,uffr,dummy,fred2bg,dummy)
              endif
            endif
          endif
        endif 
      endif
c
c.... Abfrage: Suche nach kleinstem Reduktionsfaktor!
c.... reset local elasticity matrix
      if(idtyp.eq.0) then
        h11l = e1*f
        h12l = e2*v12*f
        h22l = e2*f
        h33l = g12
        q11l = cappa*g12
        q22l = cappa*g23
      elseif(idtyp.eq.3.or.idtyp.eq.5) then
        if(idtyp.eq.3) then
        etae=fact ; etag=fact ; etaneg=fact
        elseif(idtyp.eq.5)then
        etae=min(etae,fred1e,fred2ae) ; etag=min(etag,fred1g,fred2ag)
        etaneg=min(etaneg,fred2bg,fred3g)
        endif
        if(idam.eq.0) then
          h11l = e1*f
          h12l = e2*v12*f
          h22l = e2*f
          h33l = g12
          q11l = cappa*g12
          q22l = cappa*g23
        elseif(idam.eq.1) then
          if((iff1.eq.1.or.iff2.eq.2).and.epsl(2).gt.emin) then
            h11l = e1*f
            h12l = etae*e2*v12*f      !etae*e2*v12*f
            h22l = etae*e2*f          !etae*e2*f   
            h33l = etag*g12           !etag*g12   
            q11l = cappa*g12
            q22l = etag*cappa*g23     !etag*cappa*g23 
          elseif((iff1.eq.1.or.iff2.eq.2).and.epsl(2).eq.emin) then
            h11l = e1*f
            h12l = e2*v12*f
            h22l = e2*f
            h33l = etaneg*g12
            q11l = cappa*g12
            q22l = etaneg*cappa*g23
          else
            h11l = e1*f
            h12l = e2*v12*f
            h22l = e2*f
            h33l = etaneg*g12
            q11l = cappa*g12
            q22l = etaneg*cappa*g23
          endif
        elseif(idam.eq.2) then
          h11l = fact*e1*f
          h12l = fact*e2*v12*f
          h22l = fact*e2*f
          h33l = fact*g12
          q11l = fact*cappa*g12
          q22l = fact*cappa*g23
        else
          stop
        endif
      else
        stop
      endif
c
c.... write history variables
      h2(1) = idam
      h2(2) = iff1
      h2(3) = iff2
      h2(4) = iff3
      h2(5) = emin
      if(idtyp.eq.5)then
        h2(6)  = etae
        h2(7)  = etag
        h2(8)  = etaneg
      endif
c
      return
      end
c
c
      subroutine puckpara045(i,eiff,etae,etag)
c----------------------------------------------------------------------
c     Degradation for Puck Model
c     | i |  Counter for Mode A=1 / B=2 / C=3
c----------------------------------------------------------------------
      implicit none
      integer i
      real*8 eiff,etae,etag
      real*8 ic,etar,xi
c
      if(i.eq.1)then
        ic=5.3d0 ; etar=0.03d0 ; xi=1.3d0
        etae=(1.d0-etar)/(1.d0+ic*(eiff-1.d0)**xi)+etar
        ic=0.95d0 ; etar=0.67d0 ; xi=1.17d0
        etag=1
c        etag=(1.d0-etar)/(1.d0+ic*(eiff-1.d0)**xi)+etar
      elseif(i.eq.2)then
        ic=4.d0 ; etar=0.d0 ; xi=2.d0
        etag=(1.d0-etar)/(1.d0+ic*(eiff-1.d0)**xi)+etar
      elseif(i.eq.3)then
        ic=4.d0 ; etar=0.d0 ; xi=2.d0
        etag=(1.d0-etar)/(1.d0+ic*(eiff-1.d0)**xi)+etar
      endif
c
      return
      end
c 
      subroutine cuntzepara045(idv,uffr,etae,etag,etaneg)
c----------------------------------------------------------------------
c     Degradation for Cuntze Model
c----------------------------------------------------------------------
      implicit none
      integer idv
      real*8 uffr,etae,etag,etaneg                 !etaneg bezieht sich nur auf den Zugbereich
      real*8 facte, factg, factneg
      real*8 a,b
c
c.... set constant values
      facte = 0.01                                 ! red.Faktor für E-Modul im Zugbereich
      factg = 0.50                                 ! red.Faktor für G-Modul im Zugbereich    0.05
      factneg = 0.50                               ! red.Faktor für G-Modul im Druckbereich  0.2
      a = 5                                      ! x-Wert für Geradengleichung  2.5
      b = (1+a)/2
c
      if(idv.eq.1) then
        if(a.eq.1.d0) then
          etae = facte
          etag = factg
          etaneg = factneg
        else
          etae = (facte-1.d0)/(a-1.d0)*(uffr-1.d0)+1.d0
          etag = (factg-1.d0)/(a-1.d0)*(uffr-1.d0)+1.d0
          if(etae.lt.facte) etae = facte
          if(etag.lt.factg) etag = factg
          etaneg = (etag-1)*(factneg-1)/(factg-1)+1
        endif
      elseif(idv.eq.2) then
        if(a.eq.1.d0) then
          etae = facte
          etag = factg
          etaneg = factneg
        else
          etag = (factneg-1.d0)/(a-1.d0)*(uffr-1.d0)+1.d0
          if(etag.lt.factneg) etag = factneg
        endif
      endif
c
      return 
      end
c
c
c
      subroutine hdam4515(d,h2,hdam)
c------------------------------------------------------------------
c... calculate input hdam from h2 values
c------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),h2(*),hdam(6)
c      
      h=1
c
      do i=1,6
         hdam(i)=0.d0
      enddo
c
c...  read input values
c
      iftyp = d(7)
      idtyp = d(23)
c...  Tsai-Wu criteria
      if (iftyp.eq.1) then
       if(h2(1).eq.1)  hdam(1) = hdam(1)+h
c
      elseif (iftyp.eq.2) then
c...  extended HASHIN criteria
c...  assign history variable
         idam = h2(1)
         idvt = h2(2)
         idvc = h2(3) 
c
       if(idam.eq.1.and.idvt.eq.1) hdam(1) = hdam(1) + h
       if(idam.eq.1.and.idvc.eq.1) hdam(2) = hdam(2) + h
       if(idam.eq.2) hdam(3) = hdam(3) + h
       if(idam.eq.3) hdam(4) = hdam(4) + h
       if(idam.eq.4) hdam(5) = hdam(5) + h   
c
c... Puck criteria
      elseif (iftyp.eq.3) then
c.... write history values
       idam = h2(1) 
       idva = h2(2)
       idvb = h2(3)
       idvc = h2(4)
       phic = h2(6)
       if(idam.eq.1) hdam(1) = hdam(1) + h      ! Matrix cracking
       if(idva.eq.1) hdam(2) = hdam(2) + h      ! Matrix Mode A
       if(idvb.eq.1) hdam(3) = hdam(3) + h      ! Matrix Mode B
       if(idvc.eq.1) hdam(4) = hdam(4) + h      ! Matrix Mode C
       if(idam.eq.2) hdam(5) = hdam(5) + h      ! Fibre fracture
                     hdam(6) = phic             ! Fracture angle Mode C
c
c... Cuntze criterion
      elseif (iftyp.eq.4) then
       idam = h2(1)
       iff1 = h2(2)
       iff2 = h2(3)
       iff3 = h2(4)
c
c.... write hdam for plot
       if(idam.eq.1) hdam(1) = hdam(1) + h      ! Matrix cracking
       if(iff1.eq.1) hdam(2) = hdam(2) + h      ! Matrix IFF1
       if(iff2.eq.2) hdam(3) = hdam(3) + h      ! Matrix IFF2a
       if(iff2.eq.1) hdam(4) = hdam(4) + h      ! Matrix IFF2b
       if(iff3.eq.1) hdam(5) = hdam(5) + h      ! Matrix IFF3
       if(idam.eq.2) hdam(6) = hdam(6) + h      ! Fiber fracture
      endif
c
      return
      end
c
c....End of mate15    