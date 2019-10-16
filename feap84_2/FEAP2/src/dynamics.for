c....FEAP-routines for dynamics
c-----------------------------------------------------------------------
c    dparam    set dynamic parameters
c    dsetci    set the integration constants
c    darray    form tangent matrix and residual modifications
c              for transient algorithm being used.
c    piacel    acceleration for form,acel
c    pisolv    solv for form,expl with/without damping
c    update    solution vectors at begin of time step/after iteration step
c    ploadd    calculate dynamic forces for equilibrium
c    resmod    modification of right hand side for lumped matrix
c    du_at_n12 calculate  uq(nneq) at n+1/2
c-----------------------------------------------------------------------
c
c
      subroutine dparam(ct,lct)
c-----------------------------------------------------------------------
c.... Purpose: set dynamic parameters for different algorithms
c
c     Input:
c            ct(3) -  input values for algorithm
c           lct    -  name of algorithm
c
c
c     nop = 1 Newmark
c           2 Backward-Euler
c           3 HHT-alpha
c           4 Generalized-Alpha
c           5 Energy Momentum Conserving
c           6 Implicit Composite Scheme - Bathe
c           7 Explizit Central Difference (beta=0) Taylor
c           8 Explizit Central Difference Abaqus
c
c     nrt = number of used vectors in array urate(nneq,nrt)
c     nrk,nrc,nrm vector adresses in urate for iteration
c-----------------------------------------------------------------------
c.... Algorithm  nrt  array urate
c.... Newmark -1   5  u(1)=v|n+1  u(2)=a|n+1  u(3)=v|n     u(4)=a|n  u(5)=a|n+a or v|n+a
c.... Back    -2   1  u(1)=v|n+1
c.... HHT     -3   5  u(1)=v|n+1  u(2)=a|n+1  u(3)=v|n     u(4)=a|n  u(5)=a|n+a or v|n+a
c.... G-alpha -4   5  u(1)=v|n+1  u(2)=a|n+1  u(3)=v|n     u(4)=a|n  u(5)=a|n+a or v|n+a
c.... Energy  -5   3  u(1)=v|n+1  u(2)=a|n+1  u(3)=a|n+1/2
c.... Bathe   -6      u(1)=v|n+1  u(2)=a|n+1  u(3)=u|n     u(4)=v|n  u(5)=a|n
c.... Expl    -7   4  u(1)=v|n+1  u(2)=a|n+1  u(3)=v|n     u(4)=a|n
c.... Expa    -8   5  u(1)=v|n+1  u(2)=a|n+1  u(3)=v|n-1/2 u(4)=a|n  u(5)=v|n+1/2
c-----------------------------------------------------------------------
      USE ddata
      USE iofile
      implicit double precision(a-h,o-z)
      character*4 lct
      dimension ct(3)
      theta(4) = 0.d0 ! for nop.ne.4

c     set values for algorithms

      if(lct.eq.'    '.or. lct.eq.'newm') then
c-----------------------------------------------------------------------
c....   newmark-beta method, ct(1) = beta  ;  ct(2) = gamma
c....   'nrk'        solution vector
c....   'nrc'        velocity vector v|n+1
c....   'nrm'    acceleration vector a|n+1
c....   'nrm+1'  =                   v|n
c....   'nrm+2'  =                   a|n,
c....   'nrm+3'  = help vector       => nrt=5
          nop = 1
          nrk = 0
          nrc = 1
          nrm = 2
          nrt = 5
          if(ct(1).eq.0.0d0) ct(1) = 0.25d0
          if(ct(2).eq.0.0d0) ct(2) = 0.5d0
                       write(iow,2001) ct(1),ct(2)
          if(ior.lt.0) write(*  ,2001) ct(1),ct(2)

      else if(lct.eq.'back') then
c-----------------------------------------------------------------------
c....   backward euler for first order ode (e.g. heat transfer)
c....   'nrk' solution vector
c....   'nrc' not used
c....   'nrm' first rate term
          nop = 2
          nrk = 0
          nrc = 0
          nrm = 1
          nrt = 1
c....     no theta(i) values are required (i.e., theta(1) is 1.0!)
          ct(1) = 1.
                       write(iow,2002)
          if(ior.lt.0) write(*  ,2002)

      else if(lct.eq.'hht') then
c-----------------------------------------------------------------------
c....   HHT-alpha method, ct(1) = alpha with 0<alpha<1/3
c....   then ct(1) = beta; ct(2) = gamma; ct(3) = alpha
c....   'nrk'         solution vector
c....   'nrc'         velocity vector v|n+1
c....   'nrm'     acceleration vector a|n+1
c....   'nrm+1'  =                    v|n
c....   'nrm+2'  =                    a|n,
c....   'nrm+3'  = help vector        => nrt=5
          nop = 3
          nrk = 0
          nrc = 1
          nrm = 2
          nrt = 7
          alpha = ct(1)
          if(alpha.eq.0.d0) alpha=0.05d0
c...      check limits of alpha
          if(alpha.gt.1.d0/3.d0)  stop '0<alpha<1/3 SR DPARAM'
          if(alpha.lt.0.d0)       stop '0<alpha<1/3 SR DPARAM'
          beta   = 0.25d0*(1.d0+alpha)**2
          gamma  = 0.5d0+alpha
          ct(1)  = beta
          ct(2)  = gamma
          ct(3)  = alpha
                       write(iow,2003) ct(1),ct(2),ct(3)
          if(ior.lt.0) write(*  ,2003) ct(1),ct(2),ct(3)

      else if(lct.eq.'alph') then
c-----------------------------------------------------------------------
c....   generalized-alpha method,
c       limits
c       ct(1) = rho with 0.5<rho<1,
c               (rho<0.5->alpham<0,rho>1->alpham>1/2 and alphaf>1/2)
c       alpham<1/2,alphaf<1/2
c....   then ct(1) = beta; ct(2) = gamma; ct(3) = alpham; ct(4)=alphaf
c       1) input rho .ne. 0: alpham, alphaf from rho
c       2) input rho  =   0: alpham, alphaf from input
c....   'nrk'     solution vector
c....   'nrc'     velocity vector v|n+1
c....   'nrm' acceleration vector a|n+1
c....   'nrm+1'  =                v|n
c....   'nrm+2'  =                a|n,
c....   'nrm+3'  = help vector    => nrt=5
          nop = 4
          nrk = 0
          nrc = 1
          nrm = 2
          nrt = 7
          rho = ct(1)
          if(rho.eq.0.d0) then ! direct input of alpham,alphaf
c           (rho=alpham=alphaf=0->Newmark)
            alpham= ct(2)
            alphaf= ct(3)
          else if (rho.ne.0.d0) then ! calculate alpham,alphaf from rho
c...        check limits of rho
            if(rho.lt.0.5d0) stop '0.5<rho<1) SR DPARAM'
            if(rho.gt.1.0d0) stop '0.5<rho<1) SR DPARAM'
            alpham = (2.d0*rho-1.d0)/(rho+1.0d0)
            alphaf =       rho      /(rho+1.0d0)
          end if
          beta   = 0.25d0*(1.d0-alpham+alphaf)**2
          gamma  = 0.5d0-alpham+alphaf
          ct(1)  = beta
          ct(2)  = gamma
          ct(3)  = alpham
          theta(4) = alphaf
                       write(iow,2004) rho,ct(1),ct(2),ct(3),theta(4)
          if(ior.lt.0) write(*  ,2004) rho,ct(1),ct(2),ct(3),theta(4)

      else if(lct.eq.'ener') then
c-----------------------------------------------------------------------
c....   energy momentum conserving method
c....   'nrk'       solution vector
c....   'nrc'       velocity vector v|n+1
c....   'nrm'   acceleration vector a|n+1
c....   'nrm+1' acceleration vector a|n+1/2
          nop = 5
          nrk = 0
          nrc = 1
          nrm = 2
          nrt = 3
                       write(iow,2005)
          if(ior.lt.0) write(*  ,2005)

      else if(lct.eq.'bath') then
c-----------------------------------------------------------------------
c....   Bathe composite implicit time integration procedure
c....   'nrk'       solution vector
c....   'nrc'       velocity vector v|n+gamma*dt (first substep), v|n+1 (second substep)
c....   'nrm'   acceleration vector a|n+gamma*dt (first substep), a|n+1 (second substep)
c....   'nrm+1'     velocity vector u|n for backstep
c....   'nrm+2' acceleration vector v|n
c....   'nrm+3'     velocity vector a|n
          nop = 6
          nrk = 0
          nrc = 1
          nrm = 2
          nrt = 5
          if(ct(1).eq.0.0d0) ct(1) = 2.d0-sqrt(2.d0)
          ct(2)=1.3d0
                       write(iow,2006) ct(1)
          if(ior.lt.0) write(*  ,2006) ct(1)

      else if(lct.eq.'expl') then
c-----------------------------------------------------------------------
c....   explicit Central difference method (beta=0) Taylor
c....   'nrk'       solution vector
c....   'nrc'       velocity vector v|n+1
c....   'nrm'   acceleration vector a|n+1
c....   'nrm+1'     velocity vector v|n
c....   'nrm+2' acceleration vector a|n
          nop = 7
          nrk = 0
          nrc = 1
          nrm = 2
          nrt = 4
                       write(iow,2007)
          if(ior.lt.0) write(*  ,2007)

      else if(lct.eq.'expa') then
c-----------------------------------------------------------------------
c....   explicit Central difference method (beta=0) Abaqus
c....   'nrk'       solution vector
c....   'nrc'       velocity vector v|n+1
c....   'nrm'   acceleration vector a|n+1
c....   'nrm+1'     velocity vector v|n-1/2
c....   'nrm+2' acceleration vector a|n
c....   'nrm+3'     velocity vector v|n+1/2
          nop = 8
          nrk = 0
          nrc = 1
          nrm = 2
          nrt = 5
                       write(iow,2008)
          if(ior.lt.0) write(*  ,2008)
      else
        call drawmess('No dynamic algorithm found',1,0)
      end if
c.... transfer values to 'theta' in common /ddata/
      do 100 i = 1,3
          theta(i) = ct(i)
100   continue
c
      return
2001  format(' Newmark  parameters'/'  beta = ',f9.4,' gamma = ',f9.4)
2002  format(' Backward Euler for first order systems.')
2003  format(' HHT-alpha parameters'/'  beta = ',f9.4,' gamma = ',f9.4,
     1                          ' alpha = ',f9.4)
2004  format(' Generalized-Alpha  parameters'/ '  rho   = ',f9.4,
     1       '   beta = ',f9.4,'  gamma = ',f9.4/
     2       '                    alpham = ',f9.4,' alphaf = ',f9.4)
2005  format(' Energy momentum conserving method.')
2006  format(' Implicit composite time integration procedure., Bathe'/
     1       'gamma = ',f9.4)
2007  format(' Explicit central difference method (beta=0)., Taylor')
2008  format(' Explicit central difference method (beta=0)., Abaqus')
      end
c

      subroutine dsetci
c------------------------------------------------------------------------------------------------------------------------|
c.... set the integration constants for a dynamic analysis                                                               |
c-------------------------------------------------------------------------------------------------------------------------
c    |  Newmark(1)     |Back(2)|   HHT(3) = 1 !  | GenAlph(4) = 1 !|Energy(5)| Bathe(7)               | Expl(8) | Expa(9)|
c ------------------------------------------------------------------------------------------------------------------------
c c1 | 1/(beta*dt*dt)  |1.d0/dt| 1/(beta*dt*dt)  | 1/(beta*dt*dt)  |2/(dt*dt)|4/(gamma*dt)^2          |0.5d0*dt |        |
c    |                 |       |                 |                 |         |c2*c2                   |*dt      |        |
c ------------------------------------------------------------------------------------------------------------------------
c c2 | gamm/(dt*beta)  |       | gamm/(dt*beta)  | gamm/(dt*beta)  |2/dt     |2/(gamma*dt)            |0.5*dt   |        |
c    |                 |       |                 |                 |         |(2-gamma)/((1-gamma)*dt)|         |        |
c ------------------------------------------------------------------------------------------------------------------------
c c3 | 1-1/(2*beta)    |       | 1-1/(2*beta)    | 1-1/(2*beta)    |         |                        |         |        |
c    |                 |       |                 |                 |         |(1-gamma)/(gamma*dt)    |         |        |
c ------------------------------------------------------------------------------------------------------------------------
c c4 | 1-gamm/beta     |       | 1-gamm/beta     | 1-gamm/beta     |         |                        |         |        |
c    |                 |       |                 |                 |         |-1/((1-gamma)*gamma*dt) |         |        |
c ------------------------------------------------------------------------------------------------------------------------
c c5 |(1-gamm/(2*beta))|       |(1-gamm/(2*beta))|(1-gamm/(2*beta))|         |                        |         |        |
c    | *dt             |       | *dt             | *dt             |         |                        |         |        |
c-------------------------------------------------------------------------------------------------------------------------

      USE ddata
      USE tdata
      implicit double precision (a-h,o-z)

      if(nop.eq.1.or.nop.eq.3.or.nop.eq.4) then

c....   newmark-beta,
c       HHT-alpha and generalized-alpha parameters(beta,gamma different)
        beta = theta(1)
        gamm = theta(2)
        c1 = 1.d0/(beta*dt*dt)
        c2 = gamm/(dt*beta)
        c3 = 1.d0 - 1.d0/(beta+beta)
        c4 = 1.d0 - gamm/beta
        c5 = (1.d0 - gamm/(beta+beta))*dt

      else if(nop.eq.2) then

c....   backward euler constants
        c1 = 1.d0/dt

      else if(nop.eq.5) then

c....   energy momentum conserving parameters
        c1 = 2.d0/(dt*dt)
        c2 = 2.d0/dt

      else if(nop.eq.6) then

c.... Implicit Composite Scheme - Bathe

        gamm=theta(1)

        if(int(theta(2)).eq.2) then
            !dt equals gamma*dt
            c1 = 4/(dt*dt)
            c2 = 2/(dt)
            c3 =-4/dt

        else
            !dt equals (1-gamma)*dt
            c2 = (2.d0-gamm)/(dt)
            c1 = c2*c2
            c3 = (1.d0-gamm)*(1.d0-gamm)/(gamm*dt)
            c4 = -1.d0/(gamm*dt)


        endif


      else if(nop.eq.7) then

c....   explicit central difference parameters Taylor
        c1 = 0.5d0*dt*dt
        c2 = 0.5d0*dt

c     else if(nop.eq.8) then

c....   explicit central difference parameters Abaqus

      end if
      return
      end
c


      subroutine darray(dr,al,au,ad,xml,xmu,xmd,xdl,xdu,xdd,ur,jd,
     1                  nneq,aflg,mflg,fflg,uflg)
c-----------------------------------------------------------------------
c.... Purpose: form tangent matrix and residual modifications for
c              transient algorithms being used.
c
c      Arguments
c          dr(1)     - residual vector
c          al(1)     - lower part of tangent matrix (only if uflg true)
c          au(1)     - upper part of tangent matrix
c          ad(1)     - diagonal part of tangent matrix
c          xml(1)    - lower part of "mass" matrix (only if xxxx true)
c          xmu(1)    - upper part of "mass" matrix
c          xmd(1)    - diagonal part of "mass" matrix
c          xdl(1)    - lower part of "damping" matrix (only if xxx true)
c          xdu(1)    - upper part of "damping" matrix
c          xdd(1)    - diagonal part of "damping" matrix
c          ur(nneq,1)- rate terms defined at trans
c          jd(1)     - pointers to bottom of column/rows of al,au, xml
c          csrka(1)  - pointers to diagonals CSR
c          neq       - number of active equations
c          nneq      - total number of d.o.f. in problem
c          aflg      - form modified tangent matrix if true
c          mflg      - mass is consistent if true, else lumped
c                      and stored in xmd (xmu is not used)
c          fflg      - form modified force vector if true
c          uflg      - tangent array is unsymmetric if true
c          flgda     - damping matrix is consistent if true, else lumped
c
c-----------------------------------------------------------------------
c         Standard Solver 0 and 4-8
c           TANG; CMAS; LMAS; DAMP,CONS; DAMP,LUMP
c           UTAN; UMAS;       DAMP,UCON
c
c           comment: UMAS, DAMP,UCON possible only with UTAN
c
c         Solver 1,2, 3 sparse matrix
c           TANG; CMAS; LMAS; DAMP,CONS; DAMP,LUMP
c
c.....     ww 1/2005
c          ww IBS KIT Update 12/2014 
c 
c-----------------------------------------------------------------------
c      Required operations for darray
c
c      a.) Tangent modification requires
c
c               a  <--- ci * xm + ck * xd + cj * a
c
c          where 'a' includes 'ad' and 'au' (if 'uflg' true
c          also includes 'al'), 'xm' includes 'xmd' and 'xmu',
c          'xd' includes 'xdd' and 'xdu',and 'ci','ck' and 'cj'
c          are parameters defined by the algorithm.
c
c      b.) Residual modification requires
c
c               dr <---  dr - xm * ur(nrm) - xd * ur(nrc)
c
c          where xm includes 'xmd' and 'xmu', xd includes 'xdd' and 'xdu'
c          and 'nrm','nrc' are specified columns of the ur array.
c
c-----------------------------------------------------------------------
c
c      Example algorithm:  Newmark
c
c               a  <--- a  + c1 * xm + c2 * xd
c
c               dr <--- dr - xm * ur(nrm) - xd * ur(nrc)
c
c-----------------------------------------------------------------------

      USE cdata
      USE damp1
      USE ddata
      USE iscsr
      USE ndata
      USE soltyp
      USE tdata
      implicit double precision (a-h,o-z)
c
      logical aflg,fflg,mflg,uflg
      dimension dr(*),al(*),au(*),ad(*),xml(*),xmu(*),xmd(*),
     +          xdl(*),xdu(*),xdd(*),ur(nneq,*),jd(*)
c
      if(nop.eq.7) goto 107
      if(nop.eq.8) goto 108

c-----------------------------------------------------------------------
c
c.... modify tangent matrix if 'aflg' is true
c
c-----------------------------------------------------------------------

      if (aflg) then
c....   set the integration parameters for current time step
        call dsetci
        if( nop.le.6) then
c....     newmark:                       c1*M +c2*C + K_T
c....     back euler:                    c1*M       + K_T
c....     HHT-alpha:                  c6*c1*M +c2*C + K_T
c....     generalized-alpha:          c6*c1*M +c2*C + K_T
c....     energy momentum conserving:    c1*M + K_T (no damping!)
c....     Implicit Composite Scheme:     c1*M +c2*C + K_T
          c6 = 1.d0
          if(nop.eq.3) c6 =            1.0d0 / (1.0d0-theta(3))
          if(nop.eq.4) c6 = (1.0d0-theta(3)) / (1.0d0-theta(4))


c-----------------------------------------------------------------------
c
c     modify for a mass term
c
c-----------------------------------------------------------------------


          if(istyp.eq.0) then ! standard solver
c-----------------------------------------------------------------------
c....       modify for mass terms:  diagonal terms (lumped mass)
            do n = 1,neq
              ad(n) = ad(n) + c6*c1*xmd(n)
            end do
c....       mass terms: if consistent mass modify remainder
            if (mflg) then
              do n = 1,jd(neq)
                au(n) = au(n) + c6*c1*xmu(n)
c....           for UTAN form lower part of mass too
                if(uflg) al(n) = al(n) + c6*c1*xml(n)
              end do
            end if

          else if(istyp.eq.1.or.istyp.eq.2) then ! sparse matrix solver
c-----------------------------------------------------------------------
            if (mflg) then
c....         consistent mass terms
              do n = 1,jd(neq+1)-1
                ad(n) = ad(n) + c6*c1*xmd(n)
              end do
c....         unsymmetric part not allowed
              if(uflg) stop
     +        ' Sparse Matrix Solv: UMAS not allowed in TRANS,SR darray'
            else
c....         lumped mass terms
              do n = 1,neq
                ad(jd(n)) = ad(jd(n)) + c6*c1*xmd(n)
              end do
            end if

          else if(istyp.ge.3.and.istyp.le.8) then !S-LU/Pardiso/PBCG/PGMRES
c-----------------------------------------------------------------------
            if (mflg) then
c....         consistent mass terms sym+usym
              do n = 1,jd(neq+1)-1
                ad(n) = ad(n) + c6*c1*xmd(n)
              end do
            else
c....         lumped mass terms
              call dadd_csr(neq,ad,xmd,c6*c1,csrka)
            end if
c-----------------------------------------------------------------------

          end if

        end if


c-----------------------------------------------------------------------
c
c....   modify for a damping term
c
c-----------------------------------------------------------------------

        if(nc.gt.1) then

          if(istyp.eq.0) then ! standard solver
c-----------------------------------------------------------------------
c....       modify damping terms:  diagonal terms (lumped damping)
c           call colred(xdd,-c2,neq, ad)
            call daxpty(neq,c2,xdd,ad)
c....       damping terms: if consistent damping modify remainder
            if (flgda) then
c             call colred(xdu,-c2,jd(neq),au)
              call daxpty(jd(neq),c2,xdu,au)
c....         for UTAN form lower part of damping too
c             if(uflg) call colred(xdl,-c2,jd(neq), al)
              if(uflg) call daxpty(jd(neq),c2,xdl,al)
            end if

          else if(istyp.eq.1.or.istyp.eq.2) then ! sparse matrix solver
c-----------------------------------------------------------------------
            if (flgda) then
c....         consistent damping terms
              do n = 1,jd(neq+1)-1
                ad(n) = ad(n) + c2*xdd(n)
              end do
            else
c....         diagonal damping terms
              do n = 1,neq
                ad(jd(n)) = ad(jd(n)) + c2*xdd(n)
              end do
            end if

          else if(istyp.ge.3.and.istyp.le.8) then !S-LU/Pardiso/PBCG/PGMRES
c-----------------------------------------------------------------------
            if (flgda) then
c....         consistent damping terms
              do n = 1,jd(neq+1)-1
                ad(n) = ad(n) + c2*xdd(n)
              end do
            else
c....         diagonal damping terms
              call dadd_csr(neq,ad,xdd,c2,csrka)

            end if
c-----------------------------------------------------------------------

          end if ! istyp

        end if ! nc
      end if ! aflg

c-----------------------------------------------------------------------
c
c.... modify residual vector if 'fflg' is true
c
c-----------------------------------------------------------------------

      if (fflg) then

                 imas=2
        if(mflg) imas=1

        if(nop.lt.3.or.nop.eq.6) then
c-----------------------------------------------------------------------
c....     newmark, backward euler, Bathe: modify for  mass term R = R - M*a|n+1
          call promul(xml,xmu,xmd,ur(1,nrm),dr,jd,neq,imas,2)

        else if(nop.eq.3) then
c-----------------------------------------------------------------------
c....     HHT-alpha: modify for  mass term R = R - M*c6*a|n+1
c....     set acceleration a|n+a
          c6 = 1.0d0 / (1.0d0-theta(3))
          c7 = theta(3)*c6
          do n = 1,nneq
            ur(n,6)=c7*dr(n)          ! storing Fint|t
            dr(n)=dr(n)+ur(n,7)       ! adding Fint|t-1
            ur(n,5) = c6*ur(n,2)
          end do
          call promul(xml,xmu,xmd,ur(1,5),dr,jd,neq,imas,2)

        else if(nop.eq.4) then
c-----------------------------------------------------------------------
c....     generalized alpha: modify for mass term R = R - M*(c6*a|n+1 +c7*a|n)
c....     set acceleration a|n+a
          c6 = (1.0d0-theta(3)) / (1.0d0-theta(4))
          c7 =        theta(3)  / (1.0d0-theta(4))

          c8 =        theta(4)  / (1.0d0-theta(4))

          do n = 1,nneq
            ur(n,6)=c8*dr(n)          ! storing Fint|t
            dr(n)=dr(n)+ur(n,7)       ! adding Fint|t-1
            ur(n,5) = c6*ur(n,2) + c7*ur(n,4)
          end do
          call promul(xml,xmu,xmd,ur(1,5),dr,jd,neq,imas,2)

        else if(nop.eq.5) then
c-----------------------------------------------------------------------
c....     energy momentum conserving
c....     modify for  mass term R = R - M*a|n+1/2
          call promul(xml,xmu,xmd,ur(1,nrm+1),dr,jd,neq,imas,2)
        end if
c

c-----------------------------------------------------------------------

        if(nc.gt.1) then

c-----------------------------------------------------------------------

                     imas=2
          if (flgda) imas=1

          if(nop.lt.3.or.nop.eq.6) then
c-----------------------------------------------------------------------
c....       modify for  damping term R = R - C*v|n+1
            call promul(xdl,xdu,xdd,ur(1,nrc),dr,jd,neq,imas,2)

          else if(nop.eq.3) then
c-----------------------------------------------------------------------
c....       modify for  damping term R = R - C*(v|n+1+c6*v|n)
c....       set velocity
            c6 = theta(3)/(1.0d0-theta(3))
            do n = 1,nneq
              ur(n,5) = ur(n,1) + c6*ur(n,3)
            end do
            call promul(xdl,xdu,xdd,ur(1,5),dr,jd,neq,imas,2)

          else if(nop.eq.4) then
c-----------------------------------------------------------------------
c....       modify for  damping term R = R - C*(v|n+1+c6*v|n)
c....       set velocity
            c6 = theta(4)/(1.0d0-theta(4))
            do n = 1,nneq
              ur(n,5) = ur(n,1) + c6*ur(n,3)
            end do
            call promul(xdl,xdu,xdd,ur(1,5),dr,jd,neq,imas,2)

c         else if(nop.eq.5) then
c-----------------------------------------------------------------------
c           no damping possible
          end if

        end if
      end if
      return


107   continue

c.... expl. algorithm Taylor:  modify right hand side for damping term
c-----------------------------------------------------------------------
      if(nc.gt.1) then
c     for all solvers Standard/SuperLU/Pardiso/SM/PCG/PGMRES
        do n = 1,neq
          dr(n) = dr(n) - xdd(n)*ur(n,nrc) ! ok for istyp 0,1,2
        end do
      end if
      return

108   continue

c.... expl. algorithm Abaqus: modify right hand side for damping term???
c-----------------------------------------------------------------------

      return
      end
c

      subroutine piacel(xml,dr,a,neq)
c-----------------------------------------------------------------------
c.... Purpose: compute starting acceleration
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 xml(*),dr(*),a(*)
      do 100 n = 1,neq
        if(xml(n).ne.0.0d0) then
           a(n) = dr(n)/xml(n)
        else
           a(n) = 0.0d0
        end if
100   continue
      return
      end
c

      subroutine pisolv(xmd,xdd,nc,dr,a,neq)
c-----------------------------------------------------------------------
c
c.... Purpose: solve for form,expl
c.... EXPL(RLT)   nop=6: a_i=(F_i^ext-F_i^int-C_ii*v_i^0)/(M_ii+0.5*dt*C_ii)
c.... EXPA(ABAQUS)nop=7: a_i=(F_i^ext-F_i^int)/M_ii
c
c-----------------------------------------------------------------------
      USE tdata
      implicit real*8 (a-h,o-z)
      real*8 xmd(*),xdd(*),dr(*),a(*)

      if(nc.gt.1) then ! with damping, only nop=7

        do  n = 1,neq
          xn = xmd(n) + c2*xdd(n)
          if(xn.ne.0.0d0) then
            a(n) = dr(n)/xn
          else
            a(n) = 0.0d0
          end if
        end do

      else ! without damping

        do  n = 1,neq
          xn = xmd(n)
          if(xn.ne.0.0d0) then
            a(n) = dr(n)/xn
          else
            a(n) = 0.0d0
          end if
c        write(*,*) 'R_i   ',a(n)
        end do
      end if
      return
      end


c
      subroutine update(id,f0,f,u,urate,du,nneq,nqq,fdyn,pfr,isw)
c-----------------------------------------------------------------------
c     Purpose: update solution vectors
c              to begin of a time step/ after iteration step
c
c        isw = 1 at begin of time      step
c        isw = 2 after       iteration step
c        isw = 3 to make a   back      step
c        isw = 4 set start increment at begin of time step
c
c     urate = trans
c
c.... Algorithm  array urate
c.... Algorithm  nrt array urate
c.... Newmark -1   5  u(1)=v|n+1  u(2)=a|n+1  u(3)=v|n     u(4)=a|n u(5)=a|n+a or v|n+a
c.... Back    -2   1  u(1)=v|n+1
c.... HHT     -3   5  u(1)=v|n+1  u(2)=a|n+1  u(3)=v|n     u(4)=a|n u(5)=a|n+a or v|n+a
c.... G-alpha -4   5  u(1)=v|n+1  u(2)=a|n+1  u(3)=v|n     u(4)=a|n u(5)=a|n+a or v|n+a
c.... Energy  -5   3  u(1)=v|n+1  u(2)=a|n+1  u(3)=a|n+1/2
c.... Bath    -6   5  u(1)=v|n+gt u(2)=a|n+gt u(3)=u|n     u(4)=v|n u(5)=a|n        gt=gamma*dt
c.... Expl    -7   4  u(1)=v|n+1  u(2)=a|n+1  u(3)=v|n     u(4)=a|n
c.... Expl    -8   5  u(1)=v|n+1  u(2)=a|n+1  u(3)=v|n-1/2 u(4)=a|n u(5)=v|n+1/2
c
c-----------------------------------------------------------------------
c
c..... type declaration for variables
      USE back1
      USE cdata
      USE ddata
      USE edgdat
      USE hdatam
      USE iofile
      USE psize
      USE prlod
      USE tdata
      logical fdyn,pfr
      integer j, m, n, nneq,nneq2,nqq,isw
      real*8  c6,vnorm,anorm
      real*8  ddot, det, ub,ur1,ur2

c.... type declaration for arrays
      integer id(*)
      real*8  f0(*),f(*),u(*),urate(nneq,*),du(*)

c.... intrinsic declarations
      intrinsic sqrt
c-----------------------------------------------------------------------
c
      if(isw.eq.1) then !  at time step
c
c-----------------------------------------------------------------------
        if(pfr) then

c....     compute and output norms in v and a
          vnorm = sqrt(ddot(nneq,urate(1,1),1,urate(1,1),1))
          anorm = 0.0d0
          if(nop.eq.1 .or. nop.ge.3 ) then
            anorm = sqrt(ddot(nneq,urate(1,2),1,urate(1,2),1))
          end if
                       write(iow,2000) vnorm,anorm
          if(ior.lt.0) write(*  ,2000) vnorm,anorm
        end if

        if(nop.eq.1) then
c-----------------------------------------------------------------------
c....     newmark update:
c....     save v|n and a|n for time stepping [auto] only
          do n = 1,nneq
            urate(n,3) = urate(n,1)
            urate(n,4) = urate(n,2)
          end do

c....     newmark update: set part from t|n:  v|n+1_0, a|n+1_0
          c6 = dt*c1
          do n = 1,nneq
            ur1        = - c6*urate(n,1) + c3*urate(n,2) !v = urate(n,1)
            urate(n,1) =   c4*urate(n,1) + c5*urate(n,2) !a = urate(n,2)
            urate(n,2) = ur1
          end do

        else if(nop.eq.2) then
c-----------------------------------------------------------------------
c....     backward euler update
          do n = 1,nneq
            urate(n,1) = 0.0d0
          end do

        else if(nop.eq.3) then
c-----------------------------------------------------------------------
c....     HHT-alpha update
c....     save v|n and a|n
          do n = 1,nneq
            urate(n,3) = urate(n,1)
            urate(n,4) = urate(n,2)
          end do
c....     set part from t|n:  v|n+1_0, a|n+1_0 see Newmark
          c6 = dt*c1
          do n = 1,nneq
            ur1        = - c6*urate(n,1) + c3*urate(n,2) !v = urate(n,1)
            urate(n,1) =   c4*urate(n,1) + c5*urate(n,2) !a = urate(n,2)
            urate(n,2) = ur1
            urate(n,7) = urate(n,6)
          end do

        else if(nop.eq.4) then
c-----------------------------------------------------------------------
c....     generalized-alpha update
c....     save v|n and a|n
          do n = 1,nneq
            urate(n,3) = urate(n,1)
            urate(n,4) = urate(n,2)
          end do
c....     set part from t|n:  v|n+1_0, a|n+1_0 see Newmark
          c6 = dt*c1
          do n = 1,nneq
            ur1        = - c6*urate(n,1) + c3*urate(n,2) !v = urate(n,1)
            urate(n,1) =   c4*urate(n,1) + c5*urate(n,2) !a = urate(n,2)
            urate(n,2) = ur1
            urate(n,7) = urate(n,6)
          end do

        else if(nop.eq.5) then
c-----------------------------------------------------------------------
c....     energy momentum conserving update
c....     acceleration a|n+1/2_0 = -c2*v|n
c....     acceleration a|n+1_0   = -c2*v|n - a|n
c....     velocity     v|n+1_0   = -v|n
          do n = 1,nneq
            urate(n,3) =      -c2*urate(n,1)
            urate(n,2) = -2.d0*c2*urate(n,1) - urate(n,2)
            urate(n,1) =         -urate(n,1)
          end do

        else if(nop.eq.6) then
c-----------------------------------------------------------------------
c....     composite procedure update - Bathe
c....     newmark step update: set part from t|n:  v|n+1_0, a|n+1_0
          if(int(theta(2)).eq.2)then
           do n = 1,nneq
            urate(n,4) = urate(n,1)                        !storing v|t = urate(n,4)
            urate(n,5) = urate(n,2)                        !storing a|t = urate(n,5)
            ur1        = c3*urate(n,1) - urate(n,2)
            urate(n,1) = - urate(n,1)                      !v|t+gamma*dt = urate(n,1)
            urate(n,2) = ur1                               !a|t+gamma*dt = urate(n,2)
            j=id(n)
            if(j.gt.0)then                                 !storing u|t = urate(n,3)
             urate(j,3)=u(n)
            else
             urate(neq-j,3)=u(n)
            end if
           end do

          else if(int(theta(2)).eq.1) then

           do n=1,nneq
            j=id(n)
            if(j.gt.0)then
             ur1=c3*urate(j,3)+(c2+c4)*u(n)
             ur2=c3*urate(j,4)+c4*urate(j,1)+c2*ur1
             urate(j,1)=ur1
             urate(j,2)=ur2
            else
             ur1=c3*urate(neq-j,3)+(c2+c4)*u(n)
             ur2=c3*urate(neq-j,4)+c4*urate(neq-j,1)+c2*ur1
             urate(neq-j,1)=ur1
             urate(neq-j,2)=ur2
            end if
           end do
          end if
        else if(nop.eq.7) then
c-----------------------------------------------------------------------
c....     explicit central difference update Taylor
c....     save v|n and a|n
          do n = 1,nneq
            urate(n,3) = urate(n,1) ! save v_n
            urate(n,4) = urate(n,2) ! save a_n
          end do
c....     Compute displacements  u|n+1
          do n = 1,nneq
            j = id(n)
            if(j.gt.0) then
c             u|n+1 = u|n + v|n+1/2*dt with v|n+1/2 = v|n + 1/2a|n*dt
              ub         = dt*urate(j,3) + c1*urate(j,4)
              u(n)       = u(n) + ub
            else
c             Compute values from forced inputs for fixed dof
              u(n)       = f0(n) + f(n)*prop
c             WARNING: Boundary velocity and accelerations are NOT computed.
c             Do not use consistent mass or problems with damping
c             with specified boundary motions!
            end if
          end do
c....     compute initial velocity v|n+1_0
          do n = 1,nneq
            urate(n,1) = urate(n,1)+ c2*urate(n,2)
c           compute v_n+1_0 = v|n+0.5*a|n*dt, 2nd part under isw=2
          end do

        else if(nop.eq.8) then
c-----------------------------------------------------------------------
c....     explicit central difference update Abaqus
c....     save v|n-1/2 and a|n
          do n = 1,nneq
            urate(n,3) = urate(n,5) ! save v_n-1/2
            urate(n,4) = urate(n,2) ! save a_n
          end do
c....     compute velocity v|n+1/2
          do n = 1,nneq
            urate(n,5) = urate(n,3)+ dt*urate(n,4)
c           compute v_n+1/2 = v|n-1/2 + a|n * 0.5(dt_n+dt_n+1)
          end do
c....     Compute initial displacements  u|n+1
          do n = 1,nneq
            j = id(n)
            if(j.gt.0) then
c             u|n+1 = u|n + v|n+1/2*dt_n+1
              ub         = dt*urate(j,5)
              u(n)       = u(n) + ub
            else
c             Compute values from forced inputs for fixed dof
              u(n)       = f0(n) + f(n)*prop
c             WARNING: Boundary velocity and accelerations are NOT computed.
c             Do not use consistent mass or problems with damping
c             with specified boundary motions!
            end if
          end do
        end if

c-----------------------------------------------------------------------
c
      else if(isw.eq.2) then ! at iteration step
c
c-----------------------------------------------------------------------

c....   update displacements and its increments within the time step
c....   update the velocities etc for transient analysis
        nneq2 = nneq + nneq

        if(nop.eq.7.or.nop.eq.8) goto 201 ! no update necessary
c-----------------------------------------------------------------------

c....   update displacements and its increments within the time step
        do n = 1,nneq
          j = id(n)
          if (j.gt.0) then ! no update of displacements necessary
c....       for the active degrees-of-freedom compute values from
c....       solution where 'du(j)' is and increment of 'u' for
c....       active dof 'j'.
            u(n)       = u(n)      + du(j)
            u(n+nneq)  = u(n+nneq) + du(j)
            u(n+nneq2) =             du(j)
          else
c....       for fixed degrees-of-freedom compute values from forced inputs
            ub         = f0(n) + f(n)*prop
            du(neq-j)  = ub - u(n)
            u(n+nneq2) = ub - u(n)
            u(n+nneq)  = u(n+nneq) + ub - u(n)
            u(n)       = ub
          end if
        end do
201   continue
        if(nde*ned.gt.0) then
          call edgupd(du,edge1,edge3,edge4,edge4(ne4),nde,neq,numnp)
        end if

c....   for time dependent solutions update the rate terms
        if(fdyn) then

          if(nop.eq.1.or.nop.eq.3.or.nop.eq.4.or.nop.eq.6) then
c-----------------------------------------------------------------------
c....       newmark-beta update for velocity and acceleration
c....       same for HHT-alpha and Generalized alpha and Bathe scheme
            do n = 1,nneq
              urate(n,1) = urate(n,1) + c2*du(n)
              urate(n,2) = urate(n,2) + c1*du(n)
            end do

          else if(nop.eq.2) then
c-----------------------------------------------------------------------
c....       update the backward euler solution
            do n = 1,nneq
              urate(n,1) = urate(n,1) + c1*du(n)
            end do

          else if(nop.eq.5) then
c-----------------------------------------------------------------------
c....       energy momentum conserving algorithm:
c           update velocity v|n+1,acceleration a|n+1,acceleration a|n+1/2
            do n = 1,nneq
              urate(n,1) = urate(n,1) + c2*du(n)
              urate(n,2) = urate(n,2) + c1*du(n)*2.d0
              urate(n,3) = urate(n,3) + c1*du(n)
            end do

          else if(nop.eq.7) then
c-----------------------------------------------------------------------
c....       explicit cent-diff Taylor: update for velocity, set accel.
            do n = 1,nneq
              urate(n,1) = urate(n,1)+c2*du(n)!only correct without iterat.
              urate(n,2) =               du(n)
            end do
c...        displacement is set under isw=1

          else if(nop.eq.8) then
c-----------------------------------------------------------------------
c....       explicit cent-diff Abaqus: set accel., set output velocity v_n+1
            do n = 1,nneq
              urate(n,2) = du(n)
              urate(n,1) = du(n)*0.5d0*dt + urate(n,5)
            end do
c...        displacement|n+1 and velocity|n+1/2 are set under isw=1
          end if
        end if
c....   set update flag for history variables
        hflgu  = .true.
        h3flgu = .true.
c
c-----------------------------------------------------------------------
c
      else if(isw.eq.3) then !backup solution vectors to reinitiate a step
c
c-----------------------------------------------------------------------
c
        nneq2 = nneq + nneq

        if(nop.eq.1) then
c-----------------------------------------------------------------------
c....     Newmark
c....     correct but not active
c         c6 = dt*c1
c         det= c3*c4 + c5*c6
c....     values are stored in urate(3),urate(4)

c       else if(nop.eq.3) then
c-----------------------------------------------------------------------
c....     HHT-alpha
c....     values are stored in urate(3),urate(4)

c       else if(nop.eq.4) then
c-----------------------------------------------------------------------
c....     Generalized-alpha
c....     values are stored in urate(3),urate(4)

c       else if(nop.eq.5) then
c-----------------------------------------------------------------------
c....     Energy momentum conserving
c....     no further constants
        else if(nop.eq.6) then
c-----------------------------------------------------------------------
c....   Implicit composite scheme - Bathe
        do n=1,nneq
        j=id(n)
        du(n)=0
        if(j.gt.0)then
            u(n)=urate(j,3)     !displacements are stored in urate(3)
            u(n+nneq)  = 0.
            u(n+nneq2) = 0.
        else
            u(n)=urate(neq-j,3)
            u(n+nneq)  = 0.
            u(n+nneq2) = 0.
        end if
            urate(n,1)=urate(n,4) !velocities are stored in urate(4)
            urate(n,2)=urate(n,5) !accelerations are stored in urate(5)
        end do
        else if(nop.eq.7) then
c-----------------------------------------------------------------------
c....     Explicit Taylor
          goto 307

        else if(nop.eq.8) then
c-----------------------------------------------------------------------
c....     Explicit Abaqus
          goto 308

        end if
c-----------------------------------------------------------------------
c       displacements      
        do 300 n = 1,nneq
cww       ub = u(n+nneq)
          j  = id(n)
          if (j.gt.0) then
c....       for the active degrees-of-freedom compute values from solution
c....       where 'du(j)' is and increment of 'u' for active dof 'j'.
cww         u(n)       = u(n) - ub
            u(n)       = ustore(n)  ! this is open for all dynamic parts
            u(n+nneq)  = 0.
            u(n+nneq2) = 0.
          else
c....       for fixed degrees-of-freedom compute values from forced inputs
            du(neq-j)  = 0.
            u(n+nneq2) = 0.
            u(n+nneq)  = 0.
cww         u(n)       = f0(n) + f(n)*prop ! prop is set back under BACK
            u(n)       = ustore(n)
          end if
c

          if(fdyn)then
            if (j.gt.0) then

              if(nop.eq.1) then
c-----------------------------------------------------------------------
c....           newmark-beta updates
c....           correct but not active
c               ur1 = urate(j,1) - c2*ub
c               ur2 = urate(j,2) - c1*ub
c               urate(j,1) = (c3*ur1 - c5*ur2)/det
c               urate(j,2) = (c6*ur1 + c4*ur2)/det
                urate(j,1) = urate(j,3) ! vel. v|n are saved in urate(3),see isw=1
                urate(j,2) = urate(j,4) ! acc. a|n are saved in urate(4),see isw=1

              else if(nop.eq.3.or.nop.eq.4) then
c-----------------------------------------------------------------------
c....           HHT-alpha or generalized-alpha updates
                urate(j,1) = urate(j,3) ! vel. v|n are saved in urate(3),see isw=1
                urate(j,2) = urate(j,4) ! acc. a|n are saved in urate(4),see isw=1

              else if(nop.eq.5) then
c-----------------------------------------------------------------------
c....           energy momentum conserving
                ur1 =       c2*ub -         urate(j,1)
                ur2 = -2.d0*c1*ub + 2.d0*c2*urate(j,1) - urate(j,2)
                urate(j,1) = ur1
                urate(j,2) = ur2

c             else if(nop.eq.7) then
c-----------------------------------------------------------------------
c....           explicit central difference updates Taylor: nothing to do

c             else if(nop.eq.8) then
c-----------------------------------------------------------------------
c....           explicit central difference updates Abaqus: nothing to do

              end if
            end if
          end if
300     continue
        goto 302
307     continue ! step back for central differences Taylor
c....   displacements
        do n = 1,nneq
          j = id(n)
          if(j.gt.0) then
            ub         = -dt*urate(j,3) - c1*urate(j,4)
            u(n)       = u(n) + ub
          else
            u(n)       = f0(n) + f(n)*prop ! prop is set back under BACK
          end if
        end do
c....   velocities,accelerations
        do n = 1,nneq
          urate(n,1) = urate(n,3)
          urate(n,2) = urate(n,4)
        end do
        goto 302
308     continue ! step back for central differences Abaqus
c....   displacements
        do n = 1,nneq
          j = id(n)
          if(j.gt.0) then
            ub         = -dt*urate(j,5)
            u(n)       = u(n) + ub
          else
            u(n)       = f0(n) + f(n)*prop ! prop is set back under BACK
          end if
        end do
c....   velocities,accelerations
        do n = 1,nneq
          urate(n,5) = urate(n,3)
          urate(n,2) = urate(n,4)
        end do
302     continue


      else if(isw.eq.4) then !  start increment time step
c
c-----------------------------------------------------------------------
        ddt=dt/dto

        do n = 1,nneq
          j = id(n)
          if (j.gt.0) then !
c....       for the active degrees-of-freedom compute values from
c....       solution where 'u(n+nneq)' is  (u_i-1 - u_i-2)
c....       active dof 'j'.

            u(n) = u(n) + ddt*u(n+nneq)

          else
c....       for fixed degrees-of-freedom compute values from forced inputs
            ub         = f0(n) + f(n)*prop
            u(n)       = ub
          end if
        end do

      end if
c
cww2000  format('   N o r m s   f o r   D y n a m i c s'/
cww     1   10x,'Velocity:',e13.5,' Acceleration:',e13.5)
2000  format('   Velocity norm = ',1pe15.7,
     1     ' Acceleration norm = ',1pe15.7)
c
      end
c
      subroutine ploadd(du,dr,xml,xmu,xmd,xdl,xdu,xdd,ur,jd,neq,nneq,
     + mflg,ityp)
c-----------------------------------------------------------------------
c
c.... Purpose: calculate dynamic forces in residual
c          ityp      1 at n
c                    2 at n+1/2
c########## open explicit??!!
c
c      Inputs:
c          ityp      - see above
c          du(*)     - displacement vector for ityp=2
c          dr(*)     - residual vector
c          xml(*)    - lower part of "mass" matrix
c          xmu(*)    - upper part of "mass" matrix
c          xmd(*)    - diagonal part of "mass" matrix
c          xdl(*)    - lower part of "damping" matrix
c          xdu(*)    - upper part of "damping" matrix
c          xdd(*)    - diagonal part of "damping" matrix
c          ur(nneq,1)- rate terms
c          jd(*)     - pointers to bottom of column/rows of al,au, xml
c          neq       - number of active equations
c          nneq      - total number of d.o.f. in problem
c          mflg      - mass is consistent if true, else lumped
c          flgda     - damping matrix is consistent if true, else lumped
c
C          nop       - 1 Newmark, 2 Backward-Euler, 3 HHT-alpha
c                    - 4 Generalized-Alpha, 5 Energy Momentum Conserving
c                    - 6 Explizit Central Difference (beta=0) Taylor
c                    - 7 Explizit Central Difference Abaqus
c
c      Outputs:
c          dr(1)     - modified residual vector

c
c-----------------------------------------------------------------------
c          Residual modification requires F=F_ext - F_int - M a - D v
c
c               dr <---  dr - xm * ur(nrm) - xd * ur(nrc)
c
c-----------------------------------------------------------------------

      USE damp1
      USE ddata
      USE iofile
      USE ndata
      USE tdata
      implicit double precision (a-h,o-z)
c
      logical mflg
      dimension du(*),dr(*),xml(*),xmu(*),xmd(*),xdl(*),xdu(*),xdd(*),
     +          ur(nneq,*),jd(*)

c.... at n
      if(ityp.eq.1) then ! compressed vectors neq
c....   modify for  mass term R=R-M*a^i
                  imas=2
        if (mflg) imas=1
        call promul(xml,xmu,xmd,ur(1,nrm),dr,jd,neq,imas,2)
c
c....   modify for  damping term R=R-D*v^i
        if(nc.gt.1) then
                     imas=2
          if (flgda) imas=1
          call promul(xdl,xdu,xdd,ur(1,nrc),dr,jd,neq,imas,2)
        end if

c.... at n+1/2
      else if(ityp.eq.2) then    ! compressed vectors neq
c....   v_n=ur(3), a_n=ur(4), du_n+1/2 has to be calculated before

c....   modify for  mass term R=R-M*a^i,   a at n+1/2 in ur(5)
        do n = 1,neq
          ur(n,5) = 4.d0*c1*du(n) - 2.d0*c1*dt*ur(n,3) + c3*ur(n,4)
        end do
                  imas=2
        if (mflg) imas=1
        call promul(xml,xmu,xmd,ur(1,5),dr,jd,neq,imas,2)
c
c....   modify for  damping term R=R-D*v^i,     v at n+1/2 in ur(5)
        if(nc.gt.1) then
          do n = 1,neq
            ur(n,5) = 2.d0*c2*du(n) + c4*ur(n,3) + 0.5d0*c5*ur(n,4)
          end do
                     imas=2
          if (flgda) imas=1
          call promul(xdl,xdu,xdd,ur(1,5),dr,jd,neq,imas,2)
        end if
      end if
      return
      end
c

      subroutine resmod(c,a,b,neq,iadd)
c-----------------------------------------------------------------------
c     Purpose: modification of right hand side for lumped matrix
c     c(i) = c(i) +- a(i,i)*b(i)
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension c(*),a(*),b(*)
      if(iadd.eq.1) then
        do n = 1,neq
          c(n) = c(n) + a(n) * b(n)
        end do
      else if(iadd.eq.2) then
        do n = 1,neq
          c(n) = c(n) - a(n) * b(n)
        end do
      end if
      return
      end
c
      subroutine du_at_n12(u,uq,ur,id,nneq,ityp)
c-----------------------------------------------------------------------
c
c     Purpose: calculate displacements for half step residual
c       ityp=1: calculate  uq(nneq) at n+1/2
c       ityp=2: calculate  du(neq)  at n+1/2
c
c      Inputs:
c         ur(*) - v_n(neq) in ur(3), a_n(neq) in ur(4)
c         ityp  - 1,2
c         id(*) - Equation numbers for each active dof
c         id    - Total number of characters in x-array
c         nneq  - length of vector uq, could be neq or nneq
c
c      Outputs:
c         uq(*) - u_n+1/2(ityp=1) or du_n+1/2(ityp=2)
c
c      Comment:   v_n and a_n always in ur3 and ur4!!
c
c-----------------------------------------------------------------------
      USE tdata
      implicit double precision (a-h,o-z)
      dimension u(nneq,*),uq(*),ur(nneq,*),id(*)
      c18=1.0d0/8.0d0
      c38=3.0d0*c18*dt
      c116=0.5d0*c18*dt*dt

      if(ityp.eq.1) then   ! full length uq(nneq)
        do n =1,nneq
          du    = u(n,2)
          uq(n) = u(n,1)-du
          j = id(n)
          if (j.gt.0) then
            uq(n) = uq(n) + c18*du + c38*ur(j,3) + c116*ur(j,4)
          else
cww??       uq(n) = uq(n) + 0.5d0*du
            uq(n) = uq(n) + c18*du
          end if
        end do

      else if(ityp.eq.2) then  ! compressed length uq(neq)
        call pzero(uq,nneq)
        nn=0
        do n =1,nneq
          j = id(n)
          if (j.gt.0) then
            nn    = nn+1
            du    = u(n,2)
            uq(nn) = c18*du + c38*ur(j,3) + c116*ur(j,4)
          end if
        end do
      end if

      return
      end
c

