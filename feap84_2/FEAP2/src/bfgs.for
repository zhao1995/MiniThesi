      subroutine iterat(stiff,jp,u,rsd,oldrsd,d,t,
     1      accrcy,v,w,prt,f0,f,id,nbfgs,stol,nneq)
c----------------------------------------------------------------------
c      Purpose: BFGS algorithm  for solution to equation

c          R(u) = F - P(u) = 0
c          Variable number of vectors
c          Max. no. of vectors = nbfgs ( < = 15 ) (input data)
c         BFGS algorithm (Bathe's paper, 1980)
c          Version 8/11/86 with variable number of vectors
c          Max. no. of vectors = nbfgs ( < = 15 ) (input data)
c
c     CALL sequence
c      pmacr3.f -> bfgs.f (iterat) -> formfe.f (isw = 6)
c      -> pform.f -> elmlib.f ->  elmtxx.f

c     Description
c      routine will return solution and residual to equation:
c          K(u)*du - b = 0
c      Routine evaluates RHS of equation through routine operat.
c
c      Inputs:
c       stiff
c       jp(*)       - Pointer array for row/columns of tangent
c       u(*)        - Nodal solution values
c       rsd        - Pointer for current residual
c       oldrsd(*)   - Old residual
c       accrcy
c       prt         - Print if true
c       f0
c       f
c       id(*)       - Equation numbers for each dof
c       stol        - Solution tolerance for line search

c      Outputs:
c       d(*)        - Incremental solution vector for step
c       v(*)        - BFGS vectors
c       w(*)        - BFGS vectors
c       nbfgs       - Number BFGS steps

c      Scratch:
c       t(*)        - Temporary storage vectors
c----------------------------------------------------------------------
      USE cdata
      USE iofile
      USE rdata
      implicit double precision(a-h,o-z)
c
      logical accrcy,prt
      common /iupdt/ iform
c.... u(1) ^ u(nneq) = total displ.
c.... u(1+nneq) ^ u(2*nneq) = total incr. displ.
c.... u(1+2*nneq) ^ u(3*nneq) = current incr. displ.
      dimension stiff(*),jp(*),u(*),rsd(*),oldrsd(*),d(*),
     *    v(*),w(*),f0(*),f(*),id(*),t(*)
c.... iform : no. of "formfe" routine executed by bfgs (by "operat")
      iform=0
c.... itmax : max. no. of iteration allowed
c.... s     : initial step size (=1)
c.... G (s) : scalar function for search direction
      etol = tol
      itmax = nbfgs
      nupd  = 0
      s     = 1.0d0
      g0    = 0.0d0
      g     = 0.0d0
      nneq2 = nneq + nneq
c.... initialization
c.... rsd : residual R(u+s*d)
      call pzero(oldrsd,nneq)
      call pzero(rsd,nneq)
      call pzero(d,nneq)
c     compute residual-(0) for du = 0
      call operat(rsd,d,f0,f,id,u,t,nneq)
      iform = iform + 1
c.... loop for momentum balance iteration
      do 100 i=1,itmax
      if(.not.prt) go to 30
      write(iow,2000)  i
      if(ior.lt.0) write(*,2000)  i
      do 20 l=1,neq
20    write(iow,2001) l,rsd(l),oldrsd(l),d(l)
      if(ior.lt.0) write(*,2001) l,rsd(l),oldrsd(l),d(l)
30    continue
c.... compute the search direction (d) by factorized form (v,w)
c     d : not including step-size
      call dfind(stiff,jp,d,rsd,oldrsd,nupd,g0,g,s,
     *       neq,v,w,nbfgs)
c.... do line search if necessary
      s=1.0d0
c     compute new residual (rsd)
      call operat(rsd,d,f0,f,id,u,t,nneq)
      iform = iform + 1
      g0=ddot(neq,d,1,oldrsd,1)
      g =ddot(neq,d,1,rsd,1)
c     chek : tolerance condition for line search (step-size)
      chek=stol*dabs(g0)
      if(dabs(g).gt.chek) call serch2(g0,g,f0,f,id,rsd,u,d,
     *  stol,t,neq,nneq,s)
c.... NOTE: d back from serch2 = current incremental u-(i) (including s)
c.... UPDATE u (disp., incremental disp., & current incr.)
c.... refer to update.f
      do 40 n = 1,nneq
      j = id(n)
      if(j.gt.0) then
c.... for the active degrees-of-freedom compute values from solution
c.... where 'd(j)' is the increment of 'u' for active dof 'j'.
c.... for the fixed degrees-of-freedom: do nothing
        u(n)      = u(n) + d(j)
        u(n+nneq) = u(n+nneq) + d(j)
        u(n+nneq2) = d(j)
      end if
40    continue
c.... Convergence criterion: see Bathe [1980], p.72,73
c....  (1) displ. criterion: not reliable, not used
c....  (2) force residual criterion : print out for reference
c....  (3) energy criterion : used as the only criterion
c.... compute norms
c     norm of delta = incr. u vector
cww      dnorm=sqrt(dot(d,d,neq))
c     norm of residual vector
      rnorm=sqrt(ddot(neq,rsd,1,rsd,1))
c.... check for accuracy (energy criterion): eq. (44) p.73, Bathe
      if (i.eq.1) then
        oldeng = dabs(g0*s)
        energy = oldeng
      else
        energy = dabs(g0*s)
      end if
      accrcy = energy.le.etol*oldeng.or.energy.le.etol
      if(accrcy) go to 900
      write(iow,2004) rnorm,energy
      if(ior.lt.0) write(*,2004) rnorm,energy
100   continue
900   write(iow,2002) iform,rnorm,energy
      if(ior.lt.0) write(*,2002) iform,rnorm,energy
      i = min(i,itmax)
      write(iow,2003) i
      if(ior.lt.0) write(*,2003) i
      return
c.... formats
2000  format(4x,'start of iteration no. ',i5/
     1       4x,'d.o.f.',2x,'rsd. vector',5x,'oldrsd',
     2       10x,'search vector')
2001  format(2x,i5,3d16.5)
2002  format(2x,'**macro instruction form executed by bfgs',i3,' times',
     1      /6x,'   rnorm = ',g14.5/
     2       6x,'  energy = ',g14.5)
2003  format(2x,i3,' BFGS iterations required to converge')
2004  format(6x,'   rnorm = ',g14.5/
     2       6x,'  energy = ',g14.5)
      end
c
      subroutine serch2(g0,g,f0,f,id,rsd,u,d,stol,t,nqe,nneq,step)
c----------------------------------------------------------------------
c      Purpose: makes a line search in direction d and returns the
c               steplength step
c
c      Inputs:

c      Outputs:
c----------------------------------------------------------------------
      USE cdata
      USE iofile
      implicit double precision(a-h,o-z)
      dimension f0(*),f(*),id(*),rsd(*),u(*),d(*),t(*)
      linmax=10
      sb=0.0d0
      sa=1.0d0
c.... find bracket on zero
      if(g*g0.gt.0.0d0) then
         call drawmess(' no line search - ga,gb both positive',1,0)
         return
      end if
      j = 0
      gb=g0
      ga=g
10    j = j + 1
      step = sa - ga*(sa-sb)/(ga-gb)
      g = gamma2(f0,f,id,u,rsd,d,t,step,nqe,nneq)
      gb = 0.5*gb
      if (g*ga.lt.0.0d0) then
         sb = sa
         gb = ga
      end if
      sa = step
      ga = g
      write(iow,6666) j,step,g
      if(ior.lt.0) write(*,6666) j,step,g
      if (j.ge.linmax) go to 20
      if(dabs(g).gt.stol*dabs(g0)) go to 10
      if(dabs(sb-sa).gt.stol*0.5*(sa+sb)) go to 10
20    do 30 j = 1,nqe
         d(j) = step*d(j)
30    continue
      return
6666  format('j,step,g',i5,3e15.5)
      end

c
      double precision function gamma2(f0,f,id,u,dr,du,t,s,nqe,nneq)
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      USE cdata
      USE fdata
      USE mdata
      USE ndata
      USE prlod
      USE tdata
      implicit double precision (a-h,o-z)
      logical fa,tr
      common /iupdt/ iform
      dimension f0(*),f(*),id(*),u(*),dr(*),du(*),t(*)
      data fa,tr/.false.,.true./
c.... get a search displacement
      nneq2 = nneq + nneq
      do 100 n = 1,nneq
      j = id(n)
      if(j.gt.0) then
        t(n)      = u(n) + s*du(j)
        t(n+nneq) = u(n+nneq) + s*du(j)
        t(n+nneq2)= s*du(j)
      else
        db        = f0(n) + f(n)*prop
        t(n+nneq2)= db - u(n)
        t(n+nneq) = u(n+nneq) -u(n) + db
        t(n)      = db
      end if
100   continue
c.... compute a residual
      call ploads(u,dr,prop,.false.,.false.,.false.)
      call pload(psid,gloa,glo0,dr,nneq,prop)
c.... update the residual for lumped mass inertial effects
cww   hflgu/h3flgu fuer history Werte setzen?
      call formfe(t,dr,dr,dr,fa,tr,fa,fa,6,1,numel,1)
      iform = iform + 1
      if(fl(9)) then
        call pmove (trans(nw),t,nneq)
        call pmove (massm,t(nneq+1),nqe)
        ctem = s*c2
        do 110 n = 1,nneq
          j = id(n)
          if(j.gt.0) dr(j) = dr(j) - t(nneq+j)*(t(n)+ctem*du(j))
110     continue
      end if
c.... compute the value of gamma
      gamma2 = ddot(nqe,du,1,dr,1)
      return
      end

c
      subroutine dfind(stiff,jp,d,rsd,oldrsd,nupd,g0,g,
     1     s,neq,v,w,nbfgs)
c----------------------------------------------------------------------
c      Purpose: find a new search direction using BFGS updating
c               method in factored form
c
c      Inputs:

c      Outputs:
c----------------------------------------------------------------------
      USE iofile
      implicit double precision(a-h,o-z)
      logical up
      dimension stiff(*),jp(*),d(*),rsd(*),oldrsd(*),v(*),w(*)
      maxup  = nbfgs
      condmx = 1.0d5
      np1    = neq+1
c.... delgam = delta-(i) : gamma-(i)
c     dlkdl  = delta-(i) : K-(i-1) : delta-(i)
      delgam=s*(g0-g)
      dlkdl=s*s*g0
      up=delgam.gt.0.0d0.and.dlkdl.gt.0.0d0
      if(.not.up) go to 200
c     if G(0) > G(s) & G(0) > 0
c     not for i=1 from iterat
      stcond=sqrt(delgam/dlkdl)
c....test
      if(ior.lt.0) write(*,2000) stcond
2000  format(' --->stcond ',e12.5)
      fact2=s/delgam
c.... compute updating vectors v, w and put residual into d
      do 100 i=1,neq
      aux=oldrsd(i)-rsd(i)
      v(i)=-stcond*s*oldrsd(i)-aux
      w(i)=d(i)*fact2
      d(i)=rsd(i)
100   oldrsd(i)=rsd(i)
c.... check estimate on increase of condition number
      up=up.and.stcond.lt.condmx
      if(.not.up) go to 140
c.... save updating factors, with fact2 to be included later in w
c     nupd : no. of effective update loops
      call store(v,w,nupd+1,1,nbfgs)
c.... compute search direction: d
c.... the rightmost factor : w-(i) : rsd-(i)
      coef=fact2*g
c     b-(i): p. 1616 Matthies & Strang [1979]
      do 130 i=1,neq
130   d(i)=d(i)+coef*v(i)
140   continue
      go to 300
200   do 210 i=1,neq
      d(i)=rsd(i)
210   oldrsd(i)=rsd(i)
300   continue
      if(nupd.eq.0) go to 325
c.... right half of updating p. 1616 (a)
      do 320 i=1,nupd
      ii=nupd-i+1
      call store(v,w,ii,2,nbfgs)
      coef=ddot(neq,w,1,d,1)
      do 310 j=1,neq
310   d(j)=d(j)+coef*v(j)
320   continue
c.... p. 1616 (b)
c.... forward + backward substitution : c-(i) = K-(i-1)-inv : b-(i)
325   call dasol(stiff(np1),stiff(np1),stiff,d,jp,neq,engy)
      if(up) nupd=nupd+1
      if(nupd.eq.0) go to 350
c.... left half of updating : p. 1616 (c)
      do 340 i=1,nupd
      if(i.gt.1) call store(v,w,i,2,nbfgs)
      coef=ddot(neq,v,1,d,1)
      do 330 j=1,neq
330   d(j)=d(j)+coef*w(j)
340   continue
350   nupd=mod(nupd,maxup)
      return
      end
c
c
      subroutine operat(dr,du,f0,f,id,u,t,nneq)
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      USE cdata
      USE fdata
      USE mdata
      USE ndata
      USE prlod
      USE tdata
      implicit double precision(a-h,o-z)
      logical fa,tr
      dimension dr(*),du(*),f0(*),f(*),id(*),u(*),t(*)
      common /iupdt/ iform
      data fa,tr/.false.,.true./
c.... get a search displacement
c.... similar to gamma.f
      nneq2 = nneq + nneq
      do 100 n = 1,nneq
      j = id(n)
c.... t(1) ^ t(nneq) = total displ.
c.... t(1+nneq) ^ t(2*nneq) = total incr. displ.
c.... t(1+2*nneq) ^ t(3*nneq) = current incr. displ. only
      if(j.gt.0) then
        t(n)      = u(n) + du(j)
        t(n+nneq) = u(n+nneq) + du(j)
        t(n+nneq2)= du(j)
      else
        db        = f0(n) + f(n)*prop
        t(n+nneq2)= db - u(n)
        t(n+nneq) = u(n+nneq) -u(n) + db
        t(n)      = db
      end if
100   continue
c.... compute a residual
c     external loads
      call ploads(u,dr,prop,.false.,.false.,.false.)
      call pload(psid,gloa,glo0,dr,nneq,prop)
c.... update the residual for lumped mass inertial effects
c     substracting  the elmt stress
cww   hflguh/3flgu fuer history Werte setzen?
      call formfe(t,dr,dr,dr,fa,tr,fa,fa,6,1,numel,1)
      iform = iform + 1
      if(fl(9)) then
        call pmove (trans(nw),t,nneq)
        call pmove (massm,t(nneq+1),neq)
        do 110 n = 1,nneq
          j = id(n)
          if(j.gt.0) dr(j) = dr(j) - t(nneq+j)*(t(n)+c2*du(j))
110     continue
      end if
      return
      end
c
c
      subroutine store(v,w,np,itrn,nbfgs)
c----------------------------------------------------------------------
c     purpose: store vectors v,w in bfgsbs and viceversa
c----------------------------------------------------------------------
      USE cdata
      USE iofile
      USE isbfgs1
      USE yydata 
      implicit double precision(a-h,o-z)
      dimension v(*),w(*)
c
c.... check if np > nbfgs
      if (np.gt.nbfgs) then
        write(yyy,2000) nbfgs
        call drawmess(yyy,1,0)
2000    format(2x,'**WARNING** more than',i3,'iterations in BFGS algo')
cww       stop
        return
      end if
      ns=neq+neq
      go to (1,2),itrn
1     mtp = 1 + (np - 1)*ns    ! Überschreitung bfgsbs=15*neq npmax=15 ns=2*neq ww??
      call pmove(v,bfgsbs(mtp),neq)
      mtp = mtp + neq
      call pmove(w,bfgsbs(mtp),neq)
      return
2     mtp = 1 + (np - 1)*ns
      call pmove(bfgsbs(mtp),v,neq)
      mtp = mtp + neq
      call pmove(bfgsbs(mtp),w,neq)
      return
      end
c
