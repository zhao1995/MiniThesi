
c----------------------------------------------------------------------+
c                                                                      |
c     Special Subroutines for Yield-Line Element  jw020605             |
c                                                                      |
c     Part 4: Auxiliary Algorithms                                     |
c                                                                      |
c----------------------------------------------------------------------+

c=======================================================================
c====  GA CODES                                                     ====
c=======================================================================

c=======================================================================

      subroutine gendim(vmut,nprp,nrec,npop,maxpop)

c----------------------------------------------------------------------+
c     Dimension Genetic Algorithm                                      |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------

      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

      integer nprp,nrec,npop
      integer vmut(nprp,nrec)

c---- nprp    : number of properties per individual
c---- nrec    : number of genes per individual (recombination)
c---- npop    : number of individuals (population)
c---- vmut    : allowed mutation steps
c               (0=no step,1=step allowed)

c----- Variables -------------------------------------------------------

      integer nmut,npar

c---- nmut    : number of mutation variables
c---- npar    : number of parents

c==== PROGRAM ==========================================================

c---- Number of Mutational DOFs ----------------------------------------
c     Every DOF Step per Gene Has to Be Represented

      nmut=1 ! first individual
      do i1=1,nprp
        do i2=1,nrec
          if (vmut(i1,i2).gt.0) nmut=nmut+2  ! 2 step directions
        enddo
      enddo

c---- Calculate Number of Parent Pair Constellations -------------------

c      npar=(fac(nmut)/(fac(npai)*fac(nmut-npai)))

c      print*, 'npar=',npar

      npar=0

c---- Population Number ------------------------------------------------

      npop=max(nmut,npar)
      npop=npop
c      print*, 'npop=', npop,nprp,nprp*npop

      if (npop.gt.maxpop) then
        print*, '  * Maximum Population reached!'
        npop=maxpop
      endif

c==== End of Subroutine ================================================

c  999 continue
      return
      end


c=======================================================================

      subroutine genalg(pop,vmut,step,nprp,nrec,npop,npai,nseq,
     &                  iswmut,iswrec,iswcmb,iswgen,isw)

c----------------------------------------------------------------------+
c     Genetic Mutation and Recombination                      jw050603 |
c                                                                      |
c     isw=1 input : pd, vmut, step, dmut, n***, iswmut                 |
c     isw=1 output: pd                                                 |
c     isw=2 input : pd, n***, iswrec, iswcmb, iswgen                   |
c     isw=2 output: pd                                                 |
c                                                                      |
c     subroutines:                                                     |
c       gencmb - generate combinations of parents and genes            |
c     function:                                                        |
c       fac    - calculate faculty                                     |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------

      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

c      common /comga/ lmut

      integer nprp,nrec,npop,npai,nseq
      integer iswmut,iswrec,isw
      integer vmut(nprp,nrec) !,dmut(nprp,nrec)!,vrng(2,nprp,nrec)
      real*8 pop(nprp,nrec,npop),step(nprp)

c---- nprp    : number of properties per individual
c---- nrec    : number of genes per individual (recombination)
c---- npop    : number of individuals (population)
c---- npai    : number of pair members (=2) -> number of children
cc---- npar    : number of recombinable parents (npar/npai=integer!)
c---- nseq    : number of sequence/gene crossovers (=2)(>=npai)
c---- iswmut  : switch for mutation variant (if step>0 only addition)
c               (1=increase/decrease single properties (default),
c                2=increase/decrease property group,
c                3=swap properties)                                     !only 1+2 implemented yet
c---- iswrec  : switch for gene recombination variant
c               (1=swap single genes,2=sequential gene swapping)
c---- iswcmb  : switch for parents recombination
c               (1=single appearance (default),
c                2=multiple appearance 1x per new individual,
c                3=arbitrary appearance)
c---- iswgen  : switch for gene recombination
c               (1=exchange fixed gene no for parent couple (default)
c                  makes sense for iswcmb=1,
c                2=exchange arbitrary gene no)
c---- isw     : evolution task
c               (1=predefine/mutate,2=recombine)

c---- vmut    : allowed mutation steps
c               (0=no step,1=step allowed)
c---- step    : step width per property (same for each gene)
c---- dmut    : dependence of dof steps in mutation of gene properties
c               (0=no dependence,else=scaling factor)
c---- pop     : population

c==== Variables ========================================================

      integer ipar,iseq
      integer vcmb(npai,npop),vgen(nseq-1,npop),vpop(npop)
      integer hv(max(npai,nseq))
      real*8 opop(nprp,nrec,npop)

c---- ipar    : number of parent entry
c---- iseq    : number of gene sequence entry
c---- vcmb    : vector of parent individuals to be combined
c---- vgen    : vector of genes to be swapped
c---- vpop    : vector for individuals' numbers
c---- hv      : helping vector
c---- opop    : old population

c      logical lmut

c---- lmut    : not first mutation/
c               optimal mutation of individual reached

c==== PROGRAM ==========================================================

c---- Transfer to Process ----------------------------------------------

c  100 continue
      goto (1000,2000), isw

c==== isw=1: Genesis or Mutation =======================================

 1000 continue

c---- Create New Generation --------------------------------------------
c     BAUSTELLE: dmut beachten!

      if (iswmut.eq.2) then
        i1=0
 1010   continue
        do i3=1,nprp
          if (step(i3).ge.0.d0.and.step(i3).le.0.d0) cycle
          i1=i1+1
          if (i1.gt.npop) return
          do i2=1,nrec
            if (vmut(i3,i2).eq.1.and.
     &         (step(i3).gt.0.d0.or.step(i3).lt.0.d0)) then
              pop(i3,i2,i1)=pop(i3,i2,i1)+step(i3)
            endif
          enddo
          i1=i1+1
          if (i1.gt.npop) return
          do i2=1,nrec
            if (vmut(i3,i2).eq.1.and.step(i3).lt.0.d0) then
              pop(i3,i2,i1)=pop(i3,i2,i1)-step(i3)
            endif
          enddo
        enddo
        goto 1010
      else
        i1=0
        do i2=1,nrec
          do i3=1,nprp
            if (vmut(i3,i2).eq.1.and.
     &         (step(i3).gt.0.d0.or.step(i3).lt.0.d0)) then
              i1=i1+1
              if (i1.gt.npop) return
              pop(i3,i2,i1)=pop(i3,i2,i1)+step(i3)
              i1=i1+1
              if (i1.gt.npop) return
              if (step(i3).lt.0.d0) then
                pop(i3,i2,i1)=pop(i3,i2,i1)-step(i3)
              endif
            endif
          enddo
        enddo
      endif
      return

c==== isw=2: Recombination =============================================

 2000 continue

c---- Define Combination Vector for Parents ----------------------------

      call pzero(vcmb,npai*npop)
      call pzero(vpop,npop)

      n1=int(npop/npai)
      n2=npai*int(npop/npai)

      if (iswcmb.eq.2.or.iswcmb.eq.3) then
        do i1=1,npop
          call gencmb(hv,npai,1,npop,iswcmb-1)
          do i2=1,npai
            vcmb(i2,i1)=hv(i2)
          enddo
        enddo

      else
        call gencmb(vpop,npop,1,npop,1) ! 1x each individual per population
        i4=1
        do i1=0,n1-1
          do i2=1,npai
            do i3=1,npai
              vcmb(i2,i1*npai+i3)=vpop(i4)
            enddo
            i4=i4+1
          enddo
        enddo
        if (n2.lt.npop) then
          do i1=n2+1,npop
            do i2=1,npai
              vcmb(i2,i1)=vpop(i4)
            enddo
            i4=i4+1
          enddo
        endif
      endif

c---- Rearrange vcmb ---------------------------------------------------
c     (for iswcmb=1)

      do i1=1,npop
        i3=i1
        do i2=1,npai
          i3=i3+1
 2010     continue         !while-loop!
          if (i3.gt.npai) then
            i3=i3-npai
            goto 2010
          endif
          hv(i2)=vcmb(i3,i1)
        enddo
        do i2=1,npai
          vcmb(i2,i1)=hv(i2)
        enddo
      enddo

c---- Define Combination Vector for Genes ------------------------------

      call pzero(vgen,(nseq-1)*npop)

      if (iswgen.eq.2) then
        do i1=1,npop
          call gencmb(hv,nseq,2,nrec,1)
          do i2=1,nseq
            vgen(i2,i1)=hv(i2)
          enddo
        enddo
      else
        do i1=0,n1-1
          call gencmb(hv,nseq,2,nrec,1)
          do i2=1,npai
            do i3=1,nseq-1
              vgen(i3,i1*npai+i2)=hv(i3)
            enddo
          enddo
        enddo
        if (n2.lt.npop) then
          do i1=n2+1,npop
            call gencmb(hv,nseq,2,nrec,1)
            do i2=1,nseq-1
              vgen(i2,i1)=hv(i2)
            enddo
          enddo
        endif
      endif

c---- Rearrange vgen ---------------------------------------------------

      do i1=1,npop
        do i2=1,nseq-1
          hv(i2)=vgen(i2,i1)
        enddo
        do i2=nseq-1,1,-1
          i0=0
          do i3=1,nseq-1
            if (hv(i3).gt.i0) then
              i0=hv(i3)
              i4=i3
            endif
          enddo
          vgen(i2,i1)=hv(i4)
          hv(i4)=0
        enddo
      enddo

c---- Recombination of Individuals' Genes ------------------------------

c      call pzero(opop,nprp*nrec,npop)
      opop=0

      opop=pop
      if (iswrec.eq.1) then  !single exchange
        do i1=1,npop
          ipar=1
          iseq=1
          do i2=1,nrec
            if ((i2.eq.vgen(iseq,i1)).and.(iseq.lt.nseq)) then
              ipar=ipar+1
              iseq=iseq+1
              if (ipar.gt.npai) ipar=2
              do i3=1,nprp
                pop(i3,i2,i1)=opop(i3,i2,vcmb(ipar,i1))
              enddo
            else
              do i3=1,nprp
                pop(i3,i2,i1)=opop(i3,i2,vcmb(1,i1))
              enddo
            endif
          enddo
        enddo
      else                   !sequential exchange
        do i1=1,npop
          ipar=1
          iseq=1
          do i2=1,nrec
            if ((i2.eq.vgen(iseq,i1)).and.(iseq.lt.nseq)) then
              ipar=ipar+1
              iseq=iseq+1
              if (ipar.gt.npai) ipar=1
            endif
            do i3=1,nprp
              pop(i3,i2,i1)=opop(i3,i2,vcmb(ipar,i1))
            enddo
          enddo
        enddo
      endif
c      return

c==== End of Subroutine ================================================

c9999  continue
      return
      end


c=======================================================================

      subroutine gencmb(v,nx,ns,ne,isw)

c----------------------------------------------------------------------+
c     Genetic Combination Vector                                       |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------

      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

      integer v(nx),ns,ne,isw

c---- v       : combination vector
c---- nx      : number of exchanges (combinations)
c---- ns      : start number
c---- ne      : end number
c---- isw     : switch for kind of recombination
c               (1=single appearance,2=multiple appearance)

c---- Variables --------------------------------------------------------

      integer n,ndim
      real*8 q

c---- n       : random number (integer)
c---- ndim    : number dimension
c---- q       : random number (real)

c==== PROGRAM ==========================================================

c---- Dimension Search Vector ------------------------------------------

      v=0
      n=ne
      ndim=1
  100 continue
      if (n.ge.10) then
        n=n/10
        ndim=ndim+1
        goto 100
      else
        goto 200
      endif

c---- Transfer to Process ----------------------------------------------

  200 continue
      goto (1000,2000), isw

c==== isw=1: Single Appearance of Entries ==============================

 1000 continue
      i1=1
 1100 continue
      call random_number(q)
      n=int(q*(10**ndim))
      if ((n.eq.0).or.(n.lt.ns).or.(n.gt.ne)) goto 1100
      do i2=1,nx
        if (v(i2).eq.n) goto 1100
        if (v(i2).eq.0) then
          v(i2)=n
          i1=i1+1
          exit
        endif
      enddo
      if (i1.le.nx) goto 1100
      return

c==== isw=2: Multiple Appearance of Entries ============================

 2000 continue
      i1=1
 2100 continue
      call random_number(q)
      n=int(q*(10**ndim))
      if ((n.eq.0).or.(n.lt.ns).or.(n.gt.ne)) goto 2100
      v(i1)=n
      i1=i1+1
      if (i1.le.nx) goto 2100
c      return

c==== End of Subroutine ================================================

c9999  continue
      return
      end


c=======================================================================

      integer function fac(n)

c----------------------------------------------------------------------+
c     Calculate Faculty                                                |
c----------------------------------------------------------------------+

      integer i,i1

      i=1
      do i1=1,n
        i=i*i1
      enddo

      n=i
      fac=n

c==== End of Function ==================================================

c  99  continue
      return
      end


cc=======================================================================
c
c      subroutine randomi(a,e,i)
c
cc----------------------------------------------------------------------+
cc     Create Random Number Integer                                     |
cc----------------------------------------------------------------------+
c
c      real*8 q
c      integer a,e,i,ndim
c
cc---- Dimension Search Area --------------------------------------------
c
c      n=int(e)
c      ndim=1
c   10 continue
c      if (n.ge.10) then
c        n=n/10
c        ndim=ndim+1
c        goto 10
c      endif
c
c   20 continue
cc      call random_seed
c      call random_number(q)
c      i=(q*(10**ndim))
c      if ((i.lt.a).or.(i.gt.e)) goto 20
c
cc==== End of Subroutine ================================================
c
c  99  continue
c      return
c      end
c
c
cc=======================================================================
c
c      subroutine randomr(a,e,r)
c
cc----------------------------------------------------------------------+
cc     Create Random Number Real                                        |
cc----------------------------------------------------------------------+
c
c      real*8 a,e,q,r
c      integer ndim
c
cc---- Dimension Search Area --------------------------------------------
c
c      n=int(e)
c      ndim=1
c   10 continue
c      if (n.ge.10) then
c        n=n/10
c        ndim=ndim+1
c        goto 10
c      endif
c
c   20 continue
cc      call random_seed
c      call random_number(q)
c      r=(q*(10**ndim))
c      if ((r.lt.a).or.(r.gt.e)) goto 20
c
cc==== End of Subroutine ================================================
c
c  99  continue
c      return
c      end



c=======================================================================

      subroutine randomi(a,e,i)

c----------------------------------------------------------------------+
c     Create Random Number Integer                                     |
c----------------------------------------------------------------------+

      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

      common /randomb/ rran(100000),iran,lran

c---- rran     : array of random numbers
c---- iran     : pointer
c---- lran     : new creation of rran

      logical lran
c      real*8 q
      integer a,e,i !,ndim

c==== PROGRAM ==========================================================

c---- Initialize Array -------------------------------------------------

  100 continue
      if (.not.lran) then
        print*, 'i: Have to initialize'
        iran=1
c        call random_seed
        do i1=1,100000
          call random_number(rran(i1))
        enddo
        lran=.true.
        print*, 'i: Initialized'
c      else
c        print*, '----geht so!-'
      endif

c---- Give Random Number Out Of rran -----------------------------------

  200 continue
      if (iran.eq.100000) then
        lran=.false.
        print*, 'i: New Initialization'
        goto 100
      endif

      if (int(10*rran(iran)).gt.a.and.int(10*rran(iran)).lt.e) then
        i=int(10*rran(iran))
        iran=iran+1
      else
        iran=iran+1
        goto 200
      endif

c==== End of Subroutine ================================================

c  999 continue
      return
      end


c=======================================================================

      subroutine randomr(a,e,r)

c----------------------------------------------------------------------+
c     Create Random Number Real                                        |
c----------------------------------------------------------------------+

      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

      common /randomb/ rran(100000),iran,lran

c---- rran     : array of random numbers
c---- iran     : pointer
c---- lran     : new creation of rran

      logical lran
      real*8 a,e,r
c      integer ndim


c==== PROGRAM ==========================================================

c      ifac=3

c---- Initialize Array -------------------------------------------------

  100 continue
      if (.not.lran) then
c        print*, 'r: Have to initialize'
        iran=1
c        call random_seed
        do i1=1,100000
          call random_number(rran(i1))
c          rran(i1)=int(rran(i1)*10**ifac)/(10**ifac)
c          print*, 'hab',rran(i1)
        enddo
        lran=.true.
c        print*, 'r: Initialized'
      endif

c---- Give Random Number Out Of rran -----------------------------------

  200 continue
      if (iran.eq.100000) then
        lran=.false.
c        print*, 'r: New Initialization'
        goto 100
      endif

      if (rran(iran).gt.a.and.rran(iran).lt.e) then
        r=rran(iran)
        iran=iran+1
      else
        iran=iran+1
        goto 200
      endif

c==== End of Subroutine ================================================

c  999 continue
      return
      end

c=======================================================================
c=======================================================================











c=======================================================================
c=======================================================================
c=======================================================================

      subroutine ygauss(a, nr, nc, isort, irg)

      integer  nr, nc, i, ii, j, jj
      integer zeilo, zeilu
      integer irg, isortalt(nr), isort(nr)
      real*8 a(nr, nc), anew(nr, nc), asort(nr,nc), lambda, fac, eps
c      real*8 b(nr), x(nr)

c---- a(nr,nc):      Ausgangs_Matrix
c---- nr, nc  :      Anzahl: Zeilen, Spalten
c---- b, x    :      jw
c---- isort(alt):    i-te Zeile einsortiert in neue Zeile isort(i)
c---- irg     :      Rang der sortierten Matrix

c---- i, j :      "Neben"-Laufvariablen (Zeile, Spalte)
c---- ii, jj  :  "Haupt"-Laufvariablen (Zeile, Spalte)
c---- zeilo   :      Zeilenzeiger von oben, steigt
c---- zeilu   :      Zeilenzeiger von unten, sinkt
c---- anew(nr,nc):   Arbeits_Matrix
c---- asort(nr,nc):	Arbeits_Matrix (bzgl. Sortieren)
c---- lambda  :      Multiplikationsfaktor der Zeilen untereinander
c---- fac     :      Multiplikationsfaktor fuer Normierung
c---- eps     :  Fehler des Nullwerts

c-----------------------------------------------------------------------
c     Fehler-Null-Wert festsetzen
      eps=1.d-12

c     Variablenwerte vorgeben
      irg=-1
      do i=1, nr
	 isortalt(i)=i
      enddo
      do i=1, nr
	 do j=1, nc
	   anew(i,j)=a(i,j)
	 enddo
      enddo

c-----------------------------------------------------------------------
c     GAUSS:   Sortieren, Abbruchkriterium, Eliminieren, ...
      jj=0
      do ii=1,nr

c>>>  Sortieren
 50      continue
      	  zeilu=nr
	  zeilo=ii
 	  jj=jj+1
	  do i=ii, nr
	     if (abs(anew(i,jj)).le.eps) then
	        anew(i,jj)=0.d0                      !|x|<eps then 0.d0
	        isort(zeilu)=isortalt(i)
	        do j=1, nc
		    asort(zeilu,j)=anew(i,j)
	        enddo
	        zeilu=zeilu-1                        !von unten auffüllen
	     else
	        isort(zeilo)=isortalt(i)
	        do j=1, nc
		    asort(zeilo,j)=anew(i,j)
	        enddo
	        zeilo=zeilo+1                        !von oben auffüllen
	     endif
         enddo
	  do i=ii, nr
	     do j=1, nc
	        anew(i,j)=asort(i,j)
	     enddo
	  enddo
cjw	  print *, 'sortiert'
	  do i=1, nr
	     isortalt(i)=isort(i)
cjw	     print *, isort(i)
cjw	     print *, (anew(i,j), j=1, nc)
	  enddo
c<<<  Sortieren

c>>>  Abbruchkriterium
         if (abs(anew(ii,jj)).le.eps) then
            anew(ii,jj)=0.d0
            if (jj.ge.nc) then
               goto 100				!Elimination beendet
            else
               goto 50				!Elimination fortfuehren
            endif
         else
c<<<  Abbruchkriterium

c>>>  Eliminieren, vorher Masterzeile normieren
c     Normieren
            if (anew(ii,jj).gt.1.d0.or.anew(ii,jj).lt.1.d0) then
               fac=(1.d0)/(anew(ii,jj))
               do j=1, nc
                  anew(ii,j)=anew(ii,j)*fac
               enddo
            endif
c     Eliminieren
            do i=ii+1, nr
               if (anew(i,jj).gt.0.d0.or.anew(i,jj).lt.0.d0) then
	        lambda=anew(ii,jj)/anew(i,jj)
cjw	        print *, 'lambda', i, ':', lambda
	        do j=1, nc
		    anew(i,j)=anew(ii,j)-lambda*anew(i,j)
		    if (abs(anew(i,j)).le.eps) then	! abs()<eps=0.d0
		       anew(i,j)=0.d0
		    endif
	        enddo
	        endif
	     enddo
         endif
      enddo
c<<<  Eliminieren

c     Elimintation beendet
 100  continue

c-----------------------------------------------------------------------
c     Rang bestimmen
      i=1
      j=1
 201  continue					!Pfad "<>Null" verfolgen
      if (abs(anew(i,j)).le.eps) then
         j=j+1
         if (j.gt.nc) then
            goto 202  					!Abbruch, Spalte
         endif
         goto 201
      else
         i=i+1
         if (i.gt.nr) then
            goto 202  					!Abbruch, Zeile
         endif
         goto 201
      endif
 202  continue						!Auswertung Pfad
      if (i.ge.nr) then
         irg=nr
      else
         if (j.ge.nc) then
            irg=i-1
         else
            print *, 'Rang-error'
         endif
      endif

cjw>>
c---- Kontrolle auf leere Spalten --------------------------------------

c      call plotmatrix(anew,nr,nc,'ygauss')


cjw<<
c-----------------------------------------------------------------------
c     Uebergabe von anew()->a()
      do i=1, nr
	  do j=1, nc
	     a(i,j)=anew(i,j)
	  enddo
      enddo

      end
