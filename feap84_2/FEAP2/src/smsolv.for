      subroutine mkptr2(is_csr)
c-----------------------------------------------------------------------
c
c     Purpose: Set basic values for symmetric sparse matrix solver
c
c     Inputs:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c
c-----------------------------------------------------------------------
      USE cdata
      USE iofile
      USE smpak
      USE doalloc
      implicit double precision (a-h,o-z)
      logical ldummy
c
      ldummy  = .true.
      is_csr = 1

c.... Speicher fuer SMPAK
      call ialloc(smsperm, neq,'SMSOLV-perm ',ldummy)
      call ialloc(smsperm1,neq,'SMSOLV-perm1',ldummy)
      call ralloc(smsidr1, neq,'SMSOLV-idr1 ',ldummy)

      return
      end
c
      subroutine dasol2(dr,jdiag,neq,aengy)
c----------------------------------------------------------------------
c
c.... Purpose: solve the system using SM-solver
c      input: neq           - number of equations
c      input: smsperm (neq) - permutation of rows and columns  
c      input: smsperm1(neq) - smperm^-1
c      input: jdiag         - ia
c      input: csrja         - ja
c      input: smtemp        - A
c      input: dr            - B 
c     output: dr            - x   Ax=B x=B possible  
c      input: nsp           - dim(rsp)
c      input: smsisp        - working
c      input: smirsp(nsp)   - for symbolic factorization
c     output: esp           - storage flag 
c      input: 3             - 'solv'
c     output: iflag         - info on 'solv'
c
c     scratch: smsidr1(neq)
c
c     W. Wagner BS UKA 04/09
c
c----------------------------------------------------------------------
      USE iofile
      USE iscsr
      USE smpak
      implicit double precision (a-h,o-z)
      integer esp
      dimension dr(*), jdiag(*)
      real*8,  allocatable, dimension(:) :: smstemp

      iflag = 0
      allocate(smstemp(1))

      call matcop(dr,neq,1,smsidr1)

      call sdrv(neq,smsperm,smsperm1,jdiag,csrja,smstemp,dr,dr,
     +      nsp,smsisp,smsirsp,esp,3,iflag)

      deallocate(smstemp)

      if (iflag.ne.0) then
        write(iow,*) '*** FEHLER IN SDRV (2) *** iflag=',iflag
        write(*  ,*) '*** FEHLER IN SDRV (2) *** iflag=',iflag
        call sderr(iflag,neq)
      end if
      aengy = ddot(neq,smsidr1,1,dr,1)

      return
      end
c

      subroutine dasbl2 (s,p,ld,jd,ns,aufl,bfl,b,ad,ja,iperm,lodrv)
c-----------------------------------------------------------------------
c
c      Purpose: Assemble full arrays for SM-solver
c
c      Inputs:
c         s(ns,ns) - Element array to assemble
c         p(ns)    - Element vector to assemble
c         ld(ns)   - Local to Global equation numbers for assemble
c         jd(*)    - Pointer array for upper/lower parts of A array.
c         ns       - Size of element arrays
c         aufl     - If true, assemble A array
c         bfl      - If true, assemble B vector
c         ja       - Pointer array CSR
c
c      Outputs:
c         b(*)     - Assembled right hand side B vector
c         ad(*)    - Assembled full part of A array incl. diag
c
c     W. Wagner BS UKA 04/09
c
c-----------------------------------------------------------------------
      USE cdata
      USE soltyp
      implicit double precision (a-h,o-z)
      logical aufl,bfl,col,lodrv
      dimension ad(*),b(*),s(ns,ns),p(*),ld(*),jd(*)
      dimension ja(*),iperm(*)

      if (lodrv) then
c....   *** nach erfolgter Sortierung... ***
        do i=1,ns
          ii=ld(i)
          if (ii.gt.0) then
            if (aufl) then
              do j=1,ns
                jc=ld(j)
                if (jc.gt.0) then
c...              *** wenn in umsortierter Matrix in UT ***
                  if (iperm(jc).ge.iperm(ii)) then
                    do nr=jd(ii),jd(ii+1)-1
                      if (ja(nr).eq.jc) then
                        !$OMP ATOMIC
                        ad(nr)=ad(nr)+s(j,i) ! critical OMP
                        exit
                      end if
                    end do
                  end if
                end if
              end do !(j)
            end if ! (aufl)
            if (bfl) then
              !$OMP ATOMIC
              b(ii) = b(ii) + p(i) ! critical OMP
            end if ! (bfl)
          end if ! (ii)
        end do !(i)
c
      else
c...    first time ...
        do i=1,ns
          ii=ld(i)
          if (ii.gt.0) then
            if (aufl) then
              do j=1,ns
                jc=ld(j)
                if (jc.ge.ii) then
                  do nr=jd(ii),jd(ii+1)-1
                    if (ja(nr).eq.jc) then
                      !$OMP ATOMIC
                      ad(nr)=ad(nr)+s(j,i) ! critical OMP
                      exit
                    end if
                  end do
                end if
              end do !(j)
            end if
            if (bfl) then
              !$OMP ATOMIC
              b(ii)=b(ii)+p(i) ! critical OMP
            end if
          end if
        end do !(i)
      end if
      return
      end
c


      subroutine sderr(iflag,neq)
c----------------------------------------------------------------------
c
c...  Purpose: writes error messages of SM-solver
c
c----------------------------------------------------------------------
        implicit none
        integer iflag,neq,fehler,zeile
        fehler=iflag/neq
        zeile=mod(iflag,neq)
        if (fehler.eq.2) then
          write(*,*) 'DUPLICATE ENTRY IN A'
        elseif (fehler.eq.6) then
          write(*,*) 'INSUFFICIENT STORAGE IN SSF'
        elseif  (fehler.eq.7) then
          write(*,*) 'INSUFFICIENT STORAGE IN SNF'
        elseif  (fehler.eq.8) then
          write(*,*) 'ZERO PIVOT'
        elseif  (fehler.eq.10) then
          write(*,*) 'INSUFFICIENT STORAGE IN SDRV'
        elseif  (fehler.eq.11) then
          write(*,*) 'ILLEGAL PATH SPECIFICATION'
        endif
        write(*,*) 'IN LINE',zeile
        return
        end
c

        subroutine oderr(iflag,neq)
c----------------------------------------------------------------------
c
c...  Purpose: writes error messages of SM-solver
c
c----------------------------------------------------------------------
        implicit none
        integer iflag,neq,fehler,zeile
        fehler=iflag/neq
        zeile=mod(iflag,neq)
        if  (fehler.eq.9) then
          write(*,*) 'INSUFFICIENT STORAGE IN MD'
        elseif  (fehler.eq.10) then
          write(*,*) 'INSUFFICIENT STORAGE IN ODRV'
        elseif  (fehler.eq.11) then
          write(*,*) 'ILLEGAL PATH SPECIFICATION'
        endif
        write(*,*) 'IN LINE',zeile
        return
        end
c

c----------------------------------------------------------------------
      subroutine datri2(ad,jd)
c----------------------------------------------------------------------
c
c.... Purpose: triangular decomposition using Sparse matrix solver
c
c     Comments:
c             - nneg = no. of neg. diags is calculated in solver
c             - only first call: use ODRV
c               (find a minimum degree ordering rows and columns of A)
c
c     W. Wagner BS UKA 04/09
c     update USE WW 11/14
c----------------------------------------------------------------------
      USE cdata
      USE fdata
      USE iofile
      USE iscsr
      USE ndata
      USE plong
      USE smpak
      USE soltyp
      USE doalloc
      implicit double precision (a-h,o-z)
      integer, allocatable, dimension(:) :: tempisp
      real*4  tary,tary1

      logical ldummy

c-----only once at first call ------------------------------------------
      if (.not.lodrv) then

        call etimef(tary)
        if(pfr)                write(iow,2014) tary
        if(pfr .and. ior.lt.0) write(*  ,2014) tary

        nsp=(3*neq+4*ljacsr)  ! see definition in ODRV
        call ialloc(tempisp,nsp,'SMSOLV-isp-temp',ldummy)! temporary field to determine sizes
        call ralloc(smsirsp,1,'SMSOLV-irsptemp',ldummy)

c....   find a minimum degree ordering rows and columns of A
        if(istyp.eq.1) then
c....     no permutation
          do i=1,neq
            smsperm(i) = i
            smsperm1(i)= i
          end do
          call odrv(neq,jd,csrja,ad,smsperm,smsperm1,nsp,tempisp,
     1              5,iflag) ! diagonal = 1.entry in row
        else
          call odrv(neq,jd,csrja,ad,smsperm,smsperm1,nsp,tempisp,
     1              4,iflag) ! diagonal = 1.entry in row
        end if
c....   errors
        if (iflag.ne.0) then
          write(iow,*) '*** FEHLER IN ODRV ***',iflag
          write(*  ,*) '*** FEHLER IN ODRV ***',iflag
          call oderr (iflag,neq)
        end if
c...    symbolic factorisation
        call sdrv(neq,smsperm,smsperm1,jd,csrja,ad,trans(nw),trans(nw),
     +            nsp,tempisp,smirsp,iesp,4,iflag)

c...    define new arrays smsisp smsrsp

cww     nsp  =     m(isp+neq)+4*neq+2   ! original code
cww     nsp  = tempisp(neq+1)+4*neq+2   ! accesses wrong position -> wrong size
        nsp  = tempisp(neq)  +4*neq+2   ! seems to be correct, now setting the real field to zero works
        
        call ialloc(smsisp,nsp,'SMSOLV-isp',ldummy)
        lisp = min(nsp,size(tempisp))  ! to be save, necessary?
        smsisp(1:lisp) = tempisp(1:lisp)  

cww     nsp =   (m(isp+neq-1)+3*neq+2)/2+1+2*neq+  m(isp+2*neq) ! original code
cww     nsp = (tempisp(neq-1)+3*neq+2)/2+1+2*neq+tempisp(2*neq) ! as consequence from above, but neq=1 leads to error
        nsp = (tempisp(neq  )+3*neq+2)/2+1+2*neq+tempisp(2*neq) ! original, ok for first example  
cww     nsp3= tempisp(2*neq) ! could be negative in case of  large systems!!??  

        call idealloc(tempisp,ldummy) ! deallocating temp field
        call ralloc(smsirsp,nsp,'SMSOLV-irsp',ldummy)

        lodrv=.true.

c....   Time for symbolic factorization
        call etimef(tary1)
        if (pfr)             write(iow,2015) tary1-tary
        if(pfr.and.ior.lt.0) write(*  ,2015) tary1-tary
      end if

c---- at each call ----------------------------------------------------

c.... numerical factorisation
      call sdrv(neq,smsperm,smsperm1,jd,csrja,ad,trans(nw),trans(nw),
     +            nsp,smsisp,smsirsp,iesp,6,iflag)

c.... errors
      if (iflag.ne.0) then
        write(iow,*) '*** FEHLER IN SDRV (1) *** iflag=',iflag
        write(*  ,*) '*** FEHLER IN SDRV (1) *** iflag=',iflag
        call sderr(iflag,neq)
      end if

2014  format('   Begin    symbolic factorization',31x,'t=',0pf9.2)
2015  format('   Time for symbolic factorization',31x,'t=',0pf9.4)
      return
      end

      subroutine detkt2(neq,nneg,det1)
c----------------------------------------------------------------------
c     Purpose: calculate determinant of stiffness matrix
c              for sparse-matrix solver: value det1
c      Inputs:
c         smsirsp(iposd) - Factored diagonal terms of A-array
c         neq            - Length of A array
c
c      Outputs:
c        nneg - number of negative diagonal elements
c        det1 - value of determinant
c
c     W. Wagner BS UKA 04/09
c
c----------------------------------------------------------------------
      USE smpak
      implicit double precision(a-h,o-z)
      common /posd/ iposd  ! Position of D_ii in array smsirsp
      call detkt21(smsirsp(iposd),neq,nneg,det1)
      return
      end

      subroutine detkt21(ad,neq,nneg,det1)
c----------------------------------------------------------------------
c
c     Purpose: calculate determinant of stiffness matrix
c              called from sparse-matrix solver  value det1
c
c      Inputs:
c         ad(neq)  - Factored diagonal terms of A-array
c         neq      - Length of A array
c
c     Outputs:
c        nneg - number of negative diagonal elements
c        det1 - value of determinant
c
c     W. Wagner BS UKA 04/09
c
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension ad(*)
      data zero/0.0d0/
c...... calculate determinant
        do 10 i = 1,neq
          adi=1.d0/ad(i)
          if(abs(adi).gt.1.d-30) then
            if(adi.lt.zero) then
              xx = -adi
              nneg = nneg + 1
            else
              xx = adi
            endif
            det1 = det1 + log(xx)
          endif
   10   continue
      return
      end



c----------------------------------------------------------------------



      SUBROUTINE MD(N, IA,JA, MAX, V,L, HEAD,LAST,NEXT, MARK, FLAG)
c----------------------------------------------------------------------
c
c     Purpose:
c     MD FINDS A MINIMUM DEGREE ORDERING OF THE ROWS AND COLUMNS OF A
c     SYMMETRIC MATRIX M STORED IN (IA,JA,A) FORMAT.
c
c
c     ADDITIONAL PARAMETERS
c
c     MAX  - DECLARED DIMENSION OF THE ONE-DIMENSIONAL ARRAYS V AND L;
C            MAX MUST BE AT LEAST  N+2K,  WHERE K IS THE NUMBER OF
C            NONZEROES IN THE STRICT UPPER TRIANGLE OF M
C
C     V    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = MAX
C
C     L    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = MAX
C
C     HEAD - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C     LAST - INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN THE PERMUTATION
C            OF THE ROWS AND COLUMNS OF M CORRESPONDING TO THE MINIMUM
C            DEGREE ORDERING;  DIMENSION = N
C
C     NEXT - INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN THE INVERSE OF
C            THE PERMUTATION RETURNED IN LAST;  DIMENSION = N
C
C     MARK - INTEGER ONE-DIMENSIONAL WORK ARRAY (MAY BE THE SAME AS V);
C            DIMENSION = N
C
C     FLAG - INTEGER ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -
C             0      NO ERRORS DETECTED
C             11N+1  INSUFFICIENT STORAGE IN MD
C
C
C     DEFINITIONS OF INTERNAL PARAMETERS
C
C     ---------+--------------------------------------------------------
C     V(S)     \ VALUE FIELD OF LIST ENTRY
C     ---------+--------------------------------------------------------
C     L(S)     \ LINK FIELD OF LIST ENTRY  (0 => END OF LIST)
C     ---------+--------------------------------------------------------
C     L(VI)    \ POINTER TO ELEMENT LIST OF UNELIMINATED VERTEX VI
C     ---------+--------------------------------------------------------
C     L(EJ)    \ POINTER TO BOUNDARY LIST OF ACTIVE ELEMENT EJ
C     ---------+--------------------------------------------------------
C     HEAD(D)  \ VJ => VJ HEAD OF D-LIST D
C              \  0 => NO VERTEX IN D-LIST D
C
C
C              \                  VI UNELIMINATED VERTEX
C              \          VI IN EK           \       VI NOT IN EK
C     ---------+-----------------------------+--------------------------
C     NEXT(VI) \ UNDEFINED BUT NONNEGATIVE   \ VJ => VJ NEXT IN D-LIST
C             \                              \  0 => VI TAIL OF D-LIST
C     ---------+-----------------------------+--------------------------
C     LAST(VI) \ (NOT SET UNTIL MDP)         \ -D => VI HEAD OF D-LIST D
C              \-VK => COMPUTE DEGREE        \ VJ => VJ LAST IN D-LIST
C              \ EJ => VI PROTOTYPE OF EJ    \  0 => VI NOT IN ANY D-LIST
C              \  0 => DO NOT COMPUTE DEGREE \
C     ---------+-----------------------------+--------------------------
C     MARK(VI) \ MARK(VK)                    \ NONNEGATIVE TAG < MARK(VK)
C
C
C              \                   VI ELIMINATED VERTEX
C              \      EI ACTIVE ELEMENT      \           OTHERWISE
C     ---------+-----------------------------+--------------------------
C     NEXT(VI) \ -J => VI WAS J-TH VERTEX    \ -J => VI WAS J-TH VERTEX
C              \       TO BE ELIMINATED      \       TO BE ELIMINATED
C     --------+-----------------------------+---------------------------
C     LAST(VI) \  M => SIZE OF EI = M        \ UNDEFINED
C     --------+-----------------------------+---------------------------
C     MARK(VI) \ -M => OVERLAP COUNT OF EI   \ UNDEFINED
C              \       WITH EK = M           \
C              \ OTHERWISE NONNEGATIVE TAG   \
C              \       < MARK(VK)            \
C
C-----------------------------------------------------------------------
C
        INTEGER  IA(*), JA(*),  V(*), L(*),  HEAD(*), LAST(*), NEXT(*),
     *     MARK(*),  FLAG,  TAG, DMIN, VK,EK, TAIL
        EQUIVALENCE  (VK,EK)
C
C----INITIALIZATION
        TAG = 0
        CALL MDI
     *     (N, IA,JA, MAX,V,L, HEAD,LAST,NEXT, MARK,TAG, FLAG)
        IF (FLAG.NE.0)  RETURN
C
        K = 0
        DMIN = 1
C
C----WHILE  K < N  DO
   1    IF (K.GE.N)  GO TO 4
C
C------SEARCH FOR VERTEX OF MINIMUM DEGREE
   2      IF (HEAD(DMIN).GT.0)  GO TO 3
            DMIN = DMIN + 1
            GO TO 2
C
C------REMOVE VERTEX VK OF MINIMUM DEGREE FROM DEGREE LIST
   3      VK = HEAD(DMIN)
          HEAD(DMIN) = NEXT(VK)
          IF (HEAD(DMIN).GT.0)  LAST(HEAD(DMIN)) = -DMIN
C
C------NUMBER VERTEX VK, ADJUST TAG, AND TAG VK
          K = K+1
          NEXT(VK) = -K
          LAST(EK) = DMIN - 1
          TAG = TAG + LAST(EK)
          MARK(VK) = TAG
C
C------FORM ELEMENT EK FROM UNELIMINATED NEIGHBORS OF VK
          CALL MDM
     *       (VK,TAIL, V,L, LAST,NEXT, MARK)
C
C------PURGE INACTIVE ELEMENTS AND DO MASS ELIMINATION
          CALL MDP
     *       (K,EK,TAIL, V,L, HEAD,LAST,NEXT, MARK)
C
C------UPDATE DEGREES OF UNELIMINATED VERTICES IN EK
          CALL MDU
     *       (EK,DMIN, V,L, HEAD,LAST,NEXT, MARK)
C
          GO TO 1
C
C----GENERATE INVERSE PERMUTATION FROM PERMUTATION
   4    DO 5 K=1,N
          NEXT(K) = -NEXT(K)
   5      LAST(NEXT(K)) = K
C
        RETURN
        END
c
      SUBROUTINE MDI
     *     (N, IA,JA, MAX,V,L, HEAD,LAST,NEXT, MARK,TAG, FLAG)
        INTEGER  IA(*), JA(*),  V(*), L(*),  HEAD(*), LAST(*), NEXT(*),
     *     MARK(*), TAG,  FLAG,  SFS, VI,DVI, VJ
C
C----INITIALIZE DEGREES, ELEMENT LISTS, AND DEGREE LISTS
        DO 1 VI=1,N
          MARK(VI) = 1
          L(VI) = 0
   1      HEAD(VI) = 0
        SFS = N+1
C
C----CREATE NONZERO STRUCTURE
C----FOR EACH NONZERO ENTRY A(VI,VJ) IN STRICT UPPER TRIANGLE
        DO 3 VI=1,N
          JMIN = IA(VI)
          JMAX = IA(VI+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 3
          DO 2 J=JMIN,JMAX
            VJ = JA(J)
            IF (VI.GE.VJ)  GO TO 2
              IF (SFS.GE.MAX)  GO TO 101
C
C------ENTER VJ IN ELEMENT LIST FOR VI
              MARK(VI) = MARK(VI) + 1
              V(SFS) = VJ
              L(SFS) = L(VI)
              L(VI) = SFS
              SFS = SFS+1
C
C------ENTER VI IN ELEMENT LIST FOR VJ
              MARK(VJ) = MARK(VJ) + 1
              V(SFS) = VI
              L(SFS) = L(VJ)
              L(VJ) = SFS
              SFS = SFS+1
   2        CONTINUE
   3      CONTINUE
C
C----CREATE DEGREE LISTS AND INITIALIZE MARK VECTOR
        DO 4 VI=1,N
          DVI = MARK(VI)
          NEXT(VI) = HEAD(DVI)
          HEAD(DVI) = VI
          LAST(VI) = -DVI
          IF (NEXT(VI).GT.0)  LAST(NEXT(VI)) = VI
   4      MARK(VI) = TAG
C
        RETURN
C
C ** ERROR -- INSUFFICIENT STORAGE
 101    FLAG = 9*N + VI
        RETURN
        END
c
       SUBROUTINE MDM (VK,TAIL, V,L, LAST,NEXT, MARK)
        INTEGER  VK, TAIL,  V(*), L(*),  LAST(*), NEXT(*), MARK(*),
     *     TAG, S,LS,VS,ES, B,LB,VB, BLP,BLPMAX
        EQUIVALENCE  (VS, ES)
C
C----INITIALIZE TAG AND LIST OF UNELIMINATED NEIGHBORS
        TAG = MARK(VK)
        TAIL = VK
C
C----FOR EACH VERTEX/ELEMENT VS/ES IN ELEMENT LIST OF VK
        LS = L(VK)
   1    S = LS
        IF (S.EQ.0)  GO TO 5
          LS = L(S)
          VS = V(S)
          IF (NEXT(VS).LT.0)  GO TO 2
C
C------IF VS IS UNELIMINATED VERTEX, THEN TAG AND APPEND TO LIST OF
C------UNELIMINATED NEIGHBORS
            MARK(VS) = TAG
            L(TAIL) = S
            TAIL = S
            GO TO 4
C
C------IF ES IS ACTIVE ELEMENT, THEN ...
C--------FOR EACH VERTEX VB IN BOUNDARY LIST OF ELEMENT ES
   2        LB = L(ES)
            BLPMAX = LAST(ES)
            DO 3 BLP=1,BLPMAX
              B = LB
              LB = L(B)
              VB = V(B)
C
C----------IF VB IS UNTAGGED VERTEX, THEN TAG AND APPEND TO LIST OF
C----------UNELIMINATED NEIGHBORS
              IF (MARK(VB).GE.TAG)  GO TO 3
                MARK(VB) = TAG
                L(TAIL) = B
                TAIL = B
   3          CONTINUE
C
C--------MARK ES INACTIVE
            MARK(ES) = TAG
C
   4      GO TO 1
C
C----TERMINATE LIST OF UNELIMINATED NEIGHBORS
   5    L(TAIL) = 0
C
        RETURN
        END
c
        SUBROUTINE MDP (K,EK,TAIL, V,L, HEAD,LAST,NEXT, MARK)
        INTEGER  EK, TAIL,  V(*), L(*),  HEAD(*), LAST(*), NEXT(*),
     *     MARK(*), TAG, FREE, LI,VI,LVI,EVI, S,LS,ES, ILP,ILPMAX
C
C----INITIALIZE TAG
        TAG = MARK(EK)
C
C----FOR EACH VERTEX VI IN EK
        LI = EK
        ILPMAX = LAST(EK)
        IF (ILPMAX.LE.0)  GO TO 12
        DO 11 ILP=1,ILPMAX
          I = LI
          LI = L(I)
          VI = V(LI)
C
C------REMOVE VI FROM DEGREE LIST
          IF (LAST(VI).EQ.0)  GO TO 3
            IF (LAST(VI).GT.0)  GO TO 1
              HEAD(-LAST(VI)) = NEXT(VI)
              GO TO 2
   1          NEXT(LAST(VI)) = NEXT(VI)
   2        IF (NEXT(VI).GT.0)  LAST(NEXT(VI)) = LAST(VI)
C
C------REMOVE INACTIVE ITEMS FROM ELEMENT LIST OF VI
   3      LS = VI
   4      S = LS
          LS = L(S)
          IF (LS.EQ.0)  GO TO 6
            ES = V(LS)
            IF (MARK(ES).LT.TAG)  GO TO 5
              FREE = LS
              L(S) = L(LS)
              LS = S
   5        GO TO 4
C
C------IF VI IS INTERIOR VERTEX, THEN REMOVE FROM LIST AND ELIMINATE
   6      LVI = L(VI)
          IF (LVI.NE.0)  GO TO 7
            L(I) = L(LI)
            LI = I
C
            K = K+1
            NEXT(VI) = -K
            LAST(EK) = LAST(EK) - 1
            GO TO 11
C
C------ELSE ...
C--------CLASSIFY VERTEX VI
   7        IF (L(LVI).NE.0)  GO TO 9
              EVI = V(LVI)
              IF (NEXT(EVI).GE.0)  GO TO 9
                IF (MARK(EVI).LT.0)  GO TO 8
C
C----------IF VI IS PROTOTYPE VERTEX, THEN MARK AS SUCH, INITIALIZE
C----------OVERLAP COUNT FOR CORRESPONDING ELEMENT, AND MOVE VI TO END
C----------OF BOUNDARY LIST
                  LAST(VI) = EVI
                  MARK(EVI) = -1
                  L(TAIL) = LI
                  TAIL = LI
                  L(I) = L(LI)
                  LI = I
                  GO TO 10
C
C----------ELSE IF VI IS DUPLICATE VERTEX, THEN MARK AS SUCH AND ADJUST
C----------OVERLAP COUNT FOR CORRESPONDING ELEMENT
   8              LAST(VI) = 0
                  MARK(EVI) = MARK(EVI) - 1
                  GO TO 10
C
C----------ELSE MARK VI TO COMPUTE DEGREE
   9              LAST(VI) = -EK
C
C--------INSERT EK IN ELEMENT LIST OF VI
  10        V(FREE) = EK
            L(FREE) = L(VI)
            L(VI) = FREE
  11      CONTINUE
C
C----TERMINATE BOUNDARY LIST
  12    L(TAIL) = 0
C
        RETURN
        END
C


       SUBROUTINE MDU (EK,DMIN, V,L, HEAD,LAST,NEXT, MARK)



        INTEGER  EK, DMIN,  V(*), L(*),  HEAD(*), LAST(*), NEXT(*),
     *     MARK(*),  TAG, VI,EVI,DVI, S,VS,ES, B,VB, ILP,ILPMAX,
     *     BLP,BLPMAX
        EQUIVALENCE  (VS, ES)
C
C----INITIALIZE TAG
        TAG = MARK(EK) - LAST(EK)
C
C----FOR EACH VERTEX VI IN EK
        I = EK
        ILPMAX = LAST(EK)
        IF (ILPMAX.LE.0)  GO TO 11
        DO 10 ILP=1,ILPMAX
          I = L(I)
          VI = V(I)
          IF (LAST(VI))  1, 10, 8
C
C------IF VI NEITHER PROTOTYPE NOR DUPLICATE VERTEX, THEN MERGE ELEMENTS
C------TO COMPUTE DEGREE
   1        TAG = TAG + 1
            DVI = LAST(EK)
C
C--------FOR EACH VERTEX/ELEMENT VS/ES IN ELEMENT LIST OF VI
            S = L(VI)
   2        S = L(S)
            IF (S.EQ.0)  GO TO 9
              VS = V(S)
              IF (NEXT(VS).LT.0)  GO TO 3
C
C----------IF VS IS UNELIMINATED VERTEX, THEN TAG AND ADJUST DEGREE
                MARK(VS) = TAG
                DVI = DVI + 1
                GO TO 5
C
C----------IF ES IS ACTIVE ELEMENT, THEN EXPAND
C------------CHECK FOR OUTMATCHED VERTEX
   3            IF (MARK(ES).LT.0)  GO TO 6
C
C------------FOR EACH VERTEX VB IN ES
                B = ES
                BLPMAX = LAST(ES)
                DO 4 BLP=1,BLPMAX
                  B = L(B)
                  VB = V(B)
C
C--------------IF VB IS UNTAGGED, THEN TAG AND ADJUST DEGREE
                  IF (MARK(VB).GE.TAG)  GO TO 4
                    MARK(VB) = TAG
                    DVI = DVI + 1
   4              CONTINUE
C
   5          GO TO 2
C
C------ELSE IF VI IS OUTMATCHED VERTEX, THEN ADJUST OVERLAPS BUT DO NOT
C------COMPUTE DEGREE
   6        LAST(VI) = 0
            MARK(ES) = MARK(ES) - 1
   7        S = L(S)
            IF (S.EQ.0)  GO TO 10
              ES = V(S)
              IF (MARK(ES).LT.0)  MARK(ES) = MARK(ES) - 1
              GO TO 7
C
C------ELSE IF VI IS PROTOTYPE VERTEX, THEN CALCULATE DEGREE BY
C------INCLUSION/EXCLUSION AND RESET OVERLAP COUNT
   8        EVI = LAST(VI)
            DVI = LAST(EK) + LAST(EVI) + MARK(EVI)
            MARK(EVI) = 0
C
C------INSERT VI IN APPROPRIATE DEGREE LIST
   9      NEXT(VI) = HEAD(DVI)
          HEAD(DVI) = VI
          LAST(VI) = -DVI
          IF (NEXT(VI).GT.0)  LAST(NEXT(VI)) = VI
          IF (DVI.LT.DMIN)  DMIN = DVI
C
  10      CONTINUE
C
  11    RETURN
        END
C-----------------------------------------------------------
      SUBROUTINE ODRV(N, IA,JA,A, P,IP, NSP,ISP, PATH, FLAG)
C
C     DESCRIPTION
C
C     ODRV FINDS A MINIMUM DEGREE ORDERING OF THE ROWS AND COLUMNS OF A
C     SYMMETRIC MATRIX M STORED IN (IA,JA,A) FORMAT (SEE BELOW).  FOR THE
C     REORDERED MATRIX, THE WORK AND STORAGE REQUIRED TO PERFORM GAUSSIAN
C     ELIMINATION IS (USUALLY) SIGNIFICANTLY LESS.
C
C     IF ONLY THE NONZERO ENTRIES IN THE UPPER TRIANGLE OF M ARE BEING
C     STORED, THEN ODRV SYMMETRICALLY REORDERS (IA,JA,A), (OPTIONALLY)
C     WITH THE DIAGONAL ENTRIES PLACED FIRST IN EACH ROW.  THIS IS TO
C     ENSURE THAT IF M(I,J) WILL BE IN THE UPPER TRIANGLE OF M WITH
C     RESPECT TO THE NEW ORDERING, THEN M(I,J) IS STORED IN ROW I (AND
C     THUS M(J,I) IS NOT STORED);  WHEREAS IF M(I,J) WILL BE IN THE
C     STRICT LOWER TRIANGLE OF M, THEN M(J,I) IS STORED IN ROW J (AND
C     THUS M(I,J) IS NOT STORED).
C
C
C     STORAGE OF SPARSE MATRICES
C
C     THE NONZERO ENTRIES OF THE MATRIX M ARE STORED ROW-BY-ROW IN THE
C     ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO ENTRIES IN EACH ROW,
C     WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY LIES.  THESE COLUMN
C     INDICES ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN
C     JA(K) = J.  TO IDENTIFY THE INDIVIDUAL ROWS, WE NEED TO KNOW WHERE
C     EACH ROW STARTS.  THESE ROW POINTERS ARE STORED IN THE ARRAY IA;
C     I.E., IF M(I,J) IS THE FIRST NONZERO ENTRY (STORED) IN THE I-TH ROW
C     AND  A(K) = M(I,J),  THEN  IA(I) = K.  MOREOVER, IA(N+1) POINTS TO
C     THE FIRST LOCATION FOLLOWING THE LAST ELEMENT IN THE LAST ROW.
C     THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS  IA(I+1) - IA(I),
C     THE NONZERO ENTRIES IN THE I-TH ROW ARE STORED CONSECUTIVELY IN
C
C            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1),
C
C     AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN
C
C            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1).
C
C     SINCE THE COEFFICIENT MATRIX IS SYMMETRIC, ONLY THE NONZERO ENTRIES
C     IN THE UPPER TRIANGLE NEED BE STORED.  FOR EXAMPLE, THE MATRIX
C
C             ( 1  0  2  3  0 )
C             ( 0  4  0  0  0 )
C         M = ( 2  0  5  6  0 )
C             ( 3  0  6  7  8 )
C             ( 0  0  0  8  9 )
C
C     COULD BE STORED AS
C
C            \ 1  2  3  4  5  6  7  8  9 10 11 12 13
C         ---+--------------------------------------
C         IA \ 1  4  5  8 12 14
C         JA \ 1  3  4  2  1  3  4  1  3  4  5  4  5
C          A \ 1  2  3  4  2  5  6  3  6  7  8  8  9
C
C     OR (SYMMETRICALLY) AS
C
C            \ 1  2  3  4  5  6  7  8  9
C         ---+--------------------------
C         IA \ 1  4  5  7  9 10
C         JA \ 1  3  4  2  3  4  4  5  5
C          A \ 1  2  3  4  5  6  7  8  9          .
C
C
C     PARAMETERS
C
C     N    - ORDER OF THE MATRIX
C
C     IA   - INTEGER ONE-DIMENSIONAL ARRAY CONTAINING POINTERS TO DELIMIT
C            ROWS IN JA AND A;  DIMENSION = N+1
C
C     JA   - INTEGER ONE-DIMENSIONAL ARRAY CONTAINING THE COLUMN INDICES
C            CORRESPONDING TO THE ELEMENTS OF A;  DIMENSION = NUMBER OF
C            NONZERO ENTRIES IN (THE UPPER TRIANGLE OF) M
C
C     A    - REAL ONE-DIMENSIONAL ARRAY CONTAINING THE NONZERO ENTRIES IN
C            (THE UPPER TRIANGLE OF) M, STORED BY ROWS;  DIMENSION =
C           NUMBER OF NONZERO ENTRIES IN (THE UPPER TRIANGLE OF) M
C
C     P    - INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN THE PERMUTATION
C            OF THE ROWS AND COLUMNS OF M CORRESPONDING TO THE MINIMUM
C            DEGREE ORDERING;  DIMENSION = N
C
C     IP   - INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN THE INVERSE OF
C            THE PERMUTATION RETURNED IN P;  DIMENSION = N
C
C     NSP  - DECLARED DIMENSION OF THE ONE-DIMENSIONAL ARRAY ISP;  NSP
C            MUST BE AT LEAST  3N+4K,  WHERE K IS THE NUMBER OF NONZEROES
C            IN THE STRICT UPPER TRIANGLE OF M
C
C     ISP  - INTEGER ONE-DIMENSIONAL ARRAY USED FOR WORKING STORAGE;
C            DIMENSION = NSP
C
C     PATH - INTEGER PATH SPECIFICATION;  VALUES AND THEIR MEANINGS ARE -
C             1  FIND MINIMUM DEGREE ORDERING ONLY
C             2  FIND MINIMUM DEGREE ORDERING AND REORDER SYMMETRICALLY
C                  STORED MATRIX (USED WHEN ONLY THE NONZERO ENTRIES IN
C                  THE UPPER TRIANGLE OF M ARE BEING STORED)
C             3  REORDER SYMMETRICALLY STORED MATRIX AS SPECIFIED BY
C                  INPUT PERMUTATION (USED WHEN AN ORDERING HAS ALREADY
C                  BEEN DETERMINED AND ONLY THE NONZERO ENTRIES IN THE
C                  UPPER TRIANGLE OF M ARE BEING STORED)
C             4  SAME AS 2 BUT PUT DIAGONAL ENTRIES AT START OF EACH ROW
C             5  SAME AS 3 BUT PUT DIAGONAL ENTRIES AT START OF EACH ROW
C
C     FLAG - INTEGER ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -
C               0    NO ERRORS DETECTED
C              9N+K  INSUFFICIENT STORAGE IN MD
C             10N+1  INSUFFICIENT STORAGE IN ODRV
C             11N+1  ILLEGAL PATH SPECIFICATION
C
C
C     CONVERSION FROM REAL TO DOUBLE PRECISION
C
C     CHANGE THE REAL DECLARATIONS IN ODRV AND SRO TO DOUBLE PRECISION
C     DECLARATIONS.
C
C-----------------------------------------------------------------------
C
      INTEGER  IA(*), JA(*),  P(*), IP(*),  ISP(*),  PATH,  FLAG,
     *     V, L, HEAD,  TMP, Q
      DOUBLE PRECISION  A(*)
      LOGICAL  DFLAG
C
C---- INITIALIZE ERROR FLAG AND VALIDATE PATH SPECIFICATION
        FLAG = 0
        IF (PATH.LT.1 .OR. 5.LT.PATH)  GO TO 111
C
C---- ALLOCATE STORAGE AND FIND MINIMUM DEGREE ORDERING
        IF ((PATH-1) * (PATH-2) * (PATH-4) .NE. 0)  GO TO 1
          MAX = (NSP-N)/2
          V    = 1
          L    = V     +  MAX
          HEAD = L     +  MAX
cww          NEXT = HEAD  +  N
          IF (MAX.LT.N)  GO TO 110
C
          CALL MD
     *       (N, IA,JA, MAX,ISP(V),ISP(L), ISP(HEAD),P,IP, ISP(V), FLAG)
          IF (FLAG.NE.0)  GO TO 100
C
C---- ALLOCATE STORAGE AND SYMMETRICALLY REORDER MATRIX
   1    IF ((PATH-2) * (PATH-3) * (PATH-4) * (PATH-5) .NE. 0)  GO TO 2
          TMP = (NSP+1) -      N
          Q   = TMP     - (IA(N+1)-1)
          IF (Q.LT.1)  GO TO 110
C
          DFLAG = PATH.EQ.4 .OR. PATH.EQ.5
          CALL SRO
     *       (N,  IP,  IA, JA, A,  ISP(TMP),  ISP(Q),  DFLAG)
C
   2    RETURN
C
C **  ERROR -- ERROR DETECTED IN MD
 100    RETURN
C **  ERROR -- INSUFFICIENT STORAGE
 110    FLAG = 10*N + 1
        RETURN
C **  ERROR -- ILLEGAL PATH SPECIFIED
 111    FLAG = 11*N + 1
        RETURN
        END
c
        SUBROUTINE SDRV
     *     (N, P,IP, IA,JA,A, B, Z, NSP,ISP,RSP,ESP, PATH, FLAG)
C
C     DESCRIPTION
C
C     SDRV SOLVES SPARSE SYMMETRIC POSITIVE DEFINITE SYSTEMS OF LINEAR
C     EQUATIONS.  THE SOLUTION PROCESS IS DIVIDED INTO THREE STAGES --
C
C      SSF - THE COEFFICIENT MATRIX M IS FACTORED SYMBOLICALLY TO
C            DETERMINE WHERE FILLIN WILL OCCUR DURING THE NUMERIC
C            FACTORIZATION.
C
C      SNF - M IS FACTORED NUMERICALLY INTO THE PRODUCT UT-D-U, WHERE
C            D IS DIAGONAL AND U IS UNIT UPPER TRIANGULAR.
C
C      SNS - THE LINEAR SYSTEM  MX = B  IS SOLVED USING THE UT-D-U
C            FACTORIZATION FROM SNF.
C
C     FOR SEVERAL SYSTEMS WITH THE SAME COEFFICIENT MATRIX, SSF AND SNF
C     NEED BE DONE ONLY ONCE (FOR THE FIRST SYSTEM);  THEN SNS IS DONE
C     ONCE FOR EACH ADDITIONAL RIGHT-HAND SIDE.  FOR SEVERAL SYSTEMS
C     WHOSE COEFFICIENT MATRICES HAVE THE SAME NONZERO STRUCTURE, SSF
C     NEED BE DONE ONLY ONCE (FOR THE FIRST SYSTEM);  THEN SNF AND SNS
C     ARE DONE ONCE FOR EACH ADDITIONAL SYSTEM.
C
C
C     STORAGE OF SPARSE MATRICES
C
C     THE NONZERO ENTRIES OF THE MATRIX M ARE STORED ROW-BY-ROW IN THE
C     ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO ENTRIES IN EACH ROW,
C     WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY LIES.  THESE COLUMN
C     INDICES ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN
C     JA(K) = J.  TO IDENTIFY THE INDIVIDUAL ROWS, WE NEED TO KNOW WHERE
C     EACH ROW STARTS.  THESE ROW POINTERS ARE STORED IN THE ARRAY IA;
C     I.E., IF M(I,J) IS THE FIRST NONZERO ENTRY (STORED) IN THE I-TH ROW
C     AND  A(K) = M(I,J),  THEN  IA(I) = K.  MOREOVER, IA(N+1) POINTS TO
C     THE FIRST LOCATION FOLLOWING THE LAST ELEMENT IN THE LAST ROW.
C     THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS  IA(I+1) - IA(I),
C     THE NONZERO ENTRIES IN THE I-TH ROW ARE STORED CONSECUTIVELY IN
C
C            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1),
C
C     AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN
C
C            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1).
C
C     SINCE THE COEFFICIENT MATRIX IS SYMMETRIC, ONLY THE NONZERO ENTRIES
C     IN THE UPPER TRIANGLE NEED BE STORED, FOR EXAMPLE, THE MATRIX
C
C             ( 1  0  2  3  0 )
C             ( 0  4  0  0  0 )
C         M = ( 2  0  5  6  0 )
C             ( 3  0  6  7  8 )
C             ( 0  0  0  8  9 )
C
C     COULD BE STORED AS
C
C            \ 1  2  3  4  5  6  7  8  9 10 11 12 13
C         ---+--------------------------------------
C         IA \ 1  4  5  8 12 14
C         JA \ 1  3  4  2  1  3  4  1  3  4  5  4  5
C          A \ 1  2  3  4  2  5  6  3  6  7  8  8  9
C
C     OR (SYMMETRICALLY) AS
C
C            \ 1  2  3  4  5  6  7  8  9
C         ---+--------------------------
C         IA \ 1  4  5  7  9 10
C         JA \ 1  3  4  2  3  4  4  5  5
C          A \ 1  2  3  4  5  6  7  8  9          .
C
C
C     REORDERING THE ROWS AND COLUMNS OF M
C
C     A SYMMETRIC PERMUTATION OF THE ROWS AND COLUMNS OF THE COEFFICIENT
C     MATRIX M (E.G., WHICH REDUCES FILLIN OR ENHANCES NUMERICAL
C     STABILITY) MUST BE SPECIFIED.  THE SOLUTION Z IS RETURNED IN THE
C     ORIGINAL ORDER.
C
C     TO SPECIFY THE TRIVIAL ORDERING (I.E., THE IDENTITY PERMUTATION),
C     SET  P(I) = IP(I) = I,  I=1,...,N.  IN THIS CASE, P AND IP CAN BE
C     THE SAME ARRAY.
C
C     IF A NONTRIVIAL ORDERING (I.E., NOT THE IDENTITY PERMUTATION) IS
C     SPECIFIED AND M IS STORED SYMMETRICALLY (I.E., NOT BOTH M(I,J) AND
C     M(J,I) ARE STORED FOR I NE J), THEN ODRV SHOULD BE CALLED (WITH
C     PATH = 3 OR 5) TO SYMMETRICALLY REORDER (IA,JA,A) BEFORE CALLING
C     SDRV.  THIS IS TO ENSURE THAT IF M(I,J) WILL BE IN THE UPPER
C     TRIANGLE OF M WITH RESPECT TO THE NEW ORDERING, THEN M(I,J) IS
C     STORED IN ROW I (AND THUS M(J,I) IS NOT STORED);  WHEREAS IF M(I,J)
C     WILL BE IN THE STRICT LOWER TRIANGLE OF M, THEN M(J,I) IS STORED IN
C     ROW J (AND THUS M(I,J) IS NOT STORED).
C
C
C     PARAMETERS
C
C     N    - NUMBER OF VARIABLES/EQUATIONS
C
C     P    - INTEGER ONE-DIMENSIONAL ARRAY SPECIFYING A PERMUTATION OF
C            THE ROWS AND COLUMNS OF M;  DIMENSION = N
C
C     IP   - INTEGER ONE-DIMENSIONAL ARRAY CONTAINING THE INVERSE OF THE
C            PERMUTATION SPECIFIED IN P;  I.E., IP(P(I)) = I, I=1,...,N;
C            DIMENSION = N
C
C     IA   - INTEGER ONE-DIMENSIONAL ARRAY CONTAINING POINTERS TO DELIMIT
C            ROWS IN JA AND A;  DIMENSION = N+1
C
C     JA   - INTEGER ONE-DIMENSIONAL ARRAY CONTAINING THE COLUMN INDICES
C            CORRESPONDING TO THE ELEMENTS OF A;  DIMENSION = NUMBER OF
C            NONZERO ENTRIES IN M STORED
C
C     A    - REAL ONE-DIMENSIONAL ARRAY CONTAINING THE NONZERO ENTRIES IN
C           THE COEFFICIENT MATRIX M, STORED BY ROWS;  DIMENSION =
C           NUMBER OF NONZERO ENTRIES IN M STORED
C
C     B    - REAL ONE-DIMENSIONAL ARRAY CONTAINING THE RIGHT-HAND SIDE B;
C            B AND Z CAN BE THE SAME ARRAY;  DIMENSION = N
C
C     Z    - REAL ONE-DIMENSIONAL ARRAY CONTAINING THE SOLUTION X;  Z AND
C            B CAN BE THE SAME ARRAY;  DIMENSION = N
C
C     NSP  - DECLARED DIMENSION OF THE ONE-DIMENSIONAL ARRAYS ISP AND
C            RSP;  NSP MUST BE (SUBSTANTIALLY) LARGER THAN  3N+2K,  WHERE
C            K = NUMBER OF NONZERO ENTRIES IN THE UPPER TRIANGLE OF M
C
C     ISP  - INTEGER ONE-DIMENSIONAL ARRAY USED FOR WORKING STORAGE;  ISP
C            AND RSP SHOULD BE EQUIVALENCED;  DIMENSION = NSP
C
C     RSP  - REAL ONE-DIMENSIONAL ARRAY USED FOR WORKING STORAGE;  RSP
C            AND ISP SHOULD BE EQUIVALENCED;  DIMENSION = NSP
C
C     ESP  - INTEGER VARIABLE;  IF SUFFICIENT STORAGE WAS AVAILABLE TO
C            PERFORM THE SYMBOLIC FACTORIZATION (SSF), THEN ESP IS SET TO
C            THE AMOUNT OF EXCESS STORAGE PROVIDED (NEGATIVE IF
C            INSUFFICIENT STORAGE WAS AVAILABLE TO PERFORM THE NUMERIC
C            FACTORIZATION (SNF))
C
C     PATH - INTEGER PATH SPECIFICATION;  VALUES AND THEIR MEANINGS ARE -
C             1  PERFORM SSF, SNF, AND SNS
C             2  PERFORM SNF AND SNS (ISP/RSP IS ASSUMED TO HAVE BEEN
C                  SET UP IN AN EARLIER CALL TO SDRV (FOR SSF))
C             3  PERFORM SNS ONLY (ISP/RSP IS ASSUMED TO HAVE BEEN SET
C                  UP IN AN EARLIER CALL TO SDRV (FOR SSF AND SNF))
C             4  PERFORM SSF
C             5  PERFORM SSF AND SNF
C             6  PERFORM SNF ONLY (ISP/RSP IS ASSUMED TO HAVE BEEN SET
C                  UP IN AN EARLIER CALL TO SDRV (FOR SSF))
C
C     FLAG - INTEGER ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -
C               0     NO ERRORS DETECTED
C              2N+K   DUPLICATE ENTRY IN A         --  ROW = K
C              6N+K   INSUFFICIENT STORAGE IN SSF  --  ROW = K
C              7N+1   INSUFFICIENT STORAGE IN SNF
C              8N+K   ZERO PIVOT                   --  ROW = K
C             10N+1   INSUFFICIENT STORAGE IN SDRV
C             11N+1   ILLEGAL PATH SPECIFICATION
C
C
C     CONVERSION FROM REAL TO DOUBLE PRECISION
C
C     CHANGE THE REAL DECLARATIONS IN SDRV, SNF, AND SNS TO DOUBLE
C     PRECISION DECLARATIONS;  AND CHANGE THE VALUE IN THE DATA STATEMENT
C     FOR THE INTEGER VARIABLE RATIO (IN SDRV) FROM 1 TO 2.
C
C-----------------------------------------------------------------------
C

        common /posd/ iposd ! ww Position of Diagonals in RSP

        INTEGER  P(*), IP(*),  IA(*), JA(*),  ISP(*), ESP,  PATH,  FLAG,
     *     RATIO,  Q, MARK, D, U, TMP,  UMAX, iposd
        DOUBLE PRECISION  A(*),  B(*),  Z(*),  RSP(*)
        DATA  RATIO/2/
C
C---- VALIDATE PATH SPECIFICATION
        IF (PATH.LT.1 .OR. 6.LT.PATH)  GO TO 111
c.kneb>
c        open(20,file='stiff.tem')
c
c        close(20)
c.kneb<
c.thr>
        flag=0
c.thr<
C
C---- ALLOCATE STORAGE AND FACTOR M SYMBOLICALLY TO DETERMINE FILL-IN
        IJU   = 1
        IU    = IJU     +  N
        JL    = IU      + N+1
        JU    = JL      +  N
        Q     = (NSP+1) -  N
        MARK  = Q       -  N
        JUMAX = MARK    - JU
C
        IF ((PATH-1) * (PATH-4) * (PATH-5) .NE. 0)  GO TO 1
          IF (JUMAX.LE.0)  GO TO 110
          CALL SSF
     *       (N,  P, IP,  IA, JA,  ISP(IJU), ISP(JU), ISP(IU), JUMAX,
     *        ISP(Q), ISP(MARK), ISP(JL), FLAG)
          IF (FLAG.NE.0)  GO TO 100
C
C---- ALLOCATE STORAGE AND FACTOR M NUMERICALLY
   1    IL   = JU      + ISP(IJU+(N-1))
        TMP  = ((IL-1)+(RATIO-1)) / RATIO  +  1
        D    = TMP     + N
        iposd= D
        U    = D       + N
        UMAX = (NSP+1) - U
        ESP  = UMAX    - (ISP(IU+N)-1)
        IF ((PATH-1) * (PATH-2) * (PATH-5) * (PATH-6) .NE. 0)  GO TO 2
          IF (UMAX.LE.0)  GO TO 110
          CALL SNF
     *       (N,  P, IP,  IA, JA, A,
     *        RSP(D),  ISP(IJU), ISP(JU), ISP(IU), RSP(U), UMAX,
     *        ISP(IL),  ISP(JL),  FLAG)
          IF (FLAG.NE.0)  GO TO 100
C
C---- SOLVE SYSTEM OF LINEAR EQUATIONS  MX = B
   2    IF ((PATH-1) * (PATH-2) * (PATH-3) .NE. 0)  GO TO 3
          IF (UMAX.LE.0)  GO TO 110
          CALL SNS
     *       (N,  P,  RSP(D), ISP(IJU), ISP(JU), ISP(IU), RSP(U),  Z, B,
     *        RSP(TMP))
C
   3    RETURN
C
C **  ERROR -- ERROR DETECTED IN SSF, SNF, OR SNS
 100    RETURN
C **  ERROR -- INSUFFICIENT STORAGE
 110    FLAG = 10*N + 1
        write(*,*) 'INSUFFICIENT STORAGE IN SDRV'
        RETURN
C **  ERROR -- ILLEGAL PATH SPECIFICATION
 111    FLAG = 11*N + 1
        write(*,*) 'ILLEGAL PATH SPECIFICATION'
        RETURN
        END
C-----------------------------------------------------------
       SUBROUTINE SNF
     *     (N, P,IP, IA,JA,A, D, IJU,JU,IU,U,UMAX, IL, JL, FLAG)
C
C     ADDITIONAL PARAMETERS
C
C     IL    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C     JL    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C
C     DEFINITIONS OF INTERNAL PARAMETERS (DURING K-TH STAGE OF ELIMINATION)
C
C     (D(I),I=K,N) CONTAINS THE K-TH ROW OF U (EXPANDED)
C
C     IL(I) POINTS TO THE FIRST NONZERO ELEMENT IN COLUMNS K,...,N OF
C      ROW I OF U
C
C     JL CONTAINS LISTS OF ROWS TO BE ADDED TO UNELIMINATED ROWS --
C      I GE K => JL(I) IS THE FIRST ROW TO BE ADDED TO ROW I
C      I LT K => JL(I) IS THE ROW FOLLOWING ROW I IN SOME LIST OF ROWS
C      IN EITHER CASE, JL(I) = 0 INDICATES THE END OF A LIST
C
C-----------------------------------------------------------------------
C
      USE conv
      USE iofile
        INTEGER  P(*), IP(*),  IA(*), JA(*),  IJU(*), JU(*), IU(*),
     *     UMAX,  IL(*),  JL(*),  FLAG,  VJ
        DOUBLE PRECISION  A(*),  D(*), U(*),  DK, UKIDI
C
C---- CHECK FOR SUFFICIENT STORAGE FOR U
cww     IF (IU(N+1)-1 .GT. UMAX)  GO TO 107
        IF (IU(N+1)-1 .GT. UMAX)  then
C         ERROR -- INSUFFICIENT STORAGE FOR U
          FLAG = 7*N + 1
          write(*,*) '    INSUFFICIENT STORAGE IN SNF'
          RETURN
        end if
C
C---- INITIALIZATION
      call perform(0,N,1,25)
      DO 1 K=1,N
        call perform(K,N,2,25)
        D(K) = 0
   1    JL(K) = 0
      call perform(N,N,3,25) ! not really correct

C     count number of zero and neg. diagonals
          nneg=0

C
C---- FOR EACH ROW K
        DO 11 K=1,N
C
C------ INITIALIZE K-TH ROW WITH ELEMENTS NONZERO IN ROW P(K) OF M
          JMIN = IA(P(K))
cww3      JMIN = IA(P(K))
          JMAX = IA(P(K)+1) - 1
          IF (JMIN.GT.JMAX) GO TO 5
          DO 4 J=JMIN,JMAX
            VJ = IP(JA(J))
            IF (K.LE.VJ)  D(VJ) = A(J)
   4        CONTINUE
C
C------ MODIFY K-TH ROW BY ADDING IN THOSE ROWS I WITH U(I,K) NE 0
C------ FOR EACH ROW I TO BE ADDED IN
   5      DK = D(K)
          I = JL(K)
   6      IF (I.EQ.0)  GO TO 9
          NEXTI = JL(I)
C
C-------COMPUTE MULTIPLIER AND UPDATE DIAGONAL ELEMENT
          ILI = IL(I)
          UKIDI = - U(ILI) * D(I)
          DK = DK + UKIDI * U(ILI)
          U(ILI) = UKIDI
C
C-------ADD MULTIPLE OF ROW I TO K-TH ROW ...
          JMIN = ILI     + 1
          JMAX = IU(I+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 8
          MU = IJU(I) - IU(I)
          DO 7 J=JMIN,JMAX
   7            D(JU(MU+J)) = D(JU(MU+J)) + UKIDI * U(J)
C
C-------AND ADD I TO ROW LIST FOR NEXT NONZERO ENTRY
          IL(I) = JMIN
          J = JU(MU+JMIN)
          JL(I) = JL(J)
          JL(J) = I
C
   8      I = NEXTI
          GO TO 6
C
C------ CHECK FOR ZERO PIVOT AND SAVE DIAGONAL ELEMENT
cww   9   IF (DK.EQ.0)  GO TO 108
c
cww>
   9      IF (DK.EQ.0.d0)  then
C           ERROR -- ZERO PIVOT
cww         FLAG = 8*N + K
cww         write(*,*) '   ZERO PIVOT in equation', P(K)
            ii=p(k)
            call prtdii1(ii,2)
            DK=1.d0 ! set to one, ok if no further entry ww
          end if
c
          if (DK.lt.0.d0)  then
c           negative diagonal
            ii=p(k)
            call prtdii1(ii,1)
            nneg=nneg+1
          end if
cww<
          D(K) = 1 / DK
C
C------ SAVE NONZERO ENTRIES IN K-TH ROW OF U ...
          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 11
          MU = IJU(K) - JMIN
          DO 10 J=JMIN,JMAX
            JUMUJ = JU(MU+J)
            U(J) = D(JUMUJ)
  10      D(JUMUJ) = 0
C
C------ AND ADD K TO ROW LIST FOR FIRST NONZERO ENTRY IN K-TH ROW
          IL(K) = JMIN
          I = JU(MU+JMIN)
          JL(K) = JL(I)
          JL(I) = K
  11      CONTINUE
C
          FLAG = 0

          RETURN
C
cwwC **    ERROR -- INSUFFICIENT STORAGE FOR U
cww 107      FLAG = 7*N + 1
cww          RETURN
cwwC **    ERROR -- ZERO PIVOT
cww 108      FLAG = 8*N + K
cww          RETURN
cww2005    format('   FEAP sm-solve equations: '/)
cww2006  format('  current eq. ',i6,'  max. no. of eqs. ',i6)
        END
c
      SUBROUTINE SNS (N, P, D, IJU,JU,IU,U, Z, B, TMP)
        INTEGER  P(*),  IJU(*), JU(*), IU(*)
        DOUBLE PRECISION  D(*), U(*),  Z(*), B(*),  TMP(*),  TMPK, SUM
C
C  ADDITIONAL PARAMETERS
C
C    TMP   - REAL ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C-----------------------------------------------------------------------
C
C----SET TMP TO PERMUTED B
        DO 1 K=1,N
   1      TMP(K) = B(P(K))
C
C----SOLVE  UT D Y = B  BY FORWARD SUBSTITUTION
      call perform(0,N,1,26)
        DO 3 K=1,N
          TMPK = TMP(K)
          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 3
          MU = IJU(K) - JMIN
          DO 2 J=JMIN,JMAX
   2        TMP(JU(MU+J)) = TMP(JU(MU+J)) + U(J) * TMPK
          call perform(K,N,2,26)
   3      TMP(K) = TMPK * D(K)
        call perform(N,N,3,26)
C
C----SOLVE  U X = Y  BY BACK SUBSTITUTION
        K = N
cww   call perform(0,N,1,27)
        DO 6 I=1,N
          SUM = TMP(K)
          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 5
          MU = IJU(K) - JMIN
          DO 4 J=JMIN,JMAX
   4        SUM = SUM + U(J) * TMP(JU(MU+J))
   5      TMP(K) = SUM
          Z(P(K)) = SUM
cww     call perform(I,N,1,27)
   6      K = K-1
cww     call perform(N,N,1,27)
C
        RETURN
        END
c
       SUBROUTINE SRO
     *     (N, IP, IA,JA,A, Q, R, DFLAG)
C
C  DESCRIPTION
C
C    THE NONZERO ENTRIES OF THE MATRIX M ARE ASSUMED TO BE STORED
C    SYMMETRICALLY IN (IA,JA,A) FORMAT (I.E., NOT BOTH M(I,J) AND M(J,I)
C    ARE STORED IF I NE J).
C
C    SRO DOES NOT REARRANGE THE ORDER OF THE ROWS, BUT DOES MOVE
C    NONZEROES FROM ONE ROW TO ANOTHER TO ENSURE THAT IF M(I,J) WILL BE
C    IN THE UPPER TRIANGLE OF M WITH RESPECT TO THE NEW ORDERING, THEN
C    M(I,J) IS STORED IN ROW I (AND THUS M(J,I) IS NOT STORED);  WHEREAS
C    IF M(I,J) WILL BE IN THE STRICT LOWER TRIANGLE OF M, THEN M(J,I) IS
C    STORED IN ROW J (AND THUS M(I,J) IS NOT STORED).
C
C
C  ADDITIONAL PARAMETERS
C
C    Q     - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C    R     - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = NUMBER OF
C            NONZERO ENTRIES IN THE UPPER TRIANGLE OF M
C
C    DFLAG - LOGICAL VARIABLE;  IF DFLAG = .TRUE., THEN STORE NONZERO
C            DIAGONAL ELEMENTS AT THE BEGINNING OF THE ROW
C
C-----------------------------------------------------------------------
C
        INTEGER  IP(*),  IA(*), JA(*),  Q(*), R(*)
        DOUBLE PRECISION  A(*),  AK
        LOGICAL  DFLAG
C
C
C--PHASE 1 -- FIND ROW IN WHICH TO STORE EACH NONZERO
C----INITIALIZE COUNT OF NONZEROES TO BE STORED IN EACH ROW
        DO 1 I=1,N
  1       Q(I) = 0
C
C----FOR EACH NONZERO ELEMENT A(J)
        DO 3 I=1,N
          JMIN = IA(I)
          JMAX = IA(I+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 3
          DO 2 J=JMIN,JMAX
C
C--------FIND ROW (=R(J)) AND COLUMN (=JA(J)) IN WHICH TO STORE A(J) ...
            K = JA(J)
            IF (IP(K).LT.IP(I))  JA(J) = I
            IF (IP(K).GE.IP(I))  K = I
            R(J) = K
C
C--------... AND INCREMENT COUNT OF NONZEROES (=Q(R(J)) IN THAT ROW
  2         Q(K) = Q(K) + 1
  3       CONTINUE
C
C
C--PHASE 2 -- FIND NEW IA AND PERMUTATION TO APPLY TO (JA,A)
C----DETERMINE POINTERS TO DELIMIT ROWS IN PERMUTED (JA,A)
        DO 4 I=1,N
          IA(I+1) = IA(I) + Q(I)
  4       Q(I) = IA(I+1)
C
C----DETERMINE WHERE EACH (JA(J),A(J)) IS STORED IN PERMUTED (JA,A)
C----FOR EACH NONZERO ELEMENT (IN REVERSE ORDER)
        ILAST = 0
        JMIN = IA(1)
        JMAX = IA(N+1) - 1
        J = JMAX
        DO 6 JDUMMY=JMIN,JMAX
          I = R(J)
          IF (.NOT.DFLAG .OR. JA(J).NE.I .OR. I.EQ.ILAST)  GO TO 5
C
C------IF DFLAG, THEN PUT DIAGONAL NONZERO AT BEGINNING OF ROW
            R(J) = IA(I)
            ILAST = I
            GO TO 6
C
C------PUT (OFF-DIAGONAL) NONZERO IN LAST UNUSED LOCATION IN ROW
  5         Q(I) = Q(I) - 1
            R(J) = Q(I)
C
  6       J = J-1
C
C
C--PHASE 3 -- PERMUTE (JA,A) TO UPPER TRIANGULAR FORM (WRT NEW ORDERING)
        DO 8 J=JMIN,JMAX
  7       IF (R(J).EQ.J)  GO TO 8
            K = R(J)
            R(J) = R(K)
            R(K) = K
            JAK = JA(K)
            JA(K) = JA(J)
            JA(J) = JAK
            AK = A(K)
            A(K) = A(J)
            A(J) = AK
            GO TO 7
  8       CONTINUE
C
        RETURN
        END
c
       SUBROUTINE SSF
     *     (N, P,IP, IA,JA, IJU,JU,IU,JUMAX, Q, MARK, JL, FLAG)
C
C  ADDITIONAL PARAMETERS
C
C    Q     - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C    MARK  - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C    JL    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C
C  DEFINITIONS OF INTERNAL PARAMETERS (DURING K-TH STAGE OF ELIMINATION)
C
C    Q CONTAINS AN ORDERED LINKED LIST REPRESENTATION OF THE NONZERO
C      STRUCTURE OF THE K-TH ROW OF U --
C        Q(K) IS THE FIRST COLUMN WITH A NONZERO ENTRY
C        Q(I) IS THE NEXT COLUMN WITH A NONZERO ENTRY AFTER COLUMN I
C      IN EITHER CASE, Q(I) = N+1 INDICATES THE END OF THE LIST
C
C    JL CONTAINS LISTS OF ROWS TO BE MERGED INTO UNELIMINATED ROWS --
C        I GE K => JL(I) IS THE FIRST ROW TO BE MERGED INTO ROW I
C        I LT K => JL(I) IS THE ROW FOLLOWING ROW I IN SOME LIST OF ROWS
C      IN EITHER CASE, JL(I) = 0 INDICATES THE END OF A LIST
C
C    MARK(I) IS THE LAST ROW STORED IN JU FOR WHICH U(MARK(I),I) NE 0
C
C    JUMIN AND JUPTR ARE THE INDICES IN JU OF THE FIRST AND LAST
C      ELEMENTS IN THE LAST ROW SAVED IN JU
C
C    LUK IS THE NUMBER OF NONZERO ENTRIES IN THE K-TH ROW
C
C-----------------------------------------------------------------------
C
        INTEGER  P(*), IP(*),  IA(*), JA(*),  IJU(*), JU(*), IU(*),
     *      Q(*),  MARK(*),  JL(*),  FLAG,  TAG, VJ, QM
        LOGICAL  CLIQUE
C
C----INITIALIZATION
        JUMIN = 1
        JUPTR = 0
        IU(1) = 1
        DO 1 K=1,N
          MARK(K) = 0
   1      JL(K) = 0
C
C----FOR EACH ROW K
        DO 18 K=1,N
          LUK = 0
          Q(K) = N+1
C
          TAG = MARK(K)
          CLIQUE = .FALSE.
          IF (JL(K).NE.0)  CLIQUE = JL(JL(K)).EQ.0
C
C------   INITIALIZE NONZERO STRUCTURE OF K-TH ROW TO ROW P(K) OF M
          JMIN = IA(P(K))
          JMAX = IA(P(K)+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 4
          DO 3 J=JMIN,JMAX
            VJ = IP(JA(J))
            IF (VJ.LE.K)  GO TO 3
C
            QM = K
   2        M = QM
            QM = Q(M)
            IF (QM.LT.VJ)  GO TO 2
cww         IF (QM.EQ.VJ)  GO TO 102
            IF (QM.EQ.VJ)  then
C             ERROR -- DUPLICATE ENTRY IN A
              FLAG = 2*N + P(K)
              write(*,*) 'DUPLICATE ENTRY IN A', FLAG,P(K)
              RETURN
            end if

            LUK = LUK+1
            Q(M) = VJ
            Q(VJ) = QM
            IF (MARK(VJ).NE.TAG)  CLIQUE = .FALSE.
C
   3      CONTINUE
C
C------IF EXACTLY ONE ROW IS TO BE MERGED INTO THE K-TH ROW AND THERE IS
C------A NONZERO ENTRY IN EVERY COLUMN IN THAT ROW IN WHICH THERE IS A
C------NONZERO ENTRY IN ROW P(K) OF M, THEN DO NOT COMPUTE FILL-IN, JUST
C------USE THE COLUMN INDICES FOR THE ROW WHICH WAS TO HAVE BEEN MERGED
   4      IF (.NOT.CLIQUE)  GO TO 5
          IJU(K) = IJU(JL(K)) + 1
          LUK = IU(JL(K)+1) - (IU(JL(K))+1)
          GO TO 17
C
C------MODIFY NONZERO STRUCTURE OF K-TH ROW BY COMPUTING FILL-IN
C------FOR EACH ROW I TO BE MERGED IN
   5      LMAX = 0
          IJU(K) = JUPTR
C
          I = K
   6      I = JL(I)
          IF (I.EQ.0)  GO TO 10
C
C-------- MERGE ROW I INTO K-TH ROW
          LUI = IU(I+1) - (IU(I)+1)
          JMIN = IJU(I) +  1
          JMAX = IJU(I) + LUI
          QM = K
C
          DO 8 J=JMIN,JMAX
            VJ = JU(J)
   7        M = QM
            QM = Q(M)
            IF (QM.LT.VJ)  GO TO 7
            IF (QM.EQ.VJ)  GO TO 8
            LUK = LUK+1
            Q(M) = VJ
            Q(VJ) = QM
            QM = VJ
   8      CONTINUE
C
C-------- REMEMBER LENGTH AND POSITION IN JU OF LONGEST ROW MERGED
          IF (LUI.LE.LMAX)  GO TO 9
          LMAX = LUI
          IJU(K) = JMIN
C
   9      GO TO 6
C
C------IF THE K-TH ROW IS THE SAME LENGTH AS THE LONGEST ROW MERGED,
C------THEN USE THE COLUMN INDICES FOR THAT ROW
  10      IF (LUK.EQ.LMAX)  GO TO 17
C
C------IF THE TAIL OF THE LAST ROW SAVED IN JU IS THE SAME AS THE HEAD
C------OF THE K-TH ROW, THEN OVERLAP THE TWO SETS OF COLUMN INDICES --
C--------SEARCH LAST ROW SAVED FOR FIRST NONZERO ENTRY IN K-TH ROW ...
          I = Q(K)
          IF (JUMIN.GT.JUPTR)  GO TO 12
          DO 11 JMIN=JUMIN,JUPTR
            IF (JU(JMIN)-I)  11, 13, 12
  11      CONTINUE
  12      GO TO 15
C
C--------... AND THEN TEST WHETHER TAIL MATCHES HEAD OF K-TH ROW
  13      IJU(K) = JMIN
          DO 14 J=JMIN,JUPTR
            IF (JU(J).NE.I)  GO TO 15
            I = Q(I)
            IF (I.GT.N)  GO TO 17
  14      CONTINUE
          JUPTR = JMIN - 1
C
C------SAVE NONZERO STRUCTURE OF K-TH ROW IN JU
  15      I = K
          JUMIN = JUPTR +  1
          JUPTR = JUPTR + LUK
cww       IF (JUPTR.GT.JUMAX)  GO TO 106
          IF (JUPTR.GT.JUMAX)  then
C           ERROR -- INSUFFICIENT STORAGE FOR JU
            FLAG = 6*N + K
            write(*,*) 'INSUFFICIENT STORAGE IN SSF',FLAG,K
            return
          end if
          DO 16 J=JUMIN,JUPTR
            I = Q(I)
            JU(J) = I
  16        MARK(I) = K
          IJU(K) = JUMIN
C
C------ADD K TO ROW LIST FOR FIRST NONZERO ELEMENT IN K-TH ROW
  17      IF (LUK.LE.1)  GO TO 18
            I = JU(IJU(K))
            JL(K) = JL(I)
            JL(I) = K
C
  18      IU(K+1) = IU(K) + LUK
C
        FLAG = 0
        RETURN
C
cwwC ** ERROR -- DUPLICATE ENTRY IN A
cww 102    FLAG = 2*N + P(K)
cw        RETURN
cwwC ** ERROR -- INSUFFICIENT STORAGE FOR JU
cww 106    FLAG = 6*N + K
cww        RETURN
        END
C


