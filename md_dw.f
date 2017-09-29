        PROGRAM MDYN

C    *******************************************************************
C    *******************************************************************
C    **
C    **    THIS PROGRAM RUNS UPTO EQUILIBRATION AND STORES POSITIONS  **
C    **    AND VELOCITY OF THE PARTICLES IN A FILE output.out         **
C    **                                                               **
C    *******************************************************************
C    *******************************************************************


C    *******************************************************************
C    ** FORTRAN PROGRAM TO CONDUCT MOLECULAR DYNAMICS OF ATOMS.       **
C    **                                                               **
C    ** MICROCANONICAL ENSEMBLE                                       **
C    **                                                               **
C    ** APPLICATION TO ARGON. THE LENNARD-JONES (LJ) POTENTIAL IS     **
C    ** TRUNCATED AT RCUTF AND NOT SMOOTHLY CONTINUED TO ZERO.        **
C    ** INITIALLY THE N PARTICLES ARE PLACED IN AN FCC LATTICE.       **
C    ** THE VELOCITIES ARE DRAWN FROM A BOLTZMANN DISTRIBUTION WITH   **
C    ** TEMPARATURE TREF.                                             **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** D. W. HEERMANN, "COMPUTER SIMULATION METHODS"                 **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE POSITION( )                                        **
C    **    READS IN CONFIGURATION                                     **
C    **                                                               **
C    ** INPUT PARAMETERS ARE AS FOLLOWS                               **
C    **                                                               **
C    ** N              NUMBER OF PARTICLES (MUST BE A MULTIPLE OF 4)  **
C    ** BOXL           SIDE LENGTH OF THE CUBICAL BOX IN SIGMA UNIT   **
C    ** TREF           REDUCED TEMPERATURE                            **
C    ** RCUT           CUTOFF OF THE POTENTIAL IN SIGMA UNIT          **
C    ** DT             BASIC TIME STEP                                **
C    ** ISCALE         VELOCITY SCALING AT EVERY ISCALE'TH TIME STEP  **
C    ** ISTOP          STOP OF THE SCALING OF THE VELOCITIES AT ISTOP **
C    ** NSTEP          TOTAL NO. OF TIME STEPS                        **
C    ** ISEED          SEED FOR THE RANDOM NUMBER GENERATOR           **
C    **                                                               **
C    ** DIETER W. HEERMANN                                            **
C    *******************************************************************
C        PROGRAM FROGGY

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / VX, VY, VZ
        COMMON / BLOCK3 / FX, FY, FZ
	COMMON/BLOCK4/EPSLN,SIGM
        COMMON / BLOCKF1 /  DELR , HIST
        COMMON / BLOCKF2 /  MAXBIN, STEP, ISTOP, NSKIP, IAVSTART
C        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ



C    *******************************************************************
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 NUMBER OF ATOMS                     **
C    ** REAL    DT                TIMESTEP                            **
C    ** REAL    RX(N),RY(N),RZ(N) ATOMIC POSITIONS                    **
C    ** REAL    VX(N),VY(N),VZ(N) ATOMIC VELOCITIES                   **
C    ** REAL    ACV,ACK ETC.      AVERAGE VALUE ACCUMULATORS          **
C    ** REAL    AVV,AVK ETC.      AVERAGE VALUES                      **
C    ** REAL    ACVSQ, ACKSQ ETC. AVERAGE SQUARED VALUE ACCUMULATORS  **
C    ** REAL    FLV,FLK ETC.      FLUCTUATION AVERAGES                **
C    **                                                               **
C    **                                                               **
C    ** UNITS:                                                        **
C    **                                                               **
C    ** THE PROGRAM USES LENNARD-JONES REDUCED UNITS FOR USER INPUT   **
C    ** AND OUTPUT BUT CONDUCTS SIMULATION IN A BOX OF UNIT LENGTH.   **
C    ** SUMMARY FOR BOX LENGTH L, ATOMIC MASS M, AND LENNARD-JONES    **
C    ** POTENTIAL PARAMETERS SIGMA AND EPSILON:                       **
C    **                OUR PROGRAM           LENNARD-JONES SYSTEM     **
C    ** LENGTH         L                     SIGMA                    **
C    ** MASS           M                     M                        **
C    ** ENERGY         EPSILON               EPSILON                  **
C    ** TIME           SQRT(M*L**2/EPSILON)  SQRT(M*SIGMA**2/EPSILON) **
C    ** VELOCITY       SQRT(EPSILON/M)       SQRT(EPSILON/M)          **
C    ** PRESSURE       EPSILON/L**3          EPSILON/SIGMA**3         **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c       parameter (n =864)
	parameter (N=1372)
c	parameter (N=499)
c        PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

        REAL*8        RX(N), RY(N), RZ(N)
        REAL*8        VX(N), VY(N), VZ(N)
        REAL*8        FX(N), FY(N), FZ(N)
        REAL*8        HIST(10000), GR(10000), R(10000)
	REAL*8        EPSLN(N), SIGM(N)


        INTEGER     STEP, NSTEP, IPRINT
        INTEGER     NSKIP, NCOUNTGR
        REAL*8        DENS, NORM, RCUT, DT, SIGMA
        REAL*8        V, K, E, W
	 REAL*8        EPS11, SIG11, EPS22, SIG22
C       REAL*8        OLDK, NEWK, OLDV, NEWV, NEWW
C       REAL*8        OLDVC, NEWVC, VC, KC, EC
        REAL*8        VN, KN, EN, ECN, PRES, TEMP
        REAL*8        ACV, ACK, ACE, ACEC, ACP, ACT
        REAL*8        AVV, AVK, AVE, AVEC, AVP, AVT
        REAL*8        ACVSQ, ACKSQ, ACESQ, ACECSQ, ACPSQ, ACTSQ
        REAL*8        FLV, FLK, FLE, FLEC, FLP, FLT
        REAL*8        SR3, SR9, VLRC, WLRC, PI, SIGCUB
        CHARACTER   CNFILE*30, TITLE*80
        REAL*8        FREE
	real*8        MOLEFRAC1

        REAL*8        BOXL, AA, M, EPSLON
        INTEGER     ISCALE,ISTOP
        INTEGER     ISEED
        REAL*8        KN1,TEMPV,SCFACT2,SCFACT
        REAL*8        TREF,VC
        INTEGER     NCOUNT,ISCALESTART,NREADOPT
        REAL*8        BOXL2, DELR, CONST, RLOWER, RUPPER, NIDEAL
        INTEGER     BIN , MAXBIN
        INTEGER     nlabel, natom, IAVSTART
	  INTEGER     N1,N2
        PARAMETER ( FREE = 3.0 )
        PARAMETER ( PI = 3.1415927 )

        OPEN ( UNIT = 4, FILE = 'mdinp3.inp', STATUS = 'UNKNOWN')
        OPEN ( UNIT = 7, FILE = 'rxyz.out', STATUS = 'UNKNOWN')
        OPEN ( UNIT = 8, FILE = 'vxyz.out', STATUS = 'UNKNOWN')

C    *******************************************************************

C    ** READ IN INITIAL PARAMETERS **

C        WRITE(*,'(1H1,'' **** PROGRAM FROGGY ****                  '')')
        WRITE(*,'(//, '' MOLECULAR DYNAMICS OF LENNARD-JONES ATOMS '')')
        WRITE(*,'(    '' VELOCITY VERSION OF VERLET ALGORITHM      '')')

        WRITE(*,'('' ENTER RUN TITLE                               '')')
cc        READ (*,'(A)') TITLE
        READ (4,'(A)') TITLE
        WRITE(*,'('' ENTER NUMBER OF STEPS                         '')')
cc        READ (*,*) NSTEP
        READ (4,*) NSTEP
        WRITE(*,'('' ENTER NUMBER OF STEPS BETWEEN OUTPUT LINES    '')')
cc        READ (*,*) IPRINT
        READ (4,*) IPRINT
c        WRITE(*,'('' ENTER CONFIGURATION FILENAME                 '')')
c        READ (*,'(A)') CNFILE
        WRITE(*,'('' ENTER THE FOLLOWING IN LENNARD-JONES UNITS    '')')
        WRITE(*,'('' ENTER DENSITY & TEMPERATURE                   '')')
cc        READ (*,*) DENS,TREF
        READ (4,*) DENS,TREF
        WRITE(*,'('' ENTER POTENTIAL CUTOFF DISTANCE               '')')
cc        READ (*,*) RCUT
        READ (4,*) RCUT
        WRITE(*,'('' ENTER TIME STEP AND MASS                      '')')
cc        READ (*,*) DT,M
        READ (4,*) DT,M
	WRITE(*,'('' MOLEFRACTION                                  '')')
         READ(4,*)MOLEFRAC1
        WRITE(*,'(''EPSL11 AND SIG11                               '')')
         READ(4,*)EPS11,SIG11
        WRITE(*,'(''EPS22 AND SIG22                                '')')
        READ(4,*)EPS22,SIG22
        WRITE(*,'('' ENTER ISCALE AND ISTOP                        '')')
cc        READ (*,*) ISCALE , ISTOP
        READ (4,*) ISCALE , ISTOP
        READ (4,*) ISCALESTART , NREADOPT
C        WRITE(*,'('' ENTER BOX LENGTH                              '')')
C        READ (*,*) BOXL
        READ (4,*) DELR
        READ (4,*) NSKIP
        READ (4,*) IAVSTART

        WRITE(*,'(//1X,A)') TITLE
        WRITE(*,'('' NUMBER OF ATOMS  = '',I10  )') N
        WRITE(*,'('' NUMBER OF STEPS  = '',I10  )') NSTEP
        WRITE(*,'('' OUTPUT FREQUENCY = '',I10  )') IPRINT
        WRITE(*,'('' POTENTIAL CUTOFF = '',F10.4)') RCUT
        WRITE(*,'('' DENSITY          = '',F10.4)') DENS
        WRITE(*,'('' TEMPERATURE      = '',F10.4)') TREF
        WRITE(*,'('' TIME STEP AND MASS        = '',2F10.6)') DT,M
	 WRITE(*,'('' TIME STEP AND MASS        = '',2F10.6)')MOLEFRAC1
        WRITE(*,'('' TIME STEP AND MASS        = '',2F10.6)')EPS11,SIG11
       WRITE(*,'('' TIME STEP AND MASS        = '',2F10.6)')EPS22, SIG22
        WRITE(*,'('' VELOCITY SCALING AFTER STEP NO = '',I10)')ISCALE
        WRITE(*,'('' STOP VELOCITY SCALE AFTER STEP NO ='',I10)')ISTOP

        WRITE(10,'(//1X,A)') TITLE
        WRITE(10,'('' NUMBER OF ATOMS  = '',I10  )') N
        WRITE(10,'('' NUMBER OF STEPS  = '',I10  )') NSTEP
        WRITE(10,'('' OUTPUT FREQUENCY = '',I10  )') IPRINT
        WRITE(10,'('' POTENTIAL CUTOFF = '',F10.4)') RCUT
        WRITE(10,'('' DENSITY          = '',F10.4)') DENS
        WRITE(10,'('' TEMPERATURE      = '',F10.4)') TREF
        WRITE(10,'('' TIME STEP AND MASS        = '',2F10.6)') DT,M
	 WRITE(10,'('' TIME STEP AND MASS        = '',2F10.6)')MOLEFRAC1
       WRITE(10,'('' TIME STEP AND MASS        = '',2F10.6)')EPS11,SIG11
       WRITE(10,'('' TIME STEP AND MASS        = '',2F10.6)')EPS22,SIG22
        WRITE(10,'('' VELOCITY SCALING STARTS AT = '',I10)')ISCALESTART
        WRITE(10,'('' VELOCITY SCALING AT EVERY = '',I10)')ISCALE
        WRITE(10,'('' STOP VELOCITY SCALE AFTER STEP NO ='',I10)')ISTOP
        WRITE(10,'(''NREADOPT=1 READS R AND V FROM FILE'',I10)')NREADOPT

C    ** READ CONFIGURATION INTO COMMON / BLOCK1 / VARIABLES **

C**********************************************************************
C**BELOW ARE INITIAL CONFIGURATIONS FOR POSITIONS AND INITIAL        **
C**VALUES OF MOMENTA DRAWN FROM SIMPLE BOLTZMANN DISTRIBUTION        **
C**********************************************************************

               nlabel = 0
c            iseed=4711
c          call ranset ( iseed )

C    ** CONVERT INPUT DATA TO PROGRAM UNITS **

C      SIGMA = ( DENS / FLOAT ( N ) ) ** ( 1.0 / 3.0 )
C       RCUT  = RCUT * SIGMA
C       DT    = DT * SIGMA
C       DENS  = DENS / ( SIGMA ** 3 )

       BOXL  =  ( FLOAT (N) / DENS ) ** (1.0 / 3.0 )

        WRITE(*,'('' BOX LENGTH          = '', F10.6)') BOXL
        WRITE(10,'('' BOX LENGTH          = '', F10.6)') BOXL

        BOXL2  = BOXL/2.0
        MAXBIN = INT( BOXL2/DELR ) + 1
           write(*,*)'maxbin=',maxbin

               EPSLON = 1.0
               SIGMA  = 1.0

C        CALL READCN ( CNFILE )

C     ***************************************************
C     **                                               **
C     **   GENERATION OF INITIAL CONFIGURATION BELOW   **
C     **                                               **
C     ***************************************************

        IF ( NREADOPT .EQ. 0 ) THEN


c        AA=BOXL/4.0
c        CALL CONFIGIN(AA)

        CALL FCC ( BOXL )
        CALL COMVEL ( TREF )
C        CALL MAXWELL ( TREF )

        ELSE

        READ(7,*) NATOM
        READ(8,*) NATOM
           IF(NATOM.NE.N)THEN
             WRITE(*,*)'ERROR IN READING DATA', 'N=',N,'NATOM=',NATOM
             WRITE(10,*)'ERROR IN READING DATA', 'N=',N,'NATOM=',NATOM
           ENDIF
        DO I=1,N
           READ(7,*)RX(I),RY(I),RZ(I)
           READ(8,*)VX(I),VY(I),VZ(I)
        ENDDO

        ENDIF

        DO I=1,N
           WRITE(2,*)RX(I),RY(I),RZ(I)
           WRITE(3,*)VX(I),VY(I),VZ(I)
        ENDDO
C********************** MOLEFRACTION**********************************
          N1= int(MOLEFRAC1*N)
     
           N2 =N-N1
	write(12,*)N1, N2
      DO I=1,N
        IF(I.LE.N1)THEN
      		EPSLN(I)=EPS11
          	SIGM(I)=SIG11
       ELSE
      		EPSLN(I)=EPS22
          	SIGM(I)=SIG22
       ENDIF
       ENDDO

C***********************  END OF INITIALIZATION  *********************

C    ** CALCULATE LONG-RANGE CORRECTIONS **
C    ** NOTE: SPECIFIC TO LENNARD-JONES  **

        SR3    = ( SIGMA / RCUT ) ** 3
        SR9    = SR3 ** 3
        SIGCUB = SIGMA ** 3
        VLRC = ( 8.0 /9.0 ) * PI * DENS * SIGCUB * FLOAT ( N )
     :           * ( SR9 - 3.0 * SR3 )
       WLRC = ( 16.0 / 9.0 ) * PI * DENS * SIGCUB * FLOAT ( N )
     :           * ( 2.0 * SR9 - 3.0 * SR3 )

       WLRC = wlrc/FLOAT (n)
c       write(10,*)'wlrc=',wlrc

C    ** ZERO ACCUMULATORS **

        ACV  = 0.0
        ACK  = 0.0
        ACE  = 0.0
        ACEC = 0.0
        ACP  = 0.0
        ACT  = 0.0

        ACVSQ  = 0.0
        ACKSQ  = 0.0
        ACESQ  = 0.0
        ACECSQ = 0.0
        ACPSQ  = 0.0
        ACTSQ  = 0.0

        FLV  = 0.0
        FLK  = 0.0
        FLE  = 0.0
        FLEC = 0.0
        FLP  = 0.0
        FLT  = 0.0

                   NCOUNT = 0
                   NCOUNTGR = 0

C    ** INCLUDE LONG-RANGE CORRECTIONS **

        IF ( IPRINT .LE. 0 ) IPRINT = NSTEP + 1

        WRITE(*,'(//1X,''**** START OF DYNAMICS ****'')')
        WRITE(*,10001)
        WRITE(9,10001)

C    *******************************************************************
C    ** MAIN LOOP BEGINS                                              **
C    *******************************************************************

         WRITE(9,*)'STEP,  EN,  ECN,  KN,  VN,  PRES,  TEMP'

C         DO I=1,N
C          FX(I)=0.0
C          FY(I)=0.0
C          FZ(I)=0.0
C         ENDDO

        DO I=1,MAXBIN
          HIST(I)=0.0
        ENDDO

        DO 1000 STEP = 1, NSTEP

         IF ( MOD( STEP, 200 ) .EQ. 0 ) WRITE(*,*)'STEP NO.=', STEP

C       ** IMPLEMENT ALGORITHM **

           CALL MOVEA ( DT,M )

           CALL CUBICPB ( BOXL )

           CALL FORCE ( EPSLON, SIGMA, RCUT, BOXL, V, VC, W )

           CALL MOVEB ( DT, M, K )

C       ** INCLUDE LONG-RANGE CORRECTIONS **

           V = V + VLRC
c           W = W + WLRC

           kc=k

           E    = K + V
           EC   = KC + VC
           VN   = V  / FLOAT ( N )
           KN   = K  / FLOAT ( N )
           EN   = E  / FLOAT ( N )
           ECN  = EC / FLOAT ( N )
           TEMP = 2.0 * KN / FREE
c            write(102,*)w,wlrc
           PRES = DENS * TEMP + W / (boxl**3) 

C       ** CONVERT TO LENNARD-JONES UNITS **

           PRES = PRES * SIGMA ** 3



C       ************************************************
C       **                                            **
C       **         VELOCITY SCALING DONE BELOW        **
C       **                                            **
C       ************************************************

c             GO TO 101

         IF ( STEP.GE.ISCALESTART.AND.STEP.LT.ISTOP ) THEN

         IF ( MOD( STEP, ISCALE ) .EQ. 0 ) THEN

             WRITE(*,*)'VELOCITY SCALING / ADJUSTMENT'

             KN1     = K / FLOAT (N-1)
             TEMPV   = 2.0 * KN1 / ( FREE*EPSLON )
             SCFACT2 = TREF / TEMPV
             SCFACT  = SQRT (SCFACT2)
             WRITE(*,*)' SCALE FACTOR FOR VELOCITY AFTER',
     &         STEP, 'TH STEP=', SCFACT

             DO I = 1 , N

              VX(I) = VX(I) * SCFACT
              VY(I) = VY(I) * SCFACT
              VZ(I) = VZ(I) * SCFACT

             ENDDO
             K = K * SCFACT * SCFACT

c               IF( STEP .GT. 800 ) THEN
c               IF( ABS(SCFACT-1.0).LT.0.005) ISTOP = STEP
c               ENDIF

         ENDIF

         ENDIF

101      CONTINUE
C---------------------VELOCITY SCALING OVER------------------------

C       ***********************************************************
C       **     FOR COMPUTING AVERAGES OF DIFFERENT QUANTITIES    **
C       **     INCREMENT ACCUMULATORS                            **
C       ***********************************************************

           IF ( STEP .GT. (ISTOP+IAVSTART)) THEN
                   NCOUNT = NCOUNT + 1



           ACE  = ACE  + EN
           ACEC = ACEC + ECN
           ACK  = ACK  + KN
           ACV  = ACV  + VN
           ACP  = ACP  + PRES

           ACESQ  = ACESQ  + EN  ** 2
           ACECSQ = ACECSQ + ECN ** 2
           ACKSQ  = ACKSQ  + KN  ** 2
           ACVSQ  = ACVSQ  + VN  ** 2
           ACPSQ  = ACPSQ  + PRES ** 2

         ENDIF

        IF(STEP.GT.(ISTOP+IAVSTART).AND.MOD(STEP,NSKIP).EQ.0) THEN
           NCOUNTGR = NCOUNTGR + 1

           IF ( MOD( NCOUNTGR, 1000 ) .EQ. 0 ) THEN
              WRITE(19,*)NCOUNTGR
             DO I=1,MAXBIN
              WRITE(19,'(1X,i8,6(2X,E32.9))')I,HIST(I)
             ENDDO
            ENDIF
        ENDIF

C       ** OPTIONALLY PRINT INFORMATION **

           IF ( MOD( STEP, IPRINT ) .EQ. 0 ) THEN
c           IF (step.gt.7000 ) THEN
                    nlabel = nlabel + 1

              WRITE(*,'(1X,I8,6(2X,F10.4))')
     :                 STEP, EN, ECN, KN, VN, PRES, TEMP
              WRITE(9,'(1X,I8,6(2X,E19.9))')
     :                 STEP, E, EC, K, V, PRES, TEMP
                 DO I=1,N
              WRITE(23,'(1X,2I8,6(2X,F10.4))')
     :                STEP,I,RX(I),RY(I),RZ(I)
c              WRITE(21,'(1X,2I8,6(2X,F10.4))')
c     :                STEP,nlabel,VX(I),VY(I),VZ(I)
                 ENDDO
           ENDIF
c              WRITE(9,'(1X,I8,6(2X,F10.4))')
c     :                 STEP, EN, ECN, KN, VN, PRES, TEMP
cc              WRITE(9,'(1X,I8,6(2X,E19.9))')
cc     :                 STEP, E, EC, K, V, PRES, TEMP

c           IF ( MOD( STEP, 1000 ) .EQ. 0 ) THEN
c             DO I=1,MAXBIN
c              WRITE(19,'(1X,i8,6(2X,E19.9))')I,HIST(I)
c             ENDDO
c            ENDIF
c                 DO I=1,N
c              WRITE(23,'(1X,2I8,6(2X,F10.4))')
c     :                STEP,I,RX(I),RY(I),RZ(I)
c              WRITE(21,'(1X,2I8,6(2X,F10.4))')
c     :                STEP,nlabel,VX(I),VY(I),VZ(I)
c                 ENDDO

1000    CONTINUE

C    *******************************************************************
C    ** MAIN LOOP ENDS                                                **
C    *******************************************************************

        WRITE(*,'(/1X,''**** END OF DYNAMICS **** ''//)')

C    ** OUTPUT RUN AVERAGES **

C        NORM   = FLOAT ( NSTEP )
        NORM   = FLOAT ( NCOUNT )

        AVE  = ACE  / NORM
        AVEC = ACEC / NORM
        AVK  = ACK  / NORM
        AVV  = ACV  / NORM
        AVP  = ACP  / NORM

        ACESQ  = ( ACESQ  / NORM ) - AVE  ** 2
        ACECSQ = ( ACECSQ / NORM ) - AVEC ** 2
        ACKSQ  = ( ACKSQ  / NORM ) - AVK  ** 2
        ACVSQ  = ( ACVSQ  / NORM ) - AVV  ** 2
        ACPSQ  = ( ACPSQ  / NORM ) - AVP  ** 2

        IF ( ACESQ  .GT. 0.0 ) FLE  = SQRT ( ACESQ  )
        IF ( ACECSQ .GT. 0.0 ) FLEC = SQRT ( ACECSQ )
        IF ( ACKSQ  .GT. 0.0 ) FLK  = SQRT ( ACKSQ  )
        IF ( ACVSQ  .GT. 0.0 ) FLV  = SQRT ( ACVSQ  )
        IF ( ACPSQ  .GT. 0.0 ) FLP  = SQRT ( ACPSQ  )

        AVT = AVK * 2.0 / FREE
        FLT = FLK * 2.0 / FREE

        WRITE(*,*)'NCOUNT=',NCOUNT
        WRITE(10,*)'NCOUNT=',NCOUNT
        WRITE(*,*)'NORM=',NORM
        WRITE(10,*)'NORM=',NORM
        WRITE(*,*)'NCOUNTGR=',NCOUNTGR
        WRITE(10,*)'NCOUNTGR=',NCOUNTGR

        WRITE(*,'('' AVERAGES'',6(2X,F10.5))')
     :            AVE, AVEC, AVK, AVV, AVP, AVT
        WRITE(*,'('' FLUCTS  '',6(2X,F10.5))')
     :            FLE, FLEC, FLK, FLV, FLP, FLT

        WRITE(10,*)' AVE, AVEC, AVK, AVV, AVP, AVT'

        WRITE(10,'('' AVERAGES'',6(2X,F20.5))')
     :            AVE*n, AVEC*n, AVK*n, AVV*n, AVP, AVT


        WRITE(10,*)' FLE, FLEC, FLK, FLV, FLP, FLT'

        WRITE(10,'('' FLUCTS  '',6(2X,F10.5))')
     :            FLE, FLEC, FLK, FLV, FLP, FLT

C------------------------------------------------------------------
C       ** CALCULATION OF g(r) BELOW **
           CONST=4.0*PI*DENS/3.0
           DO 10 BIN = 1, MAXBIN
           RLOWER  = FLOAT (BIN-1)*DELR
           RUPPER  = RLOWER + DELR
           NIDEAL  = CONST*( RUPPER**3 - RLOWER**3 )
C           GR(BIN) = HIST(BIN)/FLOAT(NSTEP)/FLOAT(N)/NIDEAL
           GR(BIN) = HIST(BIN)/FLOAT(NCOUNTGR)/FLOAT(N)/NIDEAL
           WRITE(22,*)RLOWER, GR(BIN)
10         CONTINUE
C------------------------------------------------------------------

C    ** WRITE OUT FINAL CONFIGURATION **

c        CALL WRITCN ( CNFILE )
        CALL WRITCN2 (  )
C             DO I = 1, N
C             WRITE (13, *) RX(I), RY(I), RZ(I)
C             WRITE (14, *) VX(I), VY(I), VZ(I)
C             ENDDO

        STOP

10001   FORMAT(//1X,'TIMESTEP  ..ENERGY..  CUTENERGY.'
     :              '  ..KINETIC.  ..POTENT..',
     :              '  .PRESSURE.  ..TEMPER..'/)
        END



        SUBROUTINE READCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / VX, VY, VZ
        COMMON / BLOCK3 / FX, FY, FZ
C        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** SUBROUTINE TO READ IN INITIAL CONFIGURATION FROM UNIT 10      **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c          parameter (n =864)
	parameter (N=1372)
c	parameter (N=499)
c        PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

        CHARACTER   CNFILE*(*)
        REAL*8        RX(N), RY(N), RZ(N)
        REAL*8        VX(N), VY(N), VZ(N)
        REAL*8        FX(N), FY(N), FZ(N)

        INTEGER     CNUNIT, NN
        PARAMETER ( CNUNIT = 10 )

C     ******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'OLD',
     :         FORM = 'UNFORMATTED' )

        READ ( CNUNIT ) NN
        IF ( NN .NE. N ) STOP ' INCORRECT NUMBER OF ATOMS '
        READ ( CNUNIT ) RX, RY, RZ
        READ ( CNUNIT ) VX, VY, VZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE WRITCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ
cc        COMMON / BLOCK2 / VX, VY, VZ
cc        COMMON / BLOCK3 / FX, FY, FZ
C        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** ROUTINE TO WRITE OUT FINAL CONFIGURATION TO UNIT 10           **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c       parameter (n =864)
	parameter (N=1372)
c	parameter (n=499)
c         PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

c        CHARACTER   CNFILE*(*)
        CHARACTER   CNFILE*30
        REAL*8        RX(N), RY(N), RZ(N)
        REAL*8        VX(N), VY(N), VZ(N)
        REAL*8        FX(N), FY(N), FZ(N)

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 11 )

C    ****************************************************************

        OPEN  ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'UNKNOWN',
     :          FORM = 'UNFORMATTED' )

        WRITE ( CNUNIT ) N

        WRITE ( CNUNIT ) RX, RY, RZ
c        WRITE ( CNUNIT ) VX, VY, VZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE KINET ( K )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / VX, VY, VZ
C        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** COMPUTES KINETIC ENERGY                                       **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c       parameter (n =864)
	parameter (N=1372)
c	parameter (N=499)
c        PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

        REAL*8        RX(N), RY(N), RZ(N)
        REAL*8        VX(N), VY(N), VZ(N)
        REAL*8        K

        INTEGER     I

C    *******************************************************************

        K  = 0.0

        DO 100 I = 1, N

           K = K + VX(I) ** 2 + VY(I) ** 2 + VZ(I) ** 2

100     CONTINUE

        K = K * 0.5

        RETURN
        END


C----------------------------------------------------------------------
        SUBROUTINE CONFIGIN(AA)

        COMMON / BLOCK1 / RX, RY, RZ
C        COMMON / BLOCK2 / VX, VY, VZ
C        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** GENERATES INITIAL POSITION                                    **
C    ** SET UP FCC LATTICE FOR THE ATOMS INSIDE THE BOX               **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c         parameter (n =864)
	parameter (N=1372)
c	parameter (n=499)
c        PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

        REAL*8        RX(N), RY(N), RZ(N)
        REAL*8        VX(N), VY(N), VZ(N)
C        REAL*8        K

        INTEGER     IJK,I,J,K,LG,L

C    *******************************************************
C    *  FOR 108 PARTICLES L=2 AND FOR 256 PARTICLES L=3    *
C    *******************************************************
        L=3
        IJK  = 0

        DO 10 LG = 0, 1
         DO 10 I = 0, L
          DO 10 J = 0, L
           DO 10 K = 0, L
            IJK=IJK+1
           RX(IJK)=I*AA+LG*AA*0.5D0
           RY(IJK)=J*AA+LG*AA*0.5D0
           RZ(IJK)=K*AA

10     CONTINUE

        DO 15 LG = 1, 2
         DO 15 I = 0, L
          DO 15 J = 0, L
           DO 15 K = 0, L
            IJK=IJK+1
           RX(IJK)=I*AA+(2-LG)*LG*AA*0.5D0
           RY(IJK)=J*AA+(LG-1)*LG*AA*0.5D0
           RZ(IJK)=K*AA+AA*0.5D0

15     CONTINUE

        RETURN
        END



C*************************************************************************
C** FICHE F.24.  INITIAL VELOCITY DISTRIBUTION                          **
C** This FORTRAN code is intended to illustrate points made in the text.**
C** To our knowledge it works correctly.                                **
C** However it is the responsibility of the user to test it, if it is   **
C** to be used in a research application.                               **
C*************************************************************************

C    *******************************************************************
C    ** CENTRE OF MASS AND ANGULAR VELOCITIES FOR LINEAR MOLECULES    **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   THE NUMBER OF MOLECULES           **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
C    ** REAL    VX(N),VY(N),VZ(N)   VELOCITIES                        **
C    ** REAL    EX(N),EY(N),EZ(N)   ORIENTATIONS                      **
C    ** REAL    OX(N),OY(N),OZ(N)   SPACE-FIXED ANGULAR VELOCITIES    **
C    ** REAL    TEMP                REDUCED TEMPERATURE               **
C    ** REAL    INERT               REDUCED MOMENT OF INERTIA         **
C    **                                                               **
C    ** SUPPLIED ROUTINES:                                            **
C    **                                                               **
C    ** SUBROUTINE COMVEL ( TEMP )                                    **
C    **    SETS THE CENTRE OF MASS VELOCITIES FOR A CONFIGURATION OF  **
C    **    LINEAR MOLECULES AT A GIVEN TEMPERATURE.                   **
C    ** SUBROUTINE ANGVEL ( TEMP, INERT )                             **
C    **    SETS THE ANGULAR VELOCITIES FOR A CONFIGURATION OF LINEAR  **
C    **    MOLECULES AT A GIVEN TEMPERATURE.                          **
C    ** REAL FUNCTION RANF ( DUMMY )                                  **
C    **    RETURNS A UNIFORM RANDOM VARIATE ON THE RANGE ZERO TO ONE  **
C    ** REAL FUNCTION GAUSS ( DUMMY )                                 **
C    **    RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             **
C    **    DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.             **
C    **                                                               **
C    ** UNITS:                                                        **
C    **                                                               **
C    ** WE ASSUME UNIT MOLECULAR MASS AND EMPLOY LENNARD-JONES UNITS  **
C    **       PROPERTY                      UNITS                     **
C    **       RX, RY, RZ           (EPSILON/M)**(1.0/2.0)             **
C    **       OX, OY, OZ           (EPSILON/M*SIGMA**2)**(1.0/2.0)    **
C    **       INERT                 M*SIGMA**2                        **
C    *******************************************************************

        SUBROUTINE COMVEL ( TEMP )

c        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / VX, VY, VZ
C        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** TRANSLATIONAL VELOCITIES FROM MAXWELL-BOLTZMANN DISTRIBUTION  **
C    **                                                               **
C    ** THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND (UNIT) MASS.**
C    ** THIS ROUTINE IS GENERAL, AND CAN BE USED FOR ATOMS, LINEAR    **
C    ** MOLECULES, AND NON-LINEAR MOLECULES.                          **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           **
C    **                                                               **
C    ** REAL FUNCTION GAUSS ( DUMMY )                                 **
C    **    RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             **
C    **    DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.             **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c         parameter (n =864)
	parameter (N=1372)
c	parameter (n=499)
c        PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

c        REAL*8        RX(N), RY(N), RZ(N)
        REAL*8        VX(N), VY(N), VZ(N)
        REAL*8        TEMP
        REAL*8        TEMPV,SC2,SC,EKIN

        REAL*8        RTEMP, SUMX, SUMY, SUMZ
        REAL*8        GAUSS, DUMMY
        INTEGER     I

C    *******************************************************************

        RTEMP = SQRT ( TEMP )

        DO 100 I = 1, N

           VX(I) = RTEMP * GAUSS ( DUMMY )
           VY(I) = RTEMP * GAUSS ( DUMMY )
           VZ(I) = RTEMP * GAUSS ( DUMMY )

100     CONTINUE

C    ** REMOVE NET MOMENTUM **

        SUMX = 0.0
        SUMY = 0.0
        SUMZ = 0.0

        DO 200 I = 1, N

           SUMX = SUMX + VX(I)
           SUMY = SUMY + VY(I)
           SUMZ = SUMZ + VZ(I)

200     CONTINUE

        SUMX = SUMX / FLOAT ( N )
        SUMY = SUMY / FLOAT ( N )
        SUMZ = SUMZ / FLOAT ( N )

           EKIN=0.0
        DO 300 I = 1, N

           VX(I) = VX(I) - SUMX
           VY(I) = VY(I) - SUMY
           VZ(I) = VZ(I) - SUMZ

           EKIN=EKIN+(VX(I)*VX(I)+VY(I)*VY(I)+VZ(I)*VZ(I))


300     CONTINUE
         write(*,*)'ekin in comvel before scaling=',0.5*ekin
         write(10,*)'ekin in comvel before scaling=',0.5*ekin

c           SC=1.0
           EKIN=0.5*EKIN
           WRITE(*,*)'VELOCITY SCALING'
           TEMPV=2.0*EKIN/(3.0*(N-1))
           SC2=TEMP/TEMPV
           SC=SQRT(SC2)
           WRITE(*,*)'SCALE FACTOR IN COMVEL IS', SC
           ekin2=0.0
        DO 400 I = 1, N

           VX(I) = VX(I)*SC
           VY(I) = VY(I)*SC
           VZ(I) = VZ(I)*SC
           EKIN2=EKIN2+(VX(I)*VX(I)+VY(I)*VY(I)+VZ(I)*VZ(I))


400     CONTINUE
         write(*,*)'ekin in comvel after scaling=',0.5*ekin2
         write(10,*)'ekin in comvel after scaling=',0.5*ekin2
        RETURN
        END


        REAL*8 FUNCTION GAUSS ( DUMMY )

C    *******************************************************************
C    ** RANDOM VARIATE FROM THE STANDARD NORMAL DISTRIBUTION.         **
C    **                                                               **
C    ** THE DISTRIBUTION IS GAUSSIAN WITH ZERO MEAN AND UNIT VARIANCE.**
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION        **
C    **    ADDISON-WESLEY), 1978                                      **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           **
C    **                                                               **
C    ** REAL*8 FUNCTION RANF ( DUMMY )                                  **
C    **    RETURNS A UNIFORM RANDOM VARIATE ON THE RANGE ZERO TO ONE  **
C    *******************************************************************

        REAL*8        A1, A3, A5, A7, A9
        PARAMETER ( A1 = 3.949846138, A3 = 0.252408784 )
        PARAMETER ( A5 = 0.076542912, A7 = 0.008355968 )
        PARAMETER ( A9 = 0.029899776                   )

        REAL*8        SUM, R, R2
        REAL*8        RANF, DUMMY
        INTEGER     I

C    *******************************************************************

        SUM = 0.0

        DO 10 I = 1, 12

           SUM = SUM + RANF ( DUMMY )

10      CONTINUE

        R  = ( SUM - 6.0 ) / 4.0
        R2 = R * R

        GAUSS = (((( A9 * R2 + A7 ) * R2 + A5 ) * R2 + A3 ) * R2 +A1 )
     :          * R

        RETURN
        END



        REAL*8 FUNCTION RANF ( DUMMY )

C    *******************************************************************
C    ** RETURNS A UNIFORM RANDOM VARIATE IN THE RANGE 0 TO 1.         **
C    **                                                               **
C    **                 ***************                               **
C    **                 **  WARNING  **                               **
C    **                 ***************                               **
C    **                                                               **
C    ** GOOD RANDOM NUMBER GENERATORS ARE MACHINE SPECIFIC.           **
C    ** PLEASE USE THE ONE RECOMMENDED FOR YOUR MACHINE.              **
C    *******************************************************************

        INTEGER     L, C, M
        PARAMETER ( L = 1029, C = 221591, M = 1048576 )

        INTEGER     SEED
        REAL*8        DUMMY
        SAVE        SEED
        DATA        SEED / 0 /

C    *******************************************************************

        SEED = MOD ( SEED * L + C, M )
        RANF = FLOAT ( SEED ) / M

        RETURN
        END


C***********************************************************************
C** FICHE F.4.  VELOCITY VERSION OF VERLET ALGORITHM                  **
C** This FORTRAN code is intended to illustrate points made in        **
C** the text. To our knowledge it works correctly.                    **
C** However it is the responsibility of the user to test it, if it is **
C** to be used in a research application.                             **
C***********************************************************************

C    *******************************************************************
C    ** TWO ROUTINES THAT TOGETHER IMPLEMENT VELOCITY VERLET METHOD.  **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** SWOPE ET AL., J. CHEM. PHYS. 76, 637, 1982.                   **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE MOVEA ( DT, M )                                    **
C    **    MOVES POSITIONS AND PARTIALLY UPDATES VELOCITIES.          **
C    ** SUBROUTINE MOVEB ( DT, M, K )                                 **
C    **    COMPLETES VELOCITY MOVE AND CALCULATES KINETIC ENERGY.     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF MOLECULES               **
C    ** REAL    DT                  TIMESTEP                          **
C    ** REAL    M                   ATOMIC MASS                       **
C    ** REAL    K                   KINETIC ENERGY                    **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
C    ** REAL    VX(N),VY(N),VZ(N)   VELOCITIES                        **
C    ** REAL    FX(N),FY(N),FZ(N)   FORCES                            **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** AT THE START OF A TIMESTEP, MOVEA IS CALLED TO ADVANCE THE    **
C    ** POSITIONS AND 'HALF-ADVANCE' THE VELOCITIES.  THEN THE FORCE  **
C    ** ROUTINE IS CALLED, AND THIS IS FOLLOWED BY MOVEB WHICH        **
C    ** COMPLETES THE ADVANCEMENT OF VELOCITIES.                      **
C    *******************************************************************



        SUBROUTINE MOVEA ( DT, M )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / VX, VY, VZ
        COMMON / BLOCK3 / FX, FY, FZ

C        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ

C    *******************************************************************
C    ** FIRST PART OF VELOCITY VERLET ALGORITHM                       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE FIRST PART OF THE ALGORITHM IS A TAYLOR SERIES WHICH      **
C    ** ADVANCES POSITIONS FROM T TO T + DT AND VELOCITIES FROM       **
C    ** T TO T + DT/2.  AFTER THIS, THE FORCE ROUTINE IS CALLED.      **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c        parameter (n =864)
	parameter (N=1372)
c	parameter (n=499)
c        PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

        REAL*8        DT, M
        REAL*8        RX(N), RY(N), RZ(N)
        REAL*8        VX(N), VY(N), VZ(N)
        REAL*8        FX(N), FY(N), FZ(N)

        INTEGER     I
        REAL*8        DT2, DTSQ2

C    *******************************************************************

        DT2   = DT / 2.0
        DTSQ2 = DT * DT2

        DO 100 I = 1, N

           RX(I) = RX(I) + DT * VX(I) + DTSQ2 * FX(I) / M
           RY(I) = RY(I) + DT * VY(I) + DTSQ2 * FY(I) / M
           RZ(I) = RZ(I) + DT * VZ(I) + DTSQ2 * FZ(I) / M
           VX(I) = VX(I) + DT2 * FX(I) / M
           VY(I) = VY(I) + DT2 * FY(I) / M
           VZ(I) = VZ(I) + DT2 * FZ(I) / M

100     CONTINUE

        RETURN
        END



        SUBROUTINE MOVEB ( DT, M, K )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / VX, VY, VZ
        COMMON / BLOCK3 / FX, FY, FZ
C        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ

C    *******************************************************************
C    ** SECOND PART OF VELOCITY VERLET ALGORITHM                      **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE SECOND PART OF THE ALGORITHM ADVANCES VELOCITIES FROM     **
C    ** T + DT/2 TO T + DT. THIS ASSUMES THAT FORCES HAVE BEEN        **
C    ** COMPUTED IN THE FORCE ROUTINE AND STORED IN FX, FY, FZ.       **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c        parameter (n =864)
	parameter (N=1372)
c	parameter (n=499)
c        PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

        REAL*8        DT, M, K
        REAL*8        RX(N), RY(N), RZ(N)
        REAL*8        VX(N), VY(N), VZ(N)
        REAL*8        FX(N), FY(N), FZ(N)

        INTEGER     I
        REAL*8        DT2

C    *******************************************************************

        DT2 = DT / 2.0

        K = 0.0

        DO 200 I = 1, N

           VX(I) = VX(I) + DT2 * FX(I) / M
           VY(I) = VY(I) + DT2 * FY(I) / M
           VZ(I) = VZ(I) + DT2 * FZ(I) / M

           K = K + VX(I) ** 2 + VY(I) ** 2 + VZ(I) ** 2

200     CONTINUE

        K = 0.5 * M * K

        RETURN
        END


C**********************************************************************
C**        FICHE F.17.  A SIMPLE LENNARD-JONES FORCE ROUTINE         **
C**        This FORTRAN code is intended to illustration             **
C**        made in the text.  To our knowledge it works correctly.   **
C**        However it is the responsibility of the user to test it,  **
C**        if it is to be used in a research application.            **
C**********************************************************************



        SUBROUTINE FORCE ( EPSLON, SIGMA, RCUT, BOX, V, VC, W )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / VX, VY, VZ
        COMMON / BLOCK3 / FX, FY, FZ
	COMMON/BLOCK4/ EPSLN,SIGM
        COMMON / BLOCKF1 /  DELR , HIST
        COMMON / BLOCKF2 /  MAXBIN, STEP, ISTOP, NSKIP, IAVSTART
cC        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ

C    *******************************************************************
C    ** FORCE CALCULATION FOR LENNARD-JONES ATOMS.                    **
C    **                                                               **
C    ** IN THIS WE AIM TO SHOW HOW THE FORCES, POTENTIAL ENERGY AND   **
C    ** VIRIAL FUNCTION ARE CALCULATED IN A FAIRLY EFFICIENT WAY.     **
C    ** UNDOUBTEDLY FURTHER IMPROVEMENT WOULD BE POSSIBLE ON SPECIFIC **
C    ** MACHINES.                                                     **
C    ** THE POTENTIAL IS V(R) = 4*EPSLON*((SIGMA/R)**12-(SIGMA/R)**6) **
C    ** WE INCLUDE SPHERICAL CUTOFF AND MINIMUM IMAGING IN CUBIC BOX. **
C    ** THE BOX LENGTH IS BOX.  THE CUTOFF IS RCUT.                   **
C    ** THE ROUTINE ACTUALLY RETURNS TWO DIFFERENT POTENTIAL ENERGIES.**
C    ** V IS CALCULATED USING THE LENNARD-JONES POTENTIAL TO BE USED  **
C    ** FOR CALCULATING THE THERMODYNAMIC INTERNAL ENERGY.            **
C    ** LONG-RANGE CORRECTIONS SHOULD BE APPLIED TO THIS OUTSIDE THE  **
C    ** ROUTINE, IN THE FORM                                          **
C    **         SR3 = ( SIGMA / RCUT ) ** 3                           **
C    **         SR9 = SR3 ** 3                                        **
C    **         DENS = REAL(N) * ( SIGMA / BOX ) ** 3                 **
C    **         VLRC = ( 8.0 /9.0 ) * PI * EPSLON * DENS * REAL ( N ) **
C    **      :           * ( SR9 - 3.0 * SR3 )                        **
C    **         WLRC = ( 16.0 / 9.0 ) * PI * EPSLON * DENS * REAL( N )**
C    **      :           * ( 2.0 * SR9 - 3.0 * SR3 )                  **
C    **         V = V + VLRC                                          **
C    **         W = W + WLRC                                          **
C    ** VC IS CALCULATED USING THE SHIFTED LENNARD-JONES POTENTIAL,   **
C    ** WITH NO DISCONTINUITY AT THE CUTOFF, TO BE USED IN ASSESSING  **
C    ** ENERGY CONSERVATION.                                          **
C    ** NO REDUCED UNITS ARE USED: FOR THIS POTENTIAL WE COULD SET    **
C    ** EPSLON = 1 AND EITHER SIGMA = 1 OR BOX = 1 TO IMPROVE SPEED.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 NUMBER OF MOLECULES                 **
C    ** REAL    RX(N),RY(N),RZ(N) MOLECULAR POSITIONS                 **
C    ** REAL    VX(N),VY(N),VZ(N) MOLECULAR VELOCITIES (NOT USED)     **
C    ** REAL    FX(N),FY(N),FZ(N) MOLECULAR FORCES                    **
C    ** REAL    SIGMA             PAIR POTENTIAL LENGTH PARAMETER     **
C    ** REAL    EPSLON            PAIR POTENTIAL ENERGY PARAMETER     **
C    ** REAL    RCUT              PAIR POTENTIAL CUTOFF               **
C    ** REAL    BOX               SIMULATION BOX LENGTH               **
C    ** REAL    V                 POTENTIAL ENERGY                    **
C    ** REAL    VC                SHIFTED POTENTIAL                   **
C    ** REAL    W                 VIRIAL FUNCTION                     **
C    ** REAL    VIJ               PAIR POTENTIAL BETWEEN I AND J      **
C    ** REAL    WIJ               NEGATIVE OF PAIR VIRIAL FUNCTION W  **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c          parameter (n =864)
	parameter (N=1372)
c	parameter (n=499)
c        PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

        REAL*8        SIGMA, EPSLON, RCUT, BOX, V, VC, W
        REAL*8        RX(N), RY(N), RZ(N)
        REAL*8        VX(N), VY(N), VZ(N)
        REAL*8        FX(N), FY(N), FZ(N)
	 REAL*8        EPSLN(N),SIGM(N)
        REAL*8        HIST(10000)

        INTEGER     I, J, NCUT
        REAL*8        BOXINV, RCUTSQ, SIGSQ, EPS4, EPS24
	REAL*8        EPSFACT4 , EPSFACT24
        REAL*8        EPSFACT , EPSFACT2, EPSFACT6, EPSFACT12
        REAL*8        SIGFACT , SIGFACT2, SIGFACT6, SIGFACT12

        REAL*8        RXI, RYI, RZI, FXI, FYI, FZI
        REAL*8        RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
        REAL*8        SR2, SR6, SR12, VIJ, WIJ, FIJ
        INTEGER     BIN, MAXBIN
        INTEGER     NSKIP, STEP, ISTOP, IAVSTART
        REAL*8        DELR,TERMEXP
        REAL*8         L RNOT,T,f
C    *******************************************************************

C    ** CALCULATE USEFUL QUANTITIES **

        BOXINV = 1.0 / BOX
        RCUTSQ = RCUT ** 2
        SIGSQ  = SIGMA ** 2
        EPS4   = EPSLON * 4.0
        EPS24  = EPSLON * 24.0
	EPS4=1.0
        EPS24=1.0
c	Z=0
C    ** ZERO FORCES, POTENTIAL, VIRIAL **

        DO 100 I = 1, N

           FX(I) = 0.0
           FY(I) = 0.0
           FZ(I) = 0.0

100     CONTINUE

        NCUT = 0
        V    = 0.0
        W    = 0.0
        L    = 2.0
        RNOT = 1.1
        T  =   11.0

C    ** OUTER LOOP BEGINS **

        DO 200 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)
           FXI = FX(I)
           FYI = FY(I)
           FZI = FZ(I)
	EPSLNI=EPSLN(I)
           SIGMAI=SIGM(I)

C       ** INNER LOOP BEGINS **

           DO 199 J = I + 1, N

              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J)
              RZIJ = RZI - RZ(J)
              RXIJ = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
              RYIJ = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
              RZIJ = RZIJ - ANINT ( RZIJ * BOXINV ) * BOX
              RIJSQ = RXIJ ** 2 + RYIJ ** 2 + RZIJ ** 2

c******************************************************************
	EPSLNIJ=SQRT(EPSLNI*EPSLN(J))
              SIGMAIJ=(SIGMAI+SIGM(J))/2.0
              EPSFACT=EPSLNIJ/EPSLN(J)
                 EPSFACT4=4.0*EPSFACT
                 EPSFACT24=24.0*EPSFACT
 
             SIGFACT=SIGMAIJ/SIGM(J)
              SIGFACT2=SIGFACT*SIGFACT
              SIGFACT6=SIGFACT2*SIGFACT2*SIGFACT2
              SIGFACT12=SIGFACT6*SIGFACT6
c********************************************************************
C      LINES BELOW ARE FOR g(r) CALCULATION

              RIJ = SQRT ( RIJSQ )
              BIN = INT ( RIJ/DELR ) + 1
        IF(STEP.GT.(ISTOP+IAVSTART) .AND. MOD(STEP,NSKIP).EQ.0) THEN
          IF ( BIN .LE. MAXBIN ) THEN
              HIST ( BIN ) = HIST (BIN) + 2.0
          ENDIF
        ENDIF
c              write(*,*)'bin=',bin
C      ----------------------------------------

              IF ( RIJSQ .LT. RCUTSQ ) THEN
c                if(rijsq.eq.0.0)rijsq=1.0

                 SR2   = SIGSQ / RIJSQ
                 SR6   = SR2 * SR2 * SR2
                 SR12  = SR6 ** 2
c                 VIJ   = SR12 - SR6
              TERMEXP = - T * (RIJ - RNOT)**2
	vexp=L*EPSLNIJ*DEXP(TERMEXP)
c	write(11,*)EPSFACT4, EPSFACT24, SIGFACT12 , SIGFACT6
c       VIJ =EPSFACT4*(SR12*SIGFACT12-SR6*SIGFACT6)+vexp
	 VIJ =EPSFACT4*(SR12*SIGFACT12-SR6*SIGFACT6)+vexp
c	VIJ =(SR12-SR6)
 		
c		write(11,*)RIJ,vij
c                V     = V + VIJ
	 V= V+VIJ
c	f=L *DEXP(TERMEXP)*2.0*T*(RIJSQ-RNOT)/3.0
        WIJ =EPSFACT24*(2*SR12*SIGFACT12-SR6*SIGFACT6)+vexp*RIJ*
     :				2.0*T*(RIJ-RNOT)
                 W     = W + WIJ
                 FIJ   = WIJ / RIJSQ
                 FXIJ  = FIJ * RXIJ
                 FYIJ  = FIJ * RYIJ
                 FZIJ  = FIJ * RZIJ
                 FXI   = FXI + FXIJ
                 FYI   = FYI + FYIJ
                 FZI   = FZI + FZIJ
                 FX(J) = FX(J) - FXIJ
                 FY(J) = FY(J) - FYIJ
                 FZ(J) = FZ(J) - FZIJ
                 NCUT  = NCUT + 1

              ENDIF

199        CONTINUE

C       ** INNER LOOP ENDS **

           FX(I) = FXI
           FY(I) = FYI
           FZ(I) = FZI

200     CONTINUE

C    ** OUTER LOOP ENDS **

C    ** CALCULATE SHIFTED POTENTIAL **

c        SR2 = SIGSQ / RCUTSQ
c       SR6 = SR2 * SR2 * SR2
c       SR12 = SR6 * SR6
c        VIJ = SR12 - SR6
c        VIJ = (SR12*SIGFACT12-SR6*SIGFACT6)+L*epslon * DEXP(TERMEXP) 
c       VC = V - FLOAT ( NCUT ) * VIJ

C    ** MULTIPLY RESULTS BY ENERGY FACTORS **

c        DO 300 I = 1, N
c
c           FX(I) = FX(I) * EPS24
c           FY(I) = FY(I) * EPS24
c           FZ(I) = FZ(I) * EPS24

c       300     CONTINUE

c        V  = V  * EPS4
c        VC = VC * EPS4
        W  = W   / 3.0

        RETURN
        END


C  *********************************************************************
C  ** FICHE F.23.  ROUTINE TO SET UP ALPHA FCC LATTICE OF LINEAR      **
C  ** MOLECULES. This FORTRAN code is intended to illustrate points   **
C  ** made in the text.  To our knowledge it works correctly.         **
C  ** However it is the responsibility of the user to test it, if it  **
C  ** is to be used in a research application.                        **
C  *********************************************************************

        SUBROUTINE FCC ( BOXL )

            COMMON/BLOCK1/RX,RY,RZ
C        COMMON / BLOCK1 / RX, RY, RZ, EX, EY, EZ

C    *******************************************************************
C    ** SETS UP THE ALPHA FCC LATTICE FOR N LINEAR MOLECULES.         **
C    **                                                               **
C    ** THE SIMULATION BOX IS A UNIT CUBE CENTRED AT THE ORIGIN.      **
C    ** N SHOULD BE AN INTEGER OF THE FORM ( 4 * ( NC ** 3 ) ),       **
C    ** WHERE NC IS THE NUMBER OF FCC UNIT CELLS IN EACH DIRECTION.   **
C    ** SEE FIGURE 5.10 FOR A DIAGRAM OF THE LATTICE AND A            **
C    ** DEFINITION OF THE FOUR ORIENTATIONAL SUBLATTICES.             **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                    NUMBER OF MOLECULES              **
C    ** REAL    RX(N),RY(N),RZ(N)    MOLECULAR POSITIONS              **
C    ** REAL    EX(N),EY(N),EZ(N)    UNIT VECTORS GIVING ORIENTATIONS **
C    ** REAL    RROOT3               1.0 / SQRT ( 3.0 )               **
C    *******************************************************************

        INTEGER     N, NC
        REAL*8        RROOT3

c        PARAMETER ( NC = 3, N = 4 * NC ** 3 )
c         parameter (nc=5, n=4*nc**3)
         parameter (nc=7, n=4*nc**3)
c        PARAMETER ( NC = 4, N = 4 * NC ** 3 )
c        PARAMETER ( NC = 12, N = 4 * NC ** 3 )
c        PARAMETER ( NC = 10, N = 4 * NC ** 3 )
        PARAMETER ( RROOT3 = 0.5773503 )

        REAL*8        RX(N), RY(N), RZ(N)
C        REAL*8        EX(N), EY(N), EZ(N)
        REAL*8        CELL, CELL2
        REAL*8        BOXL
        INTEGER     I, IX, IY, IZ, IREF, M

C    *******************************************************************

C    ** CALCULATE THE SIDE OF THE UNIT CELL **

c        CELL  = 1.0 / FLOAT ( NC )
        CELL  = BOXL / FLOAT ( NC )
        CELL2 = 0.5 * CELL
C         WRITE(7,*)'BOXL,CELL=',BOXL,CELL

C    ** BUILD THE UNIT CELL **

C    ** SUBLATTICE A **

        RX(1) =  0.0
        RY(1) =  0.0
        RZ(1) =  0.0
c        EX(1) =  RROOT3
c        EY(1) =  RROOT3
c        EZ(1) =  RROOT3

C    ** SUBLATTICE B **

        RX(2) =  CELL2
        RY(2) =  CELL2
        RZ(2) =  0.0
c        EX(2) =  RROOT3
c        EY(2) = -RROOT3
c        EZ(2) = -RROOT3

C    ** SUBLATTICE C **

        RX(3) =  0.0
        RY(3) =  CELL2
        RZ(3) =  CELL2
c        EX(3) = -RROOT3
c        EY(3) =  RROOT3
c        EZ(3) = -RROOT3

C    ** SUBLATTICE D **

        RX(4) =  CELL2
        RY(4) =  0.0
        RZ(4) =  CELL2
c        EX(4) = -RROOT3
c        EY(4) = -RROOT3
c        EZ(4) =  RROOT3

C    ** CONSTRUCT THE LATTICE FROM THE UNIT CELL **

        M = 0

        DO 99 IZ = 1, NC

           DO 98 IY = 1, NC

              DO 97 IX = 1, NC

                 DO 96 IREF = 1, 4

                    RX(IREF+M) = RX(IREF) + CELL * FLOAT ( IX - 1 )
                    RY(IREF+M) = RY(IREF) + CELL * FLOAT ( IY - 1 )
                    RZ(IREF+M) = RZ(IREF) + CELL * FLOAT ( IZ - 1 )

c                    EX(IREF+M) = EX(IREF)
c                    EY(IREF+M) = EY(IREF)
c                    EZ(IREF+M) = EZ(IREF)

96               CONTINUE

                 M = M + 4

97            CONTINUE

98         CONTINUE

99      CONTINUE

C    ** SHIFT CENTRE OF BOX TO THE ORIGIN **

        DO 100 I = 1, N

           RX(I) = RX(I) - 0.5*BOXL
           RY(I) = RY(I) - 0.5*BOXL
           RZ(I) = RZ(I) - 0.5*BOXL
c         write(7,*)rx(i),ry(i),rz(i)
100     CONTINUE

        RETURN
        END


        SUBROUTINE CUBICPB ( BOX )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** PERIODIC BOUNDARY CONDITION FOR CUBIC BOX OF LENGTH BOX       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c        parameter (n =864)
	parameter (N=1372)
c	parameter (n=499)
c        PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

        REAL*8        BOX,BOXINV
        REAL*8        RX(N), RY(N), RZ(N)

        INTEGER     I

C    *******************************************************************
        BOXINV=1.0 / BOX

        DO 100 I = 1, N

          RX(I) = RX(I) - ANINT ( RX(I) * BOXINV ) * BOX
          RY(I) = RY(I) - ANINT ( RY(I) * BOXINV ) * BOX
          RZ(I) = RZ(I) - ANINT ( RZ(I) * BOXINV ) * BOX

100     CONTINUE

        RETURN
        END


        SUBROUTINE MAXWELL ( TEMP )

        COMMON / BLOCK1 / VX, VY, VZ

C    *******************************************************************
C    ** TRANSLATIONAL VELOCITIES FROM MAXWELL-BOLTZMANN DISTRIBUTION  **
C    **                                                               **
C    ** THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND (UNIT) MASS.**
C    ** THIS ROUTINE IS GENERAL, AND CAN BE USED FOR ATOMS, LINEAR    **
C    ** MOLECULES, AND NON-LINEAR MOLECULES.                          **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           **
C    **                                                               **
C    ** REAL FUNCTION GAUSS ( DUMMY )                                 **
C    **    RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             **
C    **    DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.             **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c         parameter (n =864)
          parameter (N=1372)
c        PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

        REAL*8        VX(N), VY(N), VZ(N)
        REAL*8        TEMPV,SC2,SC,EKIN

        REAL*8        U1, U2, V1, V2, S, R
        REAL*8        RTEMP, SUMX, SUMY, SUMZ
        REAL*8        RANF
        INTEGER     I

C    *******************************************************************


        DO 100 I = 1, N, 2

 1         U1 = RANF()
           U2 = RANF()

           V1 = 2.0 * U1 - 1.0
           V2 = 2.0 * U2 - 1.0
           S  = V1 * V1 + V2 * V2

           IF ( S .GE. 1.0 ) GO TO 1

           R = -2.0 * DLOG ( S ) / S

           VX ( I )   = V1 * SQRT ( R )
           VX ( I+1 ) = V2 * SQRT( R )

           VY ( I )   = V1 * SQRT ( R )
           VY ( I+1 ) = V2 * SQRT ( R )

           VZ ( I )   = V1 * SQRT ( R )
           VZ ( I+1 ) = V2 * SQRT ( R )

100     CONTINUE

C    ** REMOVE NET MOMENTUM **

        SUMX = 0.0
        SUMY = 0.0
        SUMZ = 0.0

        DO 200 I = 1, N

           SUMX = SUMX + VX(I)
           SUMY = SUMY + VY(I)
           SUMZ = SUMZ + VZ(I)

200     CONTINUE

        SUMX = SUMX / FLOAT ( N )
        SUMY = SUMY / FLOAT ( N )
        SUMZ = SUMZ / FLOAT ( N )

           EKIN=0.0
        DO 300 I = 1, N

           VX(I) = VX(I) - SUMX
           VY(I) = VY(I) - SUMY
           VZ(I) = VZ(I) - SUMZ

           EKIN=EKIN+(VX(I)*VX(I)+VY(I)*VY(I)+VZ(I)*VZ(I))

300     CONTINUE

C           SC=1.0
           EKIN=0.5*EKIN
           WRITE(*,*)'VELOCITY SCALING'
           TEMPV=2.0*EKIN/(3.0*(N-1))
           SC2=TEMP/TEMPV
           SC=SQRT(SC2)
           WRITE(*,*)'SCALE FACTOR IN MAXWELL IS', SC

        DO 400 I = 1, N

           VX(I) = VX(I)*SC
           VY(I) = VY(I)*SC
           VZ(I) = VZ(I)*SC


400     CONTINUE

        RETURN
        END


        SUBROUTINE WRITCN2 (  )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / VX, VY, VZ
        COMMON / BLOCK3 / FX, FY, FZ
C        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** ROUTINE TO WRITE OUT FINAL CONFIGURATION TO UNIT 10           **
C    *******************************************************************

        INTEGER     N
c        PARAMETER ( N = 108 )
c        parameter (n =864)
	parameter (N=1372)
c	parameter (N=499)
c        PARAMETER ( N = 256 )
c        PARAMETER ( N = 6912 )
c        PARAMETER ( N = 4000 )

c        CHARACTER   CNFILE*(*)
c        CHARACTER   CNFILE*30
        REAL*8        RX(N), RY(N), RZ(N)
        REAL*8        VX(N), VY(N), VZ(N)
        REAL*8        FX(N), FY(N), FZ(N)

c        INTEGER     CNUNIT
c        PARAMETER ( CNUNIT = 11 )

C    ****************************************************************

c        OPEN  ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'UNKNOWN',
c     :          FORM = 'UNFORMATTED' )
        OPEN  ( UNIT = 11, FILE = 'position.out', STATUS = 'UNKNOWN')
        OPEN  ( UNIT = 12, FILE = 'velocity.out', STATUS = 'UNKNOWN')

c        WRITE ( CNUNIT ) N

c        WRITE ( CNUNIT ) RX, RY, RZ
c        WRITE ( CNUNIT ) VX, VY, VZ

        WRITE (11,*) N
        WRITE (12,*) N

        DO I=1, N
c        WRITE (*,* ) RX, RY, RZ
c        WRITE (11,* ) VX, VY, VZ
        WRITE (11,* ) RX(I), RY(I), RZ(I)
        WRITE (12,* ) VX(I), VY(I), VZ(I)
        ENDDO
c        CLOSE ( UNIT = CNUNIT )

        RETURN
        END
