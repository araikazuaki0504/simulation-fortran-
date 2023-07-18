*********************************************************************
***                                                               ***
***           TWO-DIMENSIONAL PARTICLE ELEMENT METHOD             ***
***              SPHERE MODEL  ......   I-PEM1.F                  ***
***                                                               ***
*********************************************************************
      IMPLICIT REAL*8(A-H, K , M ,O-Z)

      PARAMETER (NI=1000,NJ=13,NC=20000)
      COMMON /CON/  DT,FRI,FRW,E,EW,PO,POW,SO,G,DE,PI
      COMMON /WEP/  RR(NI),WEI(NI),PMI(NI)
      COMMON /CEL/  N,IDX,IDZ,IPZ,W,C,NCL(NC),NNCL(NI)
      COMMON /POS/  X0(NI),Z0(NI),QQ(NI)
      COMMON /VEL/  U0(NI),V0(NI),F0(NI) 
      COMMON /FOR1/ XF(NI),ZF(NI),MF(NI)
      COMMON /FOR2/ EN(NI,NJ),ES(NI,NJ),JE(NI,NJ)
      COMMON /DPM/  U(NI+3),V(NI+3),F(NI+3)       
*
      MAXSTEP=10000000
*
*---- SETTING UP THE FIRST POSITION AND INITIAL CONDITION-------
***** 初期粒子座標の設定 *****
      CALL FPOSIT(RMAX)
***** 物性値，定数の設定 *****
      CALL INMAT
*
*---- INITIAL SETTING-----------------------------------  
***** 接触力の初期化 *****
      CALL INIT
*
      T=0.0D0
*
*---- ITERATION FOR EACH STEP-----------------------------------
***** ステップン繰り返し *****
      DO 120 IT=1,MAXSTEP
        T=T+DT 
*
***** セルのクリア，粒子の格納 *****
        CALL NCEL
***** 粒子数の繰り返し *****
        DO 100 I=1,N
***** 各粒子の合力のクリア *****
          XF(I)=0.0D0 ; ZF(I)=0.0D0 ; MF(I)=0.0D0
*-----CALCULATION OF CONTACT FORCE BETWEEN PARTICLE AND WALL---
***** 境界との接触判定 *****
          CALL WCONT(I)
*-----CALCULATION OF CONTACT FORCE BETWEEN PARTICLES---------
***** 粒子間の接触判定 *****
          CALL PCONT(I,RMAX)
  100   CONTINUE 
*                                           
*-----SUPERPOSITION OF INCREMENTAL DISPLACEMENT ------------
***** 速度，変位増分の計算値より粒子を移動 *****
        CALL NPOSIT(JUDGE,IT)
*-----JUDGEMENT OF STATIC STATE
***** 静止状態の判定 *****
        IF (1-JUDGE) 200,200,199
  199   IF (MOD(IT,10000).EQ.0) THEN
          WRITE(6,1000) T,X0(n),Z0(n),U0(n),V0(n)
        ENDIF
*-----OUTPUT OF GRAPHIC DATA----------------------------------
***** グラフィック用データの出力 *****
        IF (IT.EQ.1 .OR. MOD(IT,10000).EQ.0) THEN
          CALL GFOUT(IT,T,RMAX)
        ENDIF
*
  120 CONTINUE
*
 1000 format('time=',e11.4,' x0(n)=',e11.4,' z0(n)=',e11.4,' u0(n)='
     &,e11.4,' v0(n)=',e11.4)
 1999 format(e11.4,4(',',e11.4))
*
*---------OUTPUT OF BACK UP DATA-----------------------------------
***** バックアップ用データの出力 *****
  200 CALL BFOUT(T,RMAX)
*
      open(50,file='xz-out.csv')
        write(50,*)'time,x0(i-out),z0(i-out),rr(i-out)'
      do 500 I= 1,n
        write(50,5000) t,x0(i),z0(i),rr(i)
  500 continue
        write(50,*)'time,x0(i-out),z0(i-out),rr(i-out),R=R1'
      do 501 I= 1,n
        if(rr(i) .eq. rmax) then
          write(50,5000) t,x0(i),z0(i),rr(i)
        endif
  501 continue
        write(50,*)'time,x0(i-out),z0(i-out),rr(i-out),R=R2'
      do 502 I= 1,n
        if(rr(i) .ne. rmax) then
          write(50,5000) t,x0(i),z0(i),rr(i)
        endif
  502 continue
 5000 format(e11.4,3(',',e11.4))
*
      CLOSE(10) ; CLOSE(11)
*
      STOP
      END 
C  
C******* FPOSIT *******************************************
C******* 初期粒子座標の設定 *****
C                                                             
      SUBROUTINE FPOSIT(RMAX)
      IMPLICIT REAL*8(A-H,K,M,O-Z)
      PARAMETER (NI=1000,NJ=13,NC=20000)
      COMMON /CON/DT,FRI,FRW,E,EW,PO,POW,SO,G,DE,PI
      COMMON /WEP/RR(NI),WEI(NI),PMI(NI)
      COMMON /CEL/N,IDX,IDZ,IPZ,W,C,NCL(NC),NNCL(NI)
      COMMON /POS/X0(NI),Z0(NI),QQ(NI)
***** 乱数発生の引数(II)，粒子径
***** 粒子数(N),セル数(IDX,IDZ),粒子数(IPZ),容器幅(W),セル幅(C)
      II=584287
      R1=1.0D-2 ; R2=5.0D-3
      W =6.006D-1
      IPZ=30
      A =6
C
      RMAX = R1 ; RMIN = R2
      RN   = RMAX+1.0D-5
      IPX=IDINT(W/2.0D0/RN)

      N=0
      DO 20 I=1,IPZ
        IF (MOD(I,2).EQ.0)THEN
           DX=2.0D0*RN ; IP=IPX-1
        ELSE
           DX=RN       ; IP=IPX
        ENDIF
        DO 10 J=1,IP
*****　乱数発生　*****
          CALL RANDOM(II,RU)
          IF (RU.LT.2.0D-1) GOTO 10
          N=N+1
          IF (N.GT.NI) WRITE(6,*) 
     &      'NUMBER OF PARTICLES IS MORE THAN ',NI
          X0(N)=2.0D0*RN*(J-1)+DX
          Z0(N)=2.0D0*RN*(I-1)+RN 
C
          CALL RANDOM(II,RU) 
          IF (RU.LT.5.0D-1) THEN
            RR(N)=R1
          ELSE
            RR(N)=R2
          ENDIF          
   10   CONTINUE
   20 CONTINUE 
***** 粒子数の出力　*****
***** セル幅(C)
      WRITE (6,*) 'NUMBER OF PARTICLES:',N
      C=RMIN*1.35D0
      IDX=IDINT(W/C)+1
      IDZ=IDINT(Z0(N)/C)+10
      IF (IDX*IDZ.GT.NC) THEN
        WRITE(6,*) 'NCL IS OVER!!',IDX*IDZ
        STOP
      ENDIF
      WRITE(6,*) '(IDX,IDZ)=',IDX,IDZ
*
***** 粒子中心座標の初期値 *****
***** 中心座標(X0(NI),Z0(NI)),回転変位(QQ(NI)) *****
      OPEN(20,file='xz-int.csv')
        WRITE(6 ,*) 'Initail center coordinate of particles'
        WRITE(20,*)'No,X0(ini),Z0(ini),rr(ini)'
      DO 111 I= 1,N
        WRITE(6 ,1001) i,X0(I),Z0(I),rr(i)
        WRITE(20,1111) i,X0(I),Z0(I),rr(i)
  111 CONTINUE 
***** 粒子の大きさごとに出力 *****
        WRITE(20,*)'No,X0(ini),Z0(ini),rr(ini),R=R1'
      DO 112 I= 1,N
        if(rr(i) .eq. rmax) then
           WRITE(20,1111) i,X0(I),Z0(I),rr(i)
        endif
  112 CONTINUE 
        WRITE(20,*)'No,X0(ini),Z0(ini),rr(ini),R=R2'
      DO 113 I= 1,N
        if(rr(i) .ne. rmax) then
           WRITE(20,1111) i,X0(I),Z0(I),rr(i)
        endif
  113 CONTINUE 
 1001 format('No=',i3,' X0(ini)=',f10.4,' z0(ini)=',f10.4,' r=',f10.4)
 1111 format(i3,3(',',e11.4))
      CLOSE(20)
*
      RETURN
      END
C
C******* INMAT ******************************************
C******* 物性値，定数の設定 *****
C                                              
      SUBROUTINE INMAT
      IMPLICIT REAL*8(A-H,K,M,O-Z)
      PARAMETER (NI=1000,NJ=13,NC=20000)
      COMMON /CON/DT,FRI,FRW,E,EW,PO,POW,SO,G,DE,PI
      COMMON /WEP/RR(NI),WEI(NI),PMI(NI)
      COMMON /CEL/N,IDX,IDZ,IPZ,W,C,NCL(NC),NNCL(NI)
***** 差分時間(DT),摩擦係数(FRI,FRW),弾性係数(E,EW),ﾎﾟｱｿﾝ比(PO,POW)
***** 弾性係数比(SO),G,粒子密度(DE),PI
      OPEN(15,FILE='datain.dat',STATUS='UNKNOWN')
      READ(15,*)
      READ(15,151) DT
      READ(15,152) FRI,FRW
      READ(15,152) E  ,EW 
      READ(15,152) PO ,POW
  151 FORMAT(D10.2)
  152 FORMAT(D10.2,1X,D10.2)
C      WRITE(6,*)' DT =',DT
C      WRITE(6,*)' FRI=',FRI,' FRW=',FRW
C      WRITE(6,*)' E  =',E  ,' ET =',EW 
C      WRITE(6,*)' PO =',PO ,' POW=',POW
C      DT = 5.0D-7
C      FRI= 2.5D-1 ; FRW= 1.7D-1
C      E  = 4.9D10 ; EW = 3.9D9
C      PO = 2.3D-1 ; POW= 2.5D-1

      SO=1.0D0/2.0D0/(1.0D0+PO)
      G  = 9.80665D0
      DE = 2.48D3

      PI = 3.14159D0
***** 粒子質量(WEI(NI)),慣性M(PMI(NI))
      DO 10 I=1,N
        WEI(I)=4.0D0/3.0D0*PI*RR(I)*RR(I)*RR(I)*DE
        PMI(I)=8.0D0/15.0D0*DE*PI*(RR(I))**5
   10 CONTINUE
C
      RETURN
      END
C  
C******* INIT *******************************************
C******* 接触力の初期化 *****
C                                                        
      SUBROUTINE INIT
      IMPLICIT REAL*8(A-H,K,M,O-Z)
      PARAMETER (NI=1000,NJ=13,NC=20000)
      COMMON /CEL/N,IDX,IDZ,IPZ,W,C,NCL(NC),NNCL(NI)     
      COMMON /FOR2/EN(NI,NJ),ES(NI,NJ),JE(NI,NJ)
      COMMON /DPM/U(NI+3),V(NI+3),F(NI+3)
C
      DO 3 I=1,N 
        DO 2 J=1,nj
          EN(I,J)=0.0D0 ; ES(I,J)=0.0D0  
    2   CONTINUE    
        DO 4 IJ=1,nj
         JE(I,IJ)=0 
    4   CONTINUE    
    3 CONTINUE  
*
      DO 5 I=1,N+3
        U(I)=0.0D0 ; V(I)=0.0D0 ; F(I)=0.0D0
    5 CONTINUE                                    
      RETURN                                             
      END
C  
C******* NCEL *******************************************
C******* セルのクリア，粒子の格納 *****
C                                                        
      SUBROUTINE NCEL
      IMPLICIT REAL*8(A-H,K,M,O-Z)
      PARAMETER (NI=1000,NJ=13,NC=20000)
      COMMON /CEL/N,IDX,IDZ,IPZ,W,C,NCL(NC),NNCL(NI)
      COMMON /POS/X0(NI),Z0(NI),QQ(NI)
      COMMON /VEL/U0(NI),V0(NI),F0(NI)      
      COMMON /FOR1/XF(NI),ZF(NI),MF(NI)        
      COMMON /FOR2/EN(NI,NJ),ES(NI,NJ),JE(NI,NJ)     
      COMMON /DPM/U(NI+3),V(NI+3),F(NI+3) 
C
***** セル配列(NCL(NC)),セル番号(NNCL(NI))
      DO 10 IB=1,(IDX*IDZ)
        NCL(IB)=0         
   10 CONTINUE            
      DO 20 I=1,N         
         NNCL(I)=0
         IB=IDINT(Z0(I)/C)*IDX+IDINT(X0(I)/C)+1
         NCL(IB)=I
         NNCL(I)=IB
   20 CONTINUE     
      RETURN       
      END
C  
C******* WCONT *******************************************
C******* 境界との接触判定 *****
C                                                             
      SUBROUTINE WCONT(I)
      IMPLICIT REAL*8(A-H,K,M,O-Z)
      PARAMETER (NI=1000,NJ=13,NC=20000)
      COMMON /WEP/RR(NI),WEI(NI),PMI(NI)
      COMMON /CEL/N,IDX,IDZ,IPZ,W,C,NCL(NC),NNCL(NI)
      COMMON /POS/X0(NI),Z0(NI),QQ(NI)
      COMMON /VEL/U0(NI),V0(NI),F0(NI)      
      COMMON /FOR1/XF(NI),ZF(NI),MF(NI)        
      COMMON /FOR2/EN(NI,NJ),ES(NI,NJ),JE(NI,NJ)     
      COMMON /DPM/U(NI+3),V(NI+3),F(NI+3)  
C
      XI=X0(I) ; ZI=Z0(I) ; RI=RR(I)
C
C  --- LEFT WALL ---
C
      JK=11              
      J=N+1              
      IF (XI.LT.RI) THEN 
        AS=0.0D0         
        AC=-1.0D0
        GAP=DABS(XI)
        JE(I,JK)=N+1     
        CALL ACTF(I,J,JK,AS,AC,GAP) 
      ELSE                          
        EN(I,JK)=0.0D0 ; ES(I,JK)=0.0D0 
        JE(I,JK)=0
      ENDIF             
C 
C  --- UNDER WALL ---
C
      JK=12             
      J =N+2            
      IF (ZI.LT.RI) THEN
        AS=-1.0D0       
        AC= 0.0D0
        GAP=DABS(ZI)    
        JE(I,JK)=N+2
        CALL ACTF(I,J,JK,AS,AC,GAP)
      ELSE
        EN(I,JK)=0.0D0 ; ES(I,JK)=0.0D0
        JE(I,JK)=0
      ENDIF
                             
C
C   --- RIGHT WALL ---
C
      JK=nj
      J=N+3
      IF (XI+RI.GT.W) THEN
        AS=0.0D0
        AC=1.0D0
        GAP=DABS(XI-W)
        JE(I,JK)=N+3
        CALL ACTF(I,J,JK,AS,AC,GAP)
      ELSE
        EN(I,JK)=0.0D0 ; ES(I,JK)=0.0D0
        JE(I,JK)=0
      ENDIF
      RETURN                                                  
      END
C
C  
C******* PCONT *******************************************
C******* 粒子間の接触判定 *****
C                                                             
      SUBROUTINE PCONT(I,RMAX)
      IMPLICIT REAL*8(A-H,K,M,O-Z)
      PARAMETER (NI=1000,NJ=13,NC=20000)
      COMMON /WEP/RR(NI),WEI(NI),PMI(NI)
      COMMON /CEL/N,IDX,IDZ,IPZ,W,C,NCL(NC),NNCL(NI)
      COMMON /POS/X0(NI),Z0(NI),QQ(NI)
      COMMON /VEL/U0(NI),V0(NI),F0(NI)      
      COMMON /FOR1/XF(NI),ZF(NI),MF(NI)        
      COMMON /FOR2/EN(NI,NJ),ES(NI,NJ),JE(NI,NJ)     
      COMMON /DPM/U(NI+3),V(NI+3),F(NI+3) 
C
      XI=X0(I) ; ZI=Z0(I) ; RI=RR(I)
C
      LUP=IDINT((ZI+2.0D0*RMAX)/C)
      LUN=IDINT((ZI-2.0D0*RMAX)/C)
      LLF=IDINT((XI-2.0D0*RMAX)/C)
      LRG=IDINT((XI+2.0D0*RMAX)/C)
      IF (LUP.GE.IDZ) LUP=IDZ-1   
      IF (LUN.LT.0)   LUN=0       
      IF (LLF.LT.0)   LLF=0 
      IF (LRG.GE.IDX) LRG=IDX-1   
C
      IF(LUP.LE.LUN) THEN
         WRITE(6,*) I,'LUP=',LUP,'LUN=',LUN,'C=',C,
     &              'RMAX=',RMAX,XI,ZI,IDZ

      ENDIF
C
      DO 90 LZ=LUN,LUP              
        DO 80 LX=LLF,LRG            
         IB=LZ*IDX+LX+1             
         J=NCL(IB)                  
         IF (J.LE.0) GOTO 80
         IF (J.EQ.I) GOTO 80        
         DO 11 JJ=1,10              
           IF (JE(I,JJ).EQ.J) THEN  
             JK=JJ                  
             GOTO 70                
           END IF                   
   11   CONTINUE                    
        DO 12 JJ=1,10               
          IF (JE(I,JJ).EQ.0) THEN   
            JK=JJ                   
            JE(I,JJ)=J              
            GOTO 70                 
          END IF                    
   12   CONTINUE 
   70   XJ=X0(J)                    
        ZJ=Z0(J)
        RJ=RR(J)                    
        GAP=DSQRT((XI-XJ)*(XI-XJ)+(ZI-ZJ)*(ZI-ZJ))               
        IF (GAP.LT.(RI+RJ)) THEN
          IF (I.GT.J) THEN  
            AC=(XJ-XI)/(GAP)               
            AS=(ZJ-ZI)/(GAP)               
            J0=0
           DO 555 JJ=1,10
            IF (JE(J,JJ).EQ.I) THEN
                J0=JJ
                GOTO 554
           ENDIF
  555       CONTINUE   
  554       CALL ACTF(I,J,JK,AS,AC,GAP)    
            EN(J,J0)=EN(I,JK)
            ES(J,J0)=ES(I,JK)
            J0=0
         ENDIF
        ELSE            
   85     EN(I,JK)=0.0D0 ; ES(I,JK)=0.0D0
          JE(I,JK)=0  
        ENDIF                              
   80   CONTINUE                           
   90 CONTINUE  
      RETURN                               
      END
C  
C******* NPOSIT *******************************************
C******* 速度，変位増分の計算値より粒子を移動 *****
C                              
      SUBROUTINE NPOSIT(JUDGE,IT)
      IMPLICIT REAL*8(A-H,K,M,O-Z)
      PARAMETER (NI=1000,NJ=13,NC=20000)
      COMMON /CON/DT,FRI,FRW,E,EW,PO,POW,SO,G,DE,PI
      COMMON /WEP/RR(NI),WEI(NI),PMI(NI)
      COMMON /CEL/N,IDX,IDZ,IPZ,W,C,NCL(NC),NNCL(NI)
      COMMON /POS/X0(NI),Z0(NI),QQ(NI)
      COMMON /VEL/U0(NI),V0(NI),F0(NI)      
      COMMON /FOR1/XF(NI),ZF(NI),MF(NI)        
      COMMON /DPM/U(NI+3),V(NI+3),F(NI+3) 
C
      SUM=0.0D0                                            
      DO 110 I=1,N
        BV0=V0(I) ; BU0=U0(I) ; BF0=F0(I)
        V0(I)=V0(I)+(ZF(I)-WEI(I)*G)/WEI(I)*DT
        U0(I)=U0(I)+XF(I)/WEI(I)*DT           
        F0(I)=F0(I)+MF(I)/PMI(I)*DT           
        V(I)=(1.5D0*V0(I)-0.5D0*BV0)*DT       
        U(I)=(1.5D0*U0(I)-0.5D0*BU0)*DT       
        F(I)=(1.5D0*F0(I)-0.5D0*BF0)*DT       
        Z0(I)=Z0(I)+V(I)                      
        X0(I)=X0(I)+U(I)                      
        QQ(I)=QQ(I)+F(I)
        SUM=DABS(U(I))+DABS(V(I))+SUM
  110 CONTINUE                                
      AV=SUM/DFLOAT(N)/2.0D0 
*      IF (MOD(IT,10000).EQ.0) THEN
      IF (MOD(IT,100000).EQ.0) THEN
         WRITE(6,*) 'AVERAGE= ',AV
      ENDIF                                
      IF (AV.LT.(DT*DT*G*1.0D-2)) THEN
        JUDGE=1 
      ELSE
        JUDGE=0
      ENDIF
            
      RETURN                                  
      END  
C  
C******* ACTF *******************************************   
C                                                             
      SUBROUTINE ACTF(I,J,JK,AS,AC,GAP)              
      IMPLICIT REAL*8(A-H,K,M,O-Z)
      PARAMETER (NI=1000,NJ=13,NC=20000)
      COMMON /CON/DT,FRI,FRW,E,EW,PO,POW,SO,G,DE,PI
      COMMON /WEP/RR(NI),WEI(NI),PMI(NI)
      COMMON /CEL/N,IDX,IDZ,IPZ,W,C,NCL(NC),NNCL(NI)
      COMMON /POS/X0(NI),Z0(NI),QQ(NI)
      COMMON /VEL/U0(NI),V0(NI),F0(NI)      
      COMMON /FOR1/XF(NI),ZF(NI),MF(NI)        
      COMMON /FOR2/EN(NI,NJ),ES(NI,NJ),JE(NI,NJ)            
      COMMON /DPM/U(NI+3),V(NI+3),F(NI+3) 
C
      RI=RR(I)
      IF ((J-N).LE.0) THEN 
        RJ=RR(J)
        DIS=RI+RJ-GAP
        WEI3=2.0D0*WEI(I)*WEI(J)/(WEI(I)+WEI(J)) 
      ELSE
        RJ=0.0D0
        DIS=RI-GAP
        WEI3=WEI(I)
      ENDIF 
      ENN=EN(I,JK)
      IF (ENN.LE.0.0D0) THEN
        ENN=1.0D-3
      ENDIF 
      IF ((J-N).LE.0.0D0) THEN                                       
        B1=(3.0D0/2.0D0/E*RI*RJ/(RI+RJ)*(1.0D0-PO*PO)
     &  *ENN)**(1.0D0/3.0D0)
        KNN=2.0D0/3.0D0*B1*E/(1.0D0-PO*PO)
C        KSS=KNN*SO
        KSS=8.0D0*B1*SO*E/2.0D0/(2.0D0-PO)
        VNN=DSQRT(2.0D0*WEI3*KNN)
        VSS=VNN*DSQRT(KSS/KNN)
      ELSE
        B1=((3.0D0/4.0D0*RI*((1.0D0-PO*PO)/E+(1.0D0-POW*POW)/EW))
     &  *ENN)**(1.0D0/3.0D0)
        KNN=4.0D0/3.0D0*B1*E*EW/((1.0D0-PO*PO)*EW+(1.0D0-POW*POW)*E)
C        KSS=KNN*SO
        KSS=8.0D0*B1*SO*E/(2.0D0-PO)
        VNN=DSQRT(2.0D0*WEI3*KNN)
        VSS=VNN*DSQRT(KSS/KNN)
      ENDIF
      DDT=1.0D-1*DSQRT(WEI3/KNN)
      IF (DDT.LT.1.0D-6) WRITE (6,*) 'OVER!! DDT=',DDT,I,J,JK,KNN,WEI(I)
      UN=(U(I)-U(J))*AC+(V(I)-V(J))*AS                        
      US=-(U(I)-U(J))*AS+(V(I)-V(J))*AC+(RI*F(I)+RJ*F(J))
C----------------------------------------------------------------------
      IF(EN(I,JK).EQ.0.0D0) THEN
        IF (UN.NE.0.0D0) US=US*DIS/UN
        UN=DIS
      ENDIF
C----------------------------------------------------------------------
C      EN(I,JK)=EN(I,JK)+KNN*UN 
      EN(I,JK)=KNN*DIS**1.5D0
      ES(I,JK)=ES(I,JK)+KSS*US
	
      DN=VNN*UN/DT                       
      DS=VSS*US/DT
      IF (EN(I,JK).LT.0.0D0) THEN        
        EN(I,JK)=0.0D0 ; ES(I,JK)=0.0D0  
        DN=0.0D0       ; DS=0.0D0        
        JE(I,JK)=0                       
      RETURN                             
      ELSEIF ((J-N).LE.0) THEN           
        FRC=FRI                          
      ELSE                               
        FRC=FRW
      ENDIF                              
      IF ((DABS(ES(I,JK))-FRC*EN(I,JK)).GT.0.0D0) THEN
        ES(I,JK)=FRC*DSIGN(EN(I,JK),ES(I,JK))         
        DS=0.0D0
      ENDIF              
        HN=EN(I,JK)+DN           
        HS=ES(I,JK)+DS           
        XF(I)=-HN*AC+HS*AS+XF(I) 
        ZF(I)=-HN*AS-HS*AC+ZF(I) 
        MF(I)=MF(I)-RI*HS 
        IF (JK.LE.10) THEN       
          XF(J)=HN*AC-HS*AS+XF(J)
          ZF(J)=HN*AS+HS*AC+ZF(J)
          MF(J)=MF(J)-RJ*HS
        ENDIF 
 
      RETURN                     
      END  
C   
C******* GFOUT *******************************************
C******* グラフィック用データの出力 *****
C                                                             
      SUBROUTINE GFOUT(IT,T,RMAX)
      IMPLICIT REAL*8(A-H,K,M,O-Z)
      PARAMETER (NI=1000,NJ=13,NC=20000)
      COMMON /WEP/RR(NI),WEI(NI),PMI(NI)
      COMMON /CEL/N,IDX,IDZ,IPZ,W,C,NCL(NC),NNCL(NI)
      COMMON /POS/X0(NI),Z0(NI),QQ(NI)
      COMMON /VEL/U0(NI),V0(NI),F0(NI)            
      COMMON /FOR2/EN(NI,NJ),ES(NI,NJ),JE(NI,NJ) 
C
      IF (IT.EQ.1) THEN
c     OPEN(10,FILE='GRAPH11.0D')
c     OPEN(11,FILE='GRAPH21.0D')
c       WRITE(10,*) N,T,W,RMAX
      OPEN(40,file='GRAPHxz.csv')
      OPEN(41,file='GRAPHuv.csv')
        WRITE(40,*   ) 'N,T,W,RMAX'
        WRITE(40,4000)  N,T,W,RMAX
      write(40,*)'it,x0( 1),z0(  1),rr(  1),x0(100),z0(100),rr(100)
     &,x0(200),z0(200),rr(200),x0(300),z0(300),rr(300)
     &,x0(400),z0(400),rr(400),x0(500),z0(500),rr(500)
     &,x0(600),z0(600),rr(600),x0(  n),z0(  n),rr(  n)'
***** (60,FILE='tecpxz.dat') *****
      OPEN(60,file='tecpxz.dat')
      WRITE(60,*) ' TITLE=" COORDINATE SYSTEM" '
      WRITE(60,*) ' VARIABLES="X [-]","ZX [-]" '

      ELSE
c       WRITE(10,*) (SNGL(X0(I)),SNGL(Z0(I)),SNGL(RR(I)),I=1,N)
c       WRITE(10,*) (SNGL(U0(I)),SNGL(V0(I)),SNGL(F0(I)),I=1,N)

c       WRITE(11,*) ((SNGL(ES(I,J)),SNGL(EN(I,J)),J=1,nj),I=1,N)
c       WRITE(11,*) ((JE(I,J),J=1,nj),I=1,N)

      write(40,4001) it,x0( 1),z0(  1),rr(  1),x0(100),z0(100),rr(100)
     &,x0(200),z0(200),rr(200),x0(300),z0(300),rr(300)
     &,x0(400),z0(400),rr(400),x0(500),z0(500),rr(500)
     &,x0(600),z0(600),rr(600),x0(  n),z0(  n),rr(  n)
*
      write(41,4001) it,u0( 1),v0(  1),rr(  1),u0(100),v0(100),f0(100)
     &,u0(200),v0(200),f0(200),u0(300),v0(300),f0(300)
     &,u0(400),v0(400),f0(400),u0(500),v0(500),f0(500)
     &,u0(600),v0(600),f0(600),u0(  n),v0(  n),f0(  n)

      WRITE(60,*)'ZONE T="n",I=',N,',C=blue,F=POINT'
      do 600 i=1,n
      write(60,6001) x0(i),z0(i)
  600 continue
      ENDIF
 4000 format(i10, 4(',',e11.4))
 4001 format(i10,24(',',e11.4))
 6001 format(e11.4,2x,e11.4)
      RETURN 
      END

C
C******* BFOUT ******************************************* 
C******* バックアップ用データの出力 *****
C                                                             
      SUBROUTINE BFOUT(T,RMAX)
      IMPLICIT REAL*8(A-H,K,M,O-Z)                          
      PARAMETER (NI=1000,NJ=13,NC=20000)
      COMMON /CON/DT,FRI,FRW,E,EW,PO,POW,SO,G,DE,PI
      COMMON /WEP/RR(NI),WEI(NI),PMI(NI)
      COMMON /CEL/N,IDX,IDZ,IPZ,W,C,NCL(NC),NNCL(NI)
      COMMON /POS/X0(NI),Z0(NI),QQ(NI)
      COMMON /VEL/U0(NI),V0(NI),F0(NI)      
      COMMON /FOR1/XF(NI),ZF(NI),MF(NI)        
      COMMON /FOR2/EN(NI,NJ),ES(NI,NJ),JE(NI,NJ)     
      COMMON /DPM/U(NI+3),V(NI+3),F(NI+3)       
C      
***** 配列数(NI,NJ),NC *****
***** 粒子数(N),セル数(IDX,IDZ),粒子数(IPZ) ****
***** 最大粒子径(Rmax),時刻(T),容器幅(W),セル幅(C),時間刻(DT) *****
***** 粒子密度(DE),摩擦係数(FRI,FRW) *****
***** 弾性係数(E,EW),ﾎﾟｱｿﾝ比(PO,POW),弾性係数比(SO) *****
***** 粒子質量(WEI(NI)),慣性M(PMI(NI)) *****
***** 中心座標(X0(NI),Z0(NI)),粒子半径(RR(NI)) *****
***** 並進変位(U(NI+3),V(NI+3)) *****
***** 速度(U0(NI),V0(NI),F0(NI)) *****
***** 合力(XF(NI),ZF(NI)),ﾓｰﾒﾝﾄ(MF(NI)) *****
***** 接触力(EN(NI,NJ),ES(NI,NJ)) *****
***** 接触点番号(JE(NI,NJ)) *****
      OPEN(13,FILE='BACK100.0D')
      WRITE(13,*) N,IDX,IDZ,IPZ
      WRITE(13,*) RMAX,T,W,C,DT
      WRITE(13,*) DE,FRI,FRW,G,PI
      WRITE(13,*) E,EW,PO,POW,SO 
      WRITE(13,*) (WEI(I),PMI(I),I=1,N) 
      WRITE(13,*) (X0(I),Z0(I),RR(I),I=1,N)
      WRITE(13,*) (U(I),V(I),F(I),I=1,N)
      WRITE(13,*) (U0(I),V0(I),F0(I),I=1,N)
      WRITE(13,*) ((ES(I,J),EN(I,J),J=1,13),I=1,N)
      WRITE(13,*) ((JE(I,J),J=1,13),I=1,N)
      CLOSE(13)

      RETURN ; END
C
C******* RANDOM ******************************************
C******* 乱数の発生 *****
      SUBROUTINE RANDOM(II,RU)
      IMPLICIT REAL*8(A-H,K,M,O-Z)
C
      II=II*65539
      IF (II.LT.0) II=(II+2147483647)+1
      RU=II*0.4656613D-9
      RETURN ; END
