      SUBROUTINE CP_RCVQBCML2NS(MX_ML,MY_ML,IEAS,IWES,JSOU,JNOR,BUF,
     $                          ZBED_ML,QBXC_ML,QBYC_ML)
C----------------------------------------------------------------------
C     MLの境界におけるセル中心掃流砂量をを受信する。
C----------------------------------------------------------------------
C
      use mod_comm,only: comm_model
      IMPLICIT NONE
C
      INCLUDE  'mpif.h'
      INCLUDE  'CONNEC.h'
C
      INTEGER,INTENT(INOUT)::MX_ML,MY_ML
      INTEGER,INTENT(INOUT)::IEAS,IWES,JSOU,JNOR
C
      REAL(8),INTENT(INOUT)::ZBED_ML(IWES-1:IEAS+1,JSOU-1:JNOR+1)
      REAL(8),INTENT(INOUT)::QBXC_ML(IWES-1:IEAS+1,JSOU-1:JNOR+1)
      REAL(8),INTENT(INOUT)::QBYC_ML(IWES-1:IEAS+1,JSOU-1:JNOR+1)
      REAL(8),INTENT(INOUT)::BUF(*)
C
      INTEGER::ISTAT(MPI_STATUS_SIZE)
C
      INTEGER::I,IERROR,IPARNT,IREQ,J,NCOUNT,NDATA,N2D
C
      IPARNT=IPECON(2,NRANK+1)
C
      N2D = (IEAS-IWES+3)*(JNOR-JSOU+3)
      NDATA=4*(JNOR-JSOU+3)+4*MAX(0,IEAS-IWES-1)
C
C ... QBXC_MLを受信する。
      CALL ZERCLR(QBXC_ML,N2D,0.0D0)
C
      CALL MPI_IRECV(BUF,NDATA,MPI_DOUBLE_PRECISION,IPARNT,
     &               MPI_ANY_TAG,comm_model,IREQ,IERROR)
      CALL MPI_WAIT(IREQ,ISTAT,IERROR)
C
      NCOUNT=0
C
      DO 100 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      QBXC_ML(IWES-1,J)=BUF(NCOUNT)
  100 CONTINUE
C
      DO 110 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      QBXC_ML(IWES,J)=BUF(NCOUNT)
  110 CONTINUE
C
      DO 120 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      QBXC_ML(IEAS,J)=BUF(NCOUNT)
  120 CONTINUE
C
      DO 130 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      QBXC_ML(IEAS+1,J)=BUF(NCOUNT)
  130 CONTINUE
C
      DO 140 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      QBXC_ML(I,JSOU-1)=BUF(NCOUNT)
  140 CONTINUE
C
      DO 150 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      QBXC_ML(I,JSOU)=BUF(NCOUNT)
  150 CONTINUE
C
      DO 160 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      QBXC_ML(I,JNOR)=BUF(NCOUNT)
  160 CONTINUE
C
      DO 170 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      QBXC_ML(I,JNOR+1)=BUF(NCOUNT)
  170 CONTINUE
C
C
C ... QBYC_MLを受信する。
      CALL ZERCLR(QBYC_ML,N2D,0.0D0)
C
      CALL MPI_IRECV(BUF,NDATA,MPI_DOUBLE_PRECISION,IPARNT,
     &               MPI_ANY_TAG,comm_model,IREQ,IERROR)
      CALL MPI_WAIT(IREQ,ISTAT,IERROR)
C
      NCOUNT=0
C
      DO 200 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      QBYC_ML(IWES-1,J)=BUF(NCOUNT)
  200 CONTINUE
C
      DO 210 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      QBYC_ML(IWES,J)=BUF(NCOUNT)
  210 CONTINUE
C
      DO 220 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      QBYC_ML(IEAS,J)=BUF(NCOUNT)
  220 CONTINUE
C
      DO 230 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      QBYC_ML(IEAS+1,J)=BUF(NCOUNT)
  230 CONTINUE
C
      DO 240 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      QBYC_ML(I,JSOU-1)=BUF(NCOUNT)
  240 CONTINUE
C
      DO 250 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      QBYC_ML(I,JSOU)=BUF(NCOUNT)
  250 CONTINUE
C
      DO 260 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      QBYC_ML(I,JNOR)=BUF(NCOUNT)
  260 CONTINUE
C
      DO 270 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      QBYC_ML(I,JNOR+1)=BUF(NCOUNT)
  270 CONTINUE
C
C
C ... ZBED_MLを受信する。
      CALL ZERCLR(ZBED_ML,N2D,0.0D0)
C
      CALL MPI_IRECV(BUF,NDATA,MPI_DOUBLE_PRECISION,IPARNT,
     &               MPI_ANY_TAG,comm_model,IREQ,IERROR)
      CALL MPI_WAIT(IREQ,ISTAT,IERROR)
C
      NCOUNT=0
C
      DO 300 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      ZBED_ML(IWES-1,J)=BUF(NCOUNT)
  300 CONTINUE
C
      DO 310 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      ZBED_ML(IWES,J)=BUF(NCOUNT)
  310 CONTINUE
C
      DO 320 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      ZBED_ML(IEAS,J)=BUF(NCOUNT)
  320 CONTINUE
C
      DO 330 J=JSOU-1,JNOR+1
      NCOUNT=NCOUNT+1
      ZBED_ML(IEAS+1,J)=BUF(NCOUNT)
  330 CONTINUE
C
      DO 340 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      ZBED_ML(I,JSOU-1)=BUF(NCOUNT)
  340 CONTINUE
C
      DO 350 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      ZBED_ML(I,JSOU)=BUF(NCOUNT)
  350 CONTINUE
C
      DO 360 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      ZBED_ML(I,JNOR)=BUF(NCOUNT)
  360 CONTINUE
C
      DO 370 I=IWES+1,IEAS-1
      NCOUNT=NCOUNT+1
      ZBED_ML(I,JNOR+1)=BUF(NCOUNT)
  370 CONTINUE
C
      RETURN
      END
