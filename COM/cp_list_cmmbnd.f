      SUBROUTINE CP_LIST_CMMBND(IERROR)
C
C 1. MAKE UP INDEX & LIST VECTOR FOR COMMUNICATION BOUNDARY AREA
C
C [GLOBAL ENTITIES]
C      USE CP_MODULE_DC1,ONLY : MKLIST
      USE CP_MODULE_INDCMM,ONLY : MKINDEX
      INCLUDE 'CONNEC.h'
      INCLUDE 'GLOBAL.h'
      INCLUDE 'FILE.h'
C
C [INTENT(IN)]
C
C [INTENT(OUT)]
      INTEGER  IERROR
C
C [LOCAL ENTITIES]
      INTEGER  INDDOM(6,0:NPROC-1)
      INTEGER  LX,LY,LZ
      INTEGER  MYRANK
      INTEGER  IPRD,JPRD,KPRD
      INTEGER  I,J,K,IERR,IDUM1
C
C
      IERROR=0
      IF( NPROC.LT.2 ) RETURN
C
      LX=1
      LY=1
      LZ=1
C
      IPRD=0
      JPRD=0
      KPRD=0
c       write(16,*) 'start'
cc     CALL MPI_COMM_RANK(MYRANK,CHILDCOMM,IERR)
       CALL MPI_COMM_RANK(CHILDCOMM,MYRANK,IERR)
c       write(16,*) 'mpi_comm_rank'
C
C-< 2. MAKE UP LIST VECTOR >-
C
       INDDOM=INDCOM(1:6,1:NPROC)-1
cccc     INDDOM=INDCOM-1
          write(16,*) inddom
C
c      CALL MKLIST
c     &     (MYRANK,NPROC,IPRD,JPRD,KPRD
c     &     ,INDDOM,IERROR)
c       write(16,*) 'mklist'
      IF( IERROR.NE.0 ) GOTO 9999
C
      CALL MKINDEX
     &     (MYRANK,NPROC,IPRD,JPRD,KPRD,LX,LY,LZ
     &     ,INDDOM,IERROR)
      IF( IERROR.NE.0 ) GOTO 9999
C
      RETURN
 9999 CONTINUE
      WRITE(LP,*) '(cp_list_cmmbnd)'
      IERROR=1
      END
