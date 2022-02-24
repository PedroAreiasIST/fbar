
!-------------------
!*** MLS PART
!*** MESHLESS PART
!-------------------
!------------------------------
!*** GENERAL TRIANGULAR SOLVE
!*** MULTIPLE RIGHT-HAND-SIDES
!*** CHK0
!------------------------------
  SUBROUTINE MLS_GENTRIANGSOLVE(UPPER,NRHS,N,R,B,X)
    IMPLICIT REAL(8) (a-h,o-z)
    LOGICAL::UPPER
    INTEGER::N
    REAL(8),DIMENSION(N,N)::R
    REAL(8),DIMENSION(N,NRHS)::B,X
    IF(UPPER)THEN
!-------------------------------------------
!*** SOLVE R.X=B WHERE R IS UPPER TRIANGULAR
!-------------------------------------------
       DO IR=1,NRHS
          DO ID=N,1,-1
             X(ID,IR)=B(ID,IR)
             DO JD=ID+1,N
                X(ID,IR)=X(ID,IR)-R(ID,JD)*X(JD,IR)
             END DO
             X(ID,IR)=X(ID,IR)/R(ID,ID)
          END DO
       END DO
    ELSE
!---------------------------------------------
!*** SOLVE R^T.X=B WHERE R IS UPPER TRIANGULAR
!---------------------------------------------
       DO IR=1,NRHS
          DO ID=1,N
             X(ID,IR)=B(ID,IR)
             DO JD=1,ID-1
                X(ID,IR)=X(ID,IR)-R(JD,ID)*X(JD,IR)
             END DO
             X(ID,IR)=X(ID,IR)/R(ID,ID)
          END DO
       END DO
    END IF
  END SUBROUTINE MLS_GENTRIANGSOLVE

!--------------------------
!*** PERFORMS A NAIVE QR
!*** DECOMPOSITION
!*** CHK0
!--------------------------
  SUBROUTINE MLS_PROJQR(N,U,A,UA)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER::N
    REAL(8),DIMENSION(N)::U,A,UA
    UA=U*DOTPROD(N,U,A)/DOTPROD(N,U,U)
  END SUBROUTINE MLS_PROJQR

!--------------------------
!*** PERFORMS A NAIVE QR
!*** DECOMPOSITION
!*** CHK0
!--------------------------
  SUBROUTINE MLS_NAIVEQR(N,M,AT,R)
    IMPLICIT REAL(8) (a-h,o-z)
    INTEGER::N,M
    REAL(8),DIMENSION(N)::TEMP
    REAL(8),DIMENSION(N,M)::E
    REAL(8),DIMENSION(N,M)::A
    REAL(8),DIMENSION(M,M)::R
    REAL(8),DIMENSION(N,M)::AT
!-------------------
!*** SETS R TO ZERO
!-------------------
    R=0.0D00    
    DO IM=1,M
       E(1:N,IM)=AT(1:N,IM)
       DO JM=1,IM-1
          CALL MLS_PROJQR(N,E(1:N,JM),AT(1:N,IM),TEMP(1:N))
          E(1:N,IM)=E(1:N,IM)-TEMP
       END DO
       CALL NRMALI(N,E(1:N,IM))
    END DO
    DO IM=1,M
       A(1:N,IM)=0.0D00
       DO JM=1,IM
          A(1:N,IM)=A(1:N,IM)+DOTPROD(N,E(1:N,JM),AT(1:N,IM))*E(1:N,JM)
       END DO
    END DO
    DO IM=1,M
       DO JM=IM,M
          R(IM,JM)=DOTPROD(N,E(1:N,IM),A(1:N,JM))
       END DO
    END DO
  END SUBROUTINE MLS_NAIVEQR

!--------------------
!*** FROM W, P AND B
!*** DETERMINES U2
!*** CHK0 (BOTH)
!--------------------  
  SUBROUTINE MLS_DETERMU2(M,N,W,P,U2)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER::M,N
    REAL(8),DIMENSION(N)::W
    REAL(8),DIMENSION(M,N)::B,U1,U2
    REAL(8),DIMENSION(M,N)::P
    REAL(8),DIMENSION(N,M)::AT
    REAL(8),DIMENSION(M,M)::R,ATEMP
    DO IM=1,M
       DO IN=1,N
          AT(IN,IM)=SQRT(W(IN))*P(IM,IN)
       END DO
    END DO
    DO IN=1,N
       DO IM=1,M
          B(IM,IN)=P(IM,IN)*W(IN)
       END DO
    END DO
!    GOTO 133
    CALL MLS_NAIVEQR(N,M,AT,R)
    CALL MLS_GENTRIANGSOLVE(.FALSE.,N,M,R,B,U1)
    CALL MLS_GENTRIANGSOLVE(.TRUE.,N,M,R,U1,U2)
    RETURN
133 CONTINUE
    DO IM=1,M
       DO JM=1,M
          ATEMP(IM,JM)=0.0D00
          DO IN=1,N
             ATEMP(IM,JM)=ATEMP(IM,JM)+P(IM,IN)*W(IN)*P(JM,IN)
          END DO
       END DO
    END DO
    CALL MATINV(M,DET,ATEMP,ATEMP,.FALSE.)
    CALL MATMAT(M,M,N,ATEMP,B,U2,0)
  END SUBROUTINE MLS_DETERMU2

!---------------------------------------
!*** DETERMINES GENERIC SHAPE FUNCTIONS
!*** AND DERIVATIVES
!*** CHK0
!---------------------------------------
  SUBROUTINE MLS_SFDER(NDI,N,M,POLYN,DPOLYN,DPOLYN2,U2,FF,DFF,DFF2)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER::M,N
    REAL(8),DIMENSION(M)::POLYN
    REAL(8),DIMENSION(NDI,M)::DPOLYN
    REAL(8),DIMENSION(NDI,NDI,M)::DPOLYN2
    REAL(8),DIMENSION(M,N)::U2
    REAL(8),DIMENSION(N)::FF
    REAL(8),DIMENSION(NDI,N)::DFF
    REAL(8),DIMENSION(NDI,NDI,N)::DFF2
    DO IN=1,N
       FF(IN)=DOTPROD(M,POLYN(1:M),U2(1:M,IN))
       DO ID=1,NDI
          DFF(ID,IN)=DOTPROD(M,DPOLYN(ID,1:M),U2(1:M,IN))
          DO JD=1,NDI
             DFF2(ID,JD,IN)=DOTPROD(M,DPOLYN2(ID,JD,1:M),U2(1:M,IN))
          END DO
       END DO
    END DO
  END SUBROUTINE MLS_SFDER

  SUBROUTINE MLS_POLYNQUAD1D(D,X,XBAR,POLYN,DPOLYN,DPOLYN2)
    IMPLICIT NONE
    DOUBLE PRECISION V(40),D,X(1),XBAR(1),POLYN(4),DPOLYN(1,4),DPOLYN2(1,1,4)
    V(34)=1D0/D
    V(33)=X(1)-XBAR(1)
    V(32)=1D0/D**2
    V(35)=2D0*V(32)
    V(31)=1D0/D**3
    V(23)=(V(33)*V(33))
    POLYN(1)=1D0
    POLYN(2)=V(33)*V(34)
    POLYN(3)=V(23)*V(32)
    POLYN(4)=V(31)*V(33)**3
    DPOLYN(1,1)=0D0
    DPOLYN(1,2)=V(34)
    DPOLYN(1,3)=V(33)*V(35)
    DPOLYN(1,4)=3D0*V(23)*V(31)
    DPOLYN2(1,1,1)=0D0
    DPOLYN2(1,1,2)=0D0
    DPOLYN2(1,1,3)=V(35)
    DPOLYN2(1,1,4)=6D0*V(31)*V(33)
  END SUBROUTINE MLS_POLYNQUAD1D

  SUBROUTINE MLS_POLYNQUAD2D(D,X,XBAR,POLYN,DPOLYN,DPOLYN2)
    IMPLICIT NONE
    DOUBLE PRECISION V(124),D,X(2),XBAR(2),POLYN(10),DPOLYN(2,10),DPOLYN2(2,2,10)
    V(115)=1D0/D
    V(113)=X(2)-XBAR(2)
    V(112)=X(1)-XBAR(1)
    V(111)=1D0/D**2
    V(116)=V(111)*V(113)
    V(110)=1D0/D**3
    V(119)=3D0*V(110)
    V(114)=V(110)*V(113)
    V(91)=2D0*V(112)
    V(87)=(V(112)*V(112))
    V(117)=V(110)*V(87)
    V(97)=2D0*V(113)
    V(89)=(V(113)*V(113))
    V(118)=V(110)*V(89)
    V(95)=V(114)*V(91)
    V(102)=2D0*V(111)
    V(104)=2D0*V(114)
    V(105)=V(110)*V(91)
    POLYN(1)=1D0
    POLYN(2)=V(112)*V(115)
    POLYN(3)=V(113)*V(115)
    POLYN(4)=V(111)*V(87)
    POLYN(5)=V(111)*V(89)
    POLYN(6)=V(112)*V(116)
    POLYN(7)=V(110)*V(112)**3
    POLYN(8)=V(110)*V(113)**3
    POLYN(9)=V(114)*V(87)
    POLYN(10)=V(112)*V(118)
    DPOLYN(1,1)=0D0
    DPOLYN(1,2)=V(115)
    DPOLYN(1,3)=0D0
    DPOLYN(1,4)=V(111)*V(91)
    DPOLYN(1,5)=0D0
    DPOLYN(1,6)=V(116)
    DPOLYN(1,7)=3D0*V(117)
    DPOLYN(1,8)=0D0
    DPOLYN(1,9)=V(95)
    DPOLYN(1,10)=V(118)
    DPOLYN(2,1)=0D0
    DPOLYN(2,2)=0D0
    DPOLYN(2,3)=V(115)
    DPOLYN(2,4)=0D0
    DPOLYN(2,5)=V(111)*V(97)
    DPOLYN(2,6)=V(111)*V(112)
    DPOLYN(2,7)=0D0
    DPOLYN(2,8)=3D0*V(118)
    DPOLYN(2,9)=V(117)
    DPOLYN(2,10)=V(95)
    DPOLYN2(1,1,1)=0D0
    DPOLYN2(1,1,2)=0D0
    DPOLYN2(1,1,3)=0D0
    DPOLYN2(1,1,4)=V(102)
    DPOLYN2(1,1,5)=0D0
    DPOLYN2(1,1,6)=0D0
    DPOLYN2(1,1,7)=V(119)*V(91)
    DPOLYN2(1,1,8)=0D0
    DPOLYN2(1,1,9)=V(104)
    DPOLYN2(1,1,10)=0D0
    DPOLYN2(1,2,1)=0D0
    DPOLYN2(1,2,2)=0D0
    DPOLYN2(1,2,3)=0D0
    DPOLYN2(1,2,4)=0D0
    DPOLYN2(1,2,5)=0D0
    DPOLYN2(1,2,6)=V(111)
    DPOLYN2(1,2,7)=0D0
    DPOLYN2(1,2,8)=0D0
    DPOLYN2(1,2,9)=V(105)
    DPOLYN2(1,2,10)=V(104)
    DPOLYN2(2,1,1)=0D0
    DPOLYN2(2,1,2)=0D0
    DPOLYN2(2,1,3)=0D0
    DPOLYN2(2,1,4)=0D0
    DPOLYN2(2,1,5)=0D0
    DPOLYN2(2,1,6)=V(111)
    DPOLYN2(2,1,7)=0D0
    DPOLYN2(2,1,8)=0D0
    DPOLYN2(2,1,9)=V(105)
    DPOLYN2(2,1,10)=V(104)
    DPOLYN2(2,2,1)=0D0
    DPOLYN2(2,2,2)=0D0
    DPOLYN2(2,2,3)=0D0
    DPOLYN2(2,2,4)=0D0
    DPOLYN2(2,2,5)=V(102)
    DPOLYN2(2,2,6)=0D0
    DPOLYN2(2,2,7)=0D0
    DPOLYN2(2,2,8)=V(119)*V(97)
    DPOLYN2(2,2,9)=0D0
    DPOLYN2(2,2,10)=V(105)
  END SUBROUTINE MLS_POLYNQUAD2D

  SUBROUTINE MLS_POLYNQUAD3D(D,X,XBAR,POLYN,DPOLYN,DPOLYN2)
    IMPLICIT NONE
    DOUBLE PRECISION v(340),d,x(3),xbar(3),polyn(20),dpolyn(3,20),dpolyn2(3,3,20)
    v(334)=1d0/d
    v(333)=x(3)-xbar(3)
    v(332)=x(2)-xbar(2)
    v(331)=x(1)-xbar(1)
    v(330)=1d0/d**2
    v(329)=1d0/d**3
    v(335)=3d0*v(329)
    v(295)=2d0*v(331)
    v(286)=(v(331)*v(331))
    v(316)=v(329)*v(332)
    v(305)=2d0*v(332)
    v(289)=(v(332)*v(332))
    v(313)=2d0*v(333)
    v(309)=v(329)*v(333)
    v(300)=v(309)*v(332)
    v(292)=(v(333)*v(333))
    v(297)=v(330)*v(332)
    v(298)=v(330)*v(333)
    v(301)=v(295)*v(316)
    v(302)=v(295)*v(309)
    v(303)=v(289)*v(329)
    v(304)=v(292)*v(329)
    v(307)=v(330)*v(331)
    v(311)=v(286)*v(329)
    v(312)=v(305)*v(309)
    v(318)=2d0*v(330)
    v(320)=v(305)*v(329)
    v(321)=2d0*v(309)
    v(322)=v(295)*v(329)
    v(324)=v(329)*v(331)
    polyn(1)=1d0
    polyn(2)=v(331)*v(334)
    polyn(3)=v(332)*v(334)
    polyn(4)=v(333)*v(334)
    polyn(5)=v(286)*v(330)
    polyn(6)=v(289)*v(330)
    polyn(7)=v(292)*v(330)
    polyn(8)=v(307)*v(332)
    polyn(9)=v(307)*v(333)
    polyn(10)=v(298)*v(332)
    polyn(11)=v(329)*v(331)**3
    polyn(12)=v(329)*v(332)**3
    polyn(13)=v(329)*v(333)**3
    polyn(14)=v(300)*v(331)
    polyn(15)=v(311)*v(332)
    polyn(16)=v(311)*v(333)
    polyn(17)=v(289)*v(324)
    polyn(18)=v(303)*v(333)
    polyn(19)=v(292)*v(324)
    polyn(20)=v(304)*v(332)
    dpolyn(1,1)=0d0
    dpolyn(1,2)=v(334)
    dpolyn(1,3)=0d0
    dpolyn(1,4)=0d0
    dpolyn(1,5)=v(295)*v(330)
    dpolyn(1,6)=0d0
    dpolyn(1,7)=0d0
    dpolyn(1,8)=v(297)
    dpolyn(1,9)=v(298)
    dpolyn(1,10)=0d0
    dpolyn(1,11)=3d0*v(311)
    dpolyn(1,12)=0d0
    dpolyn(1,13)=0d0
    dpolyn(1,14)=v(300)
    dpolyn(1,15)=v(301)
    dpolyn(1,16)=v(302)
    dpolyn(1,17)=v(303)
    dpolyn(1,18)=0d0
    dpolyn(1,19)=v(304)
    dpolyn(1,20)=0d0
    dpolyn(2,1)=0d0
    dpolyn(2,2)=0d0
    dpolyn(2,3)=v(334)
    dpolyn(2,4)=0d0
    dpolyn(2,5)=0d0
    dpolyn(2,6)=v(305)*v(330)
    dpolyn(2,7)=0d0
    dpolyn(2,8)=v(307)
    dpolyn(2,9)=0d0
    dpolyn(2,10)=v(298)
    dpolyn(2,11)=0d0
    dpolyn(2,12)=3d0*v(303)
    dpolyn(2,13)=0d0
    dpolyn(2,14)=v(309)*v(331)
    dpolyn(2,15)=v(311)
    dpolyn(2,16)=0d0
    dpolyn(2,17)=v(301)
    dpolyn(2,18)=v(312)
    dpolyn(2,19)=0d0
    dpolyn(2,20)=v(304)
    dpolyn(3,1)=0d0
    dpolyn(3,2)=0d0
    dpolyn(3,3)=0d0
    dpolyn(3,4)=v(334)
    dpolyn(3,5)=0d0
    dpolyn(3,6)=0d0
    dpolyn(3,7)=v(313)*v(330)
    dpolyn(3,8)=0d0
    dpolyn(3,9)=v(307)
    dpolyn(3,10)=v(297)
    dpolyn(3,11)=0d0
    dpolyn(3,12)=0d0
    dpolyn(3,13)=3d0*v(304)
    dpolyn(3,14)=v(316)*v(331)
    dpolyn(3,15)=0d0
    dpolyn(3,16)=v(311)
    dpolyn(3,17)=0d0
    dpolyn(3,18)=v(303)
    dpolyn(3,19)=v(302)
    dpolyn(3,20)=v(312)
    dpolyn2(1,1,1)=0d0
    dpolyn2(1,1,2)=0d0
    dpolyn2(1,1,3)=0d0
    dpolyn2(1,1,4)=0d0
    dpolyn2(1,1,5)=v(318)
    dpolyn2(1,1,6)=0d0
    dpolyn2(1,1,7)=0d0
    dpolyn2(1,1,8)=0d0
    dpolyn2(1,1,9)=0d0
    dpolyn2(1,1,10)=0d0
    dpolyn2(1,1,11)=v(295)*v(335)
    dpolyn2(1,1,12)=0d0
    dpolyn2(1,1,13)=0d0
    dpolyn2(1,1,14)=0d0
    dpolyn2(1,1,15)=v(320)
    dpolyn2(1,1,16)=v(321)
    dpolyn2(1,1,17)=0d0
    dpolyn2(1,1,18)=0d0
    dpolyn2(1,1,19)=0d0
    dpolyn2(1,1,20)=0d0
    dpolyn2(1,2,1)=0d0
    dpolyn2(1,2,2)=0d0
    dpolyn2(1,2,3)=0d0
    dpolyn2(1,2,4)=0d0
    dpolyn2(1,2,5)=0d0
    dpolyn2(1,2,6)=0d0
    dpolyn2(1,2,7)=0d0
    dpolyn2(1,2,8)=v(330)
    dpolyn2(1,2,9)=0d0
    dpolyn2(1,2,10)=0d0
    dpolyn2(1,2,11)=0d0
    dpolyn2(1,2,12)=0d0
    dpolyn2(1,2,13)=0d0
    dpolyn2(1,2,14)=v(309)
    dpolyn2(1,2,15)=v(322)
    dpolyn2(1,2,16)=0d0
    dpolyn2(1,2,17)=v(320)
    dpolyn2(1,2,18)=0d0
    dpolyn2(1,2,19)=0d0
    dpolyn2(1,2,20)=0d0
    dpolyn2(1,3,1)=0d0
    dpolyn2(1,3,2)=0d0
    dpolyn2(1,3,3)=0d0
    dpolyn2(1,3,4)=0d0
    dpolyn2(1,3,5)=0d0
    dpolyn2(1,3,6)=0d0
    dpolyn2(1,3,7)=0d0
    dpolyn2(1,3,8)=0d0
    dpolyn2(1,3,9)=v(330)
    dpolyn2(1,3,10)=0d0
    dpolyn2(1,3,11)=0d0
    dpolyn2(1,3,12)=0d0
    dpolyn2(1,3,13)=0d0
    dpolyn2(1,3,14)=v(316)
    dpolyn2(1,3,15)=0d0
    dpolyn2(1,3,16)=v(322)
    dpolyn2(1,3,17)=0d0
    dpolyn2(1,3,18)=0d0
    dpolyn2(1,3,19)=v(321)
    dpolyn2(1,3,20)=0d0
    dpolyn2(2,1,1)=0d0
    dpolyn2(2,1,2)=0d0
    dpolyn2(2,1,3)=0d0
    dpolyn2(2,1,4)=0d0
    dpolyn2(2,1,5)=0d0
    dpolyn2(2,1,6)=0d0
    dpolyn2(2,1,7)=0d0
    dpolyn2(2,1,8)=v(330)
    dpolyn2(2,1,9)=0d0
    dpolyn2(2,1,10)=0d0
    dpolyn2(2,1,11)=0d0
    dpolyn2(2,1,12)=0d0
    dpolyn2(2,1,13)=0d0
    dpolyn2(2,1,14)=v(309)
    dpolyn2(2,1,15)=v(322)
    dpolyn2(2,1,16)=0d0
    dpolyn2(2,1,17)=v(320)
    dpolyn2(2,1,18)=0d0
    dpolyn2(2,1,19)=0d0
    dpolyn2(2,1,20)=0d0
    dpolyn2(2,2,1)=0d0
    dpolyn2(2,2,2)=0d0
    dpolyn2(2,2,3)=0d0
    dpolyn2(2,2,4)=0d0
    dpolyn2(2,2,5)=0d0
    dpolyn2(2,2,6)=v(318)
    dpolyn2(2,2,7)=0d0
    dpolyn2(2,2,8)=0d0
    dpolyn2(2,2,9)=0d0
    dpolyn2(2,2,10)=0d0
    dpolyn2(2,2,11)=0d0
    dpolyn2(2,2,12)=3d0*v(320)
    dpolyn2(2,2,13)=0d0
    dpolyn2(2,2,14)=0d0
    dpolyn2(2,2,15)=0d0
    dpolyn2(2,2,16)=0d0
    dpolyn2(2,2,17)=v(322)
    dpolyn2(2,2,18)=v(321)
    dpolyn2(2,2,19)=0d0
    dpolyn2(2,2,20)=0d0
    dpolyn2(2,3,1)=0d0
    dpolyn2(2,3,2)=0d0
    dpolyn2(2,3,3)=0d0
    dpolyn2(2,3,4)=0d0
    dpolyn2(2,3,5)=0d0
    dpolyn2(2,3,6)=0d0
    dpolyn2(2,3,7)=0d0
    dpolyn2(2,3,8)=0d0
    dpolyn2(2,3,9)=0d0
    dpolyn2(2,3,10)=v(330)
    dpolyn2(2,3,11)=0d0
    dpolyn2(2,3,12)=0d0
    dpolyn2(2,3,13)=0d0
    dpolyn2(2,3,14)=v(324)
    dpolyn2(2,3,15)=0d0
    dpolyn2(2,3,16)=0d0
    dpolyn2(2,3,17)=0d0
    dpolyn2(2,3,18)=v(320)
    dpolyn2(2,3,19)=0d0
    dpolyn2(2,3,20)=v(321)
    dpolyn2(3,1,1)=0d0
    dpolyn2(3,1,2)=0d0
    dpolyn2(3,1,3)=0d0
    dpolyn2(3,1,4)=0d0
    dpolyn2(3,1,5)=0d0
    dpolyn2(3,1,6)=0d0
    dpolyn2(3,1,7)=0d0
    dpolyn2(3,1,8)=0d0
    dpolyn2(3,1,9)=v(330)
    dpolyn2(3,1,10)=0d0
    dpolyn2(3,1,11)=0d0
    dpolyn2(3,1,12)=0d0
    dpolyn2(3,1,13)=0d0
    dpolyn2(3,1,14)=v(316)
    dpolyn2(3,1,15)=0d0
    dpolyn2(3,1,16)=v(322)
    dpolyn2(3,1,17)=0d0
    dpolyn2(3,1,18)=0d0
    dpolyn2(3,1,19)=v(321)
    dpolyn2(3,1,20)=0d0
    dpolyn2(3,2,1)=0d0
    dpolyn2(3,2,2)=0d0
    dpolyn2(3,2,3)=0d0
    dpolyn2(3,2,4)=0d0
    dpolyn2(3,2,5)=0d0
    dpolyn2(3,2,6)=0d0
    dpolyn2(3,2,7)=0d0
    dpolyn2(3,2,8)=0d0
    dpolyn2(3,2,9)=0d0
    dpolyn2(3,2,10)=v(330)
    dpolyn2(3,2,11)=0d0
    dpolyn2(3,2,12)=0d0
    dpolyn2(3,2,13)=0d0
    dpolyn2(3,2,14)=v(324)
    dpolyn2(3,2,15)=0d0
    dpolyn2(3,2,16)=0d0
    dpolyn2(3,2,17)=0d0
    dpolyn2(3,2,18)=v(320)
    dpolyn2(3,2,19)=0d0
    dpolyn2(3,2,20)=v(321)
    dpolyn2(3,3,1)=0d0
    dpolyn2(3,3,2)=0d0
    dpolyn2(3,3,3)=0d0
    dpolyn2(3,3,4)=0d0
    dpolyn2(3,3,5)=0d0
    dpolyn2(3,3,6)=0d0
    dpolyn2(3,3,7)=v(318)
    dpolyn2(3,3,8)=0d0
    dpolyn2(3,3,9)=0d0
    dpolyn2(3,3,10)=0d0
    dpolyn2(3,3,11)=0d0
    dpolyn2(3,3,12)=0d0
    dpolyn2(3,3,13)=v(313)*v(335)
    dpolyn2(3,3,14)=0d0
    dpolyn2(3,3,15)=0d0
    dpolyn2(3,3,16)=0d0
    dpolyn2(3,3,17)=0d0
    dpolyn2(3,3,18)=0d0
    dpolyn2(3,3,19)=v(322)
    dpolyn2(3,3,20)=v(320)
  END SUBROUTINE MLS_POLYNQUAD3D

!---------------------------------
!*** MLS_POLYNQUADBASE DETERMINES
!*** POLY. BASE + DERIVATIVES
!*** 1D TO 3D
!*** CHK0
!---------------------------------  
  SUBROUTINE MLS_POLYNCUBBASE(NDI,DTEMP,X,XBAR,POLYN,DPOLYN,DPOLYN2,M)
    IMPLICIT REAL(8) (a-h,o-z)
!*** CHECK THIS PARAMETER FOR UPDATES (MPOL)
    REAL(8)::D,DTEMP
    INTEGER,PARAMETER::MPOL=20
    REAL(8),DIMENSION(NDI)::X,XBAR
    REAL(8),DIMENSION(MPOL)::POLYN
    REAL(8),DIMENSION(NDI,MPOL)::DPOLYN
    REAL(8),DIMENSION(NDI,NDI,MPOL)::DPOLYN2
    D=DTEMP
    SELECT CASE(NDI)
    CASE(1)
       M=4
       CALL MLS_POLYNQUAD1D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M),DPOLYN2(1:NDI,1:NDI,1:M))
    CASE(2)
       M=10
       CALL MLS_POLYNQUAD2D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M),DPOLYN2(1:NDI,1:NDI,1:M))
    CASE(3)      
       M=20      
       CALL MLS_POLYNQUAD3D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M),DPOLYN2(1:NDI,1:NDI,1:M))
    CASE DEFAULT
       STOP "CASE NOT FOUND IN MLS_POLYNCUBBASE"
    END SELECT
  END SUBROUTINE MLS_POLYNCUBBASE

!---------------------------------
!*** MLS_POLYNQUADBASE DETERMINES
!*** POLY. BASE + DERIVATIVES
!*** 1D TO 3D
!*** CHK0
!---------------------------------  
  SUBROUTINE MLS_POLYNQBASE(NDI,DTEMP,X,XBAR,POLYN,DPOLYN,DPOLYN2,M)
    IMPLICIT REAL(8) (a-h,o-z)
!*** CHECK THIS PARAMETER FOR UPDATES (MPOL)
    REAL(8)::D,DTEMP
    INTEGER,PARAMETER::MPOL=20
    REAL(8),DIMENSION(NDI)::X,XBAR
    REAL(8),DIMENSION(MPOL)::POLYN
    REAL(8),DIMENSION(NDI,MPOL)::DPOLYN
    REAL(8),DIMENSION(NDI,NDI,MPOL)::DPOLYN2
    D=DTEMP
    SELECT CASE(NDI)
    CASE(1)
       M=3
       CALL MLS_POLYNQUAD1D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M),DPOLYN2(1:NDI,1:NDI,1:M))
    CASE(2)
       M=6
       CALL MLS_POLYNQUAD2D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M),DPOLYN2(1:NDI,1:NDI,1:M))
    CASE(3)      
       M=10       
       CALL MLS_POLYNQUAD3D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M),DPOLYN2(1:NDI,1:NDI,1:M))
    CASE DEFAULT
       STOP "CASE NOT FOUND IN MLS_POLYNQBASE"
    END SELECT
  END SUBROUTINE MLS_POLYNQBASE

!---------------------------------
!*** MLS_POLYNQUADBASE DETERMINES
!*** POLY. BASE + DERIVATIVES
!*** 1D TO 3D
!*** CHK0
!---------------------------------  
  SUBROUTINE MLS_POLYNLINBASE(NDI,DTEMP,X,XBAR,POLYN,DPOLYN,DPOLYN2,M)
    IMPLICIT REAL(8) (a-h,o-z)
!*** CHECK THIS PARAMETER FOR UPDATES (MPOL)
    REAL(8)::D,DTEMP
    INTEGER,PARAMETER::MPOL=20
    REAL(8),DIMENSION(NDI)::X,XBAR
    REAL(8),DIMENSION(MPOL)::POLYN
    REAL(8),DIMENSION(NDI,MPOL)::DPOLYN
    REAL(8),DIMENSION(NDI,NDI,MPOL)::DPOLYN2
    D=DTEMP
    SELECT CASE(NDI)
    CASE(1)
       M=2
       CALL MLS_POLYNQUAD1D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M),DPOLYN2(1:NDI,1:NDI,1:M))
    CASE(2)
       M=3
       CALL MLS_POLYNQUAD2D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M),DPOLYN2(1:NDI,1:NDI,1:M))
    CASE(3)      
       M=4       
       CALL MLS_POLYNQUAD3D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M),DPOLYN2(1:NDI,1:NDI,1:M))
    CASE DEFAULT
       STOP "CASE NOT FOUND IN MLS_POLYNLINBASE"
    END SELECT
  END SUBROUTINE MLS_POLYNLINBASE

!----------------------------
!*** WEIGHT FUNCTION FOR MLS
!*** NDI
!*** CHK0
!----------------------------
  SUBROUTINE MLS_WEIGHT(NDI,D,TOL,XI,X,W)
    IMPLICIT REAL(8) (A-H,O-Z)
    REAL(8),DIMENSION(NDI)::XI,X
    REAL(8),PARAMETER::PI=4.0D00*ATAN(1.0D00)
    REAL(8),PARAMETER::BOTTOM=1.0D-16
    S=RNORM2(NDI,XI-X)
!!    WMAX=1.0D00/(TOL**2.0D00+TOL**4.0D00)
!!    WBOT=BOTTOM*WMAX
!!    W=MAX((1.0D00/((S*S)/(D*D)+TOL*TOL))-(1.0D00/(1.0D00+TOL*TOL)),WBOT)
    W=(1.0D00/((S*S)/(D*D)+TOL*TOL))-(1.0D00/(1.0D00+TOL*TOL))
!!    W=MAX(EXP(-S/(1.0D-6*D)),1.0d-15)
!!    W=SQRT(1.0d00/(TOL*PI))*EXP(-((S/D)**2)/TOL)
!!    W=1.0d00/((S/D)+TOL)-1.0d00/(1.0d00+TOL)
!!    W=MAX(W,0.0D00)
    W=MAX(W,0.0d00)
!    W=1.0d00
!    w=1.0d00
!!    W=(1.0D00/((S*S)/(D*D)+TOL*TOL))-(1.0D00/(1.0D00+TOL*TOL))
!!    W=MAX(W,1.0d-20)
!    w=1.0d00
  END SUBROUTINE MLS_WEIGHT

!-------------------------------------------------
!*** DETERMINES U2 EVALUATED IN A GIVEN COORDINATE
!    CALL U2ATACOORDINATE(D,TOL,NDI,N,XN,X,U2)
!-------------------------------------------------
  SUBROUTINE MLS_U2ATACOORDINATE(IMESHLESS,D,TOL,NDI,N,XBAR,XN,X,U2)
    IMPLICIT REAL(8) (A-H,O-Z)
    REAL(8)::D,TOL
    INTEGER,PARAMETER::MPOL=20
    REAL(8),DIMENSION(N)::W
    REAL(8),DIMENSION(MPOL,N)::P
    REAL(8),DIMENSION(MPOL,N)::U2
    REAL(8),DIMENSION(NDI,N)::XN
    REAL(8),DIMENSION(NDI)::X,XBAR
    REAL(8),DIMENSION(3,MPOL)::TRASH
    REAL(8),DIMENSION(3,3,MPOL)::TRASH2
!-------------------------------
!*** NOW WE MUST DEFINE MATRIX P
!*** DERIVATIVES ARE DISCARDED
!-------------------------------
    SELECT CASE(IMESHLESS)
    CASE(-1)
       DO IN=1,N
          CALL MLS_POLYNLINBASE(NDI,D,XN(1:NDI,IN),XBAR(1:NDI),P(1:MPOL,IN),TRASH(1:NDI,1:MPOL),TRASH2(1:NDI,1:NDI,1:MPOL),M)
       END DO
    CASE(-2)
       DO IN=1,N
          CALL MLS_POLYNQBASE(NDI,D,XN(1:NDI,IN),XBAR(1:NDI),P(1:MPOL,IN),TRASH(1:NDI,1:MPOL),TRASH2(1:NDI,1:NDI,1:MPOL),M)
       END DO
    CASE(-3)
       DO IN=1,N
          CALL MLS_POLYNCUBBASE(NDI,D,XN(1:NDI,IN),XBAR(1:NDI),P(1:MPOL,IN),TRASH(1:NDI,1:MPOL),TRASH2(1:NDI,1:NDI,1:MPOL),M)
       END DO
    CASE DEFAULT
       STOP "WRONG REQUEST TO U2ATACOORDINATE"
    END SELECT
!------------------------
!*** NOW WE MUST DEFINE W
!------------------------
    DO IN=1,N
       CALL MLS_WEIGHT(NDI,D,TOL,XN(1:NDI,IN),X(1:NDI),W(IN))
    END DO
!------------------------------
!*** FINALIZE IT WITH SPECIFICS
!------------------------------
    CALL MLS_DETERMU2(M,N,W(1:N),P(1:M,1:N),U2(1:M,1:N))
  END SUBROUTINE MLS_U2ATACOORDINATE

!-------------------------
!*** GETS INFLUENCE NODES
!*** FOR A GIVEN POINT
!*** chk0
!-------------------------
  SUBROUTINE MLS_GETSNODES(RMAX,X,NDI,NTOT,XTOT,N,LISTN)
    IMPLICIT REAL(8)(A-H,O-Z)
    REAL(8)::RMAX
    REAL(8),DIMENSION(NDI,NTOT)::XTOT
    REAL(8),DIMENSION(NDI)::X
    INTEGER::NDI
    INTEGER,DIMENSION(:),ALLOCATABLE::LISTN
    N=0
    DO ITOT=1,NTOT
       IF(RNORM2(NDI,XTOT(1:NDI,ITOT)-X(1:NDI)).LE.RMAX)THEN
          N=N+1
       END IF
    END DO
    ALLOCATE(LISTN(N))
!--------------------------------------
!*** NOW INSERTS THE NODES IN THE LIST
!--------------------------------------
    N=0
    DO ITOT=1,NTOT
       IF(RNORM2(NDI,XTOT(1:NDI,ITOT)-X(1:NDI)).LE.RMAX)THEN
          N=N+1
          LISTN(N)=ITOT
       END IF
    END DO
  END SUBROUTINE MLS_GETSNODES

!---------------------------
!*** GETS PRECISELY N NODES
!---------------------------
  SUBROUTINE MLS_GETSNNODES(RMAX,X,NDI,NTOT,XTOT,N,LISTN)
    IMPLICIT REAL(8)(A-H,O-Z)
    REAL(8)::RMAX
    REAL(8),DIMENSION(NDI,NTOT)::XTOT
    REAL(8),DIMENSION(NDI)::X
    REAL(8),DIMENSION(NTOT)::RANKING
    INTEGER,DIMENSION(NTOT)::PER
    INTEGER::NDI
    INTEGER,DIMENSION(:),ALLOCATABLE::LISTN
    IF(N.LE.0)THEN
       CALL MLS_GETSNODES(RMAX,X,NDI,NTOT,XTOT,N,LISTN)
    ELSE
       DO INO=1,NTOT
          RANKING(INO)=RNORM2(NDI,X(1:NDI)-XTOT(1:NDI,INO))
       END DO
!--------------------
!*** SORT BY RANKING
!--------------------
       CALL SORT(NTOT,RANKING,PER)
       ALLOCATE(LISTN(N))
       DO I=1,N
          LISTN(I)=PER(I)
       END DO
    END IF
  END SUBROUTINE MLS_GETSNNODES

  SUBROUTINE MLS_CAUCHYGREEN(NDI,N,XTOT,XDEF,DFF,F,CAUCHYGREEN)
    IMPLICIT REAL(8) (a-h,o-z)
    REAL(8),DIMENSION(NDI,N)::XTOT,XDEF
    REAL(8),DIMENSION(NDI,NDI)::C,F,FO,CO
    REAL(8),DIMENSION(NDI,N)::DFF
    REAL(8),DIMENSION(NDI*(NDI+1)/2)::CAUCHYGREEN
!-------------------------
!*** DEFORMATION GRADIENT
!*** BEFORE CORRECTION
!-------------------------
    F=0.0D00
    DO IN=1,N
       DO JD=1,NDI
          DO ID=1,NDI
             F(ID,JD)=F(ID,JD)+DFF(JD,IN)*XDEF(ID,IN)
          END DO
       END DO
    END DO
!------------------------------
!*** RIGHT CAUCHY-GREEN TENSOR
!------------------------------
    CALL MATMAT(NDI,NDI,NDI,F,F,C,2)
    DO ID=1,NDI*(NDI+1)/2
       CALL APOMAT(I1,I2,ID,NDI)
       CAUCHYGREEN(ID)=C(I1,I2)
    END DO
  END SUBROUTINE MLS_CAUCHYGREEN

  SUBROUTINE MLS_STRAIN(NDI,N,XTOT,XDEF,DFF,F,STRAIN)
    IMPLICIT REAL(8) (a-h,o-z)
    REAL(8),DIMENSION(NDI,N)::XTOT,XDEF
    REAL(8),DIMENSION(NDI,NDI)::C,F,FO,CO
    REAL(8),DIMENSION(NDI,N)::DFF
    REAL(8),DIMENSION(NDI*(NDI+1)/2)::STRAIN,STRAINO
!-------------------------
!*** DEFORMATION GRADIENT
!*** BEFORE CORRECTION
!-------------------------
    F=0.0D00
    DO IN=1,N
       DO JD=1,NDI
          DO ID=1,NDI
             F(ID,JD)=F(ID,JD)+DFF(JD,IN)*XDEF(ID,IN)
          END DO
       END DO
    END DO
!------------------------------
!*** RIGHT CAUCHY-GREEN TENSOR
!------------------------------
    CALL MATMAT(NDI,NDI,NDI,F,F,C,2)
    DO ID=1,NDI*(NDI+1)/2
       CALL APOMAT(I1,I2,ID,NDI)
       STRAIN(ID)=0.5D00*(C(I1,I2)-DELTAK(I1,I2))
    END DO
    DO ID=NDI+1,NDI*(NDI+1)/2
       STRAIN(ID)=2.0D00*STRAIN(ID)
    END DO
  END SUBROUTINE MLS_STRAIN

!----------------------------------------------------
!*** DETERMINES THE DEFORMATION GRADIENT
!*** AND GREEN-LAGRANGE STRAIN IN ENGINEERING FORM
!*** MLS STUFF
!*** CHK0,1
!
! REAL(8),DIMENSION(NDI,NDNOS)::DER
! REAL(8),DIMENSION(NDI,NDI,NDNOS)::DER2
!----------------------------------------------------
  SUBROUTINE MLS_STRAIN2(NDI,N,GRAD,GRAD2,E,E2)
    IMPLICIT REAL(8) (a-h,o-z)
!----------------------------------------------------
    REAL(8),DIMENSION(NDI,NDI)::GRAD
    REAL(8),DIMENSION(NDI,NDI,NDI)::GRAD2
!----------------------------------------------------
    REAL(8),DIMENSION(NDI*(NDI+1)/2)::E
    REAL(8),DIMENSION(NDI*(NDI+1)/2,NDI)::E2
!----------------------------------------------------
    NVOIGT=NDI*(NDI+1)/2
    DO IJ=1,NVOIGT
       CALL APOMAT(I,J,IJ,NDI)
!E
       E(IJ)=0.0D00
       DO K=1,NDI
          E(IJ)=E(IJ)+0.5D00*GRAD(K,I)*GRAD(K,J)
       END DO
       E(IJ)=E(IJ)-0.5D00*DELTAK(I,J)
!E2
       DO L=1,NDI
          E2(IJ,L)=0.0D00
          DO K=1,NDI
             E2(IJ,L)=E2(IJ,L)+0.5D00*GRAD2(K,I,L)*GRAD(K,J)+0.5D00*GRAD(K,I)*GRAD2(K,J,L)
          END DO
       END DO
    END DO
  END SUBROUTINE MLS_STRAIN2

  SUBROUTINE MLS_2DFORCE(NK,F,S,NUCLEUS)
    IMPLICIT NONE
    DOUBLE PRECISION V(27),NK(2),F(2,2),S(3),NUCLEUS(2)
    V(22)=NK(2)*S(2)
    V(21)=NK(1)*S(1)
    NUCLEUS(1)=(F(1,2)*NK(1)+F(1,1)*NK(2))*S(3)+F(1,1)*V(21)+F(1,2)*V(22)
    NUCLEUS(2)=(F(2,2)*NK(1)+F(2,1)*NK(2))*S(3)+F(2,1)*V(21)+F(2,2)*V(22)
  END SUBROUTINE MLS_2DFORCE

  SUBROUTINE MLS_2D(NK,NL,F,S,DS,KERNEL)
    IMPLICIT NONE
    DOUBLE PRECISION V(67),NK(2),NL(2),F(2,2),S(3),DS(3,3),KERNEL(2,2)
    V(53)=NK(2)*(NL(2)*S(2)+NL(1)*S(3))+NK(1)*(NL(1)*S(1)+NL(2)*S(3))
    V(52)=F(2,2)*NK(1)+F(2,1)*NK(2)
    V(51)=F(1,2)*NK(1)+F(1,1)*NK(2)
    V(50)=F(2,2)*NK(2)
    V(49)=F(1,2)*NK(2)
    V(48)=F(2,1)*NK(1)
    V(47)=F(1,1)*NK(1)
    V(46)=F(2,2)*NL(1)+F(2,1)*NL(2)
    V(45)=F(1,2)*NL(1)+F(1,1)*NL(2)
    V(44)=F(2,2)*NL(2)
    V(60)=DS(2,2)*V(44)
    V(59)=DS(1,2)*V(44)+DS(1,3)*V(46)
    V(43)=F(1,2)*NL(2)
    V(54)=DS(1,2)*V(43)+DS(1,3)*V(45)
    V(42)=F(2,1)*NL(1)
    V(62)=DS(3,1)*V(42)+DS(3,2)*V(44)+DS(3,3)*V(46)
    V(61)=DS(2,1)*V(42)+DS(2,3)*V(46)
    V(58)=DS(1,1)*V(42)
    V(41)=F(1,1)*NL(1)
    V(56)=DS(3,1)*V(41)+DS(3,2)*V(43)+DS(3,3)*V(45)
    V(55)=DS(2,1)*V(41)+DS(2,3)*V(45)
    V(38)=V(47)*V(58)
    V(39)=V(49)*V(60)
    V(57)=V(38)+V(39)
    KERNEL(1,1)=V(53)+V(47)*(DS(1,1)*V(41)+V(54))+V(49)*(DS(2,2)*V(43)+V(55))+V(51)*V(56)
    KERNEL(1,2)=V(57)+V(47)*V(59)+V(49)*V(61)+V(51)*V(62)
    KERNEL(2,1)=V(48)*V(54)+V(50)*V(55)+V(52)*V(56)+V(57)
    KERNEL(2,2)=V(53)+V(48)*(V(58)+V(59))+V(50)*(V(60)+V(61))+V(52)*V(62)
  END SUBROUTINE MLS_2D

  SUBROUTINE MLS_3DFORCE(NK,F,S,NUCLEUS)
    IMPLICIT NONE
    DOUBLE PRECISION V(51),NK(3),F(3,3),S(6),NUCLEUS(3)
    V(46)=NK(3)*S(3)
    V(45)=NK(2)*S(2)
    V(44)=NK(1)*S(1)
    NUCLEUS(1)=(F(1,2)*NK(1)+F(1,1)*NK(2))*S(4)+(F(1,3)*NK(1)+F(1,1)*NK(3))*S(5)+(F(1,3)*NK(2)+F(1,2)*NK(3))*S(6)+F(1,1)*V&
         &(44)+F(1,2)*V(45)+F(1,3)*V(46)
    NUCLEUS(2)=(F(2,2)*NK(1)+F(2,1)*NK(2))*S(4)+(F(2,3)*NK(1)+F(2,1)*NK(3))*S(5)+(F(2,3)*NK(2)+F(2,2)*NK(3))*S(6)+F(2,1)*V&
         &(44)+F(2,2)*V(45)+F(2,3)*V(46)
    NUCLEUS(3)=(F(3,2)*NK(1)+F(3,1)*NK(2))*S(4)+(F(3,3)*NK(1)+F(3,1)*NK(3))*S(5)+(F(3,3)*NK(2)+F(3,2)*NK(3))*S(6)+F(3,1)*V&
         &(44)+F(3,2)*V(45)+F(3,3)*V(46)
  END SUBROUTINE MLS_3DFORCE

  SUBROUTINE MLS_3D(NK,NL,F,S,DS,KERNEL)
    IMPLICIT NONE
    DOUBLE PRECISION V(173),NK(3),NL(3),F(3,3),S(6),DS(6,6),KERNEL(3,3)
    V(144)=NK(1)*(NL(1)*S(1)+NL(2)*S(4)+NL(3)*S(5))+NK(3)*(NL(3)*S(3)+NL(1)*S(5)+NL(2)*S(6))+NK(2)*(NL(2)*S(2)+NL(1)*S(4)&
         &+NL(3)*S(6))
    V(143)=F(3,3)*NK(2)+F(3,2)*NK(3)
    V(142)=F(2,3)*NK(2)+F(2,2)*NK(3)
    V(141)=F(1,3)*NK(2)+F(1,2)*NK(3)
    V(140)=F(3,3)*NK(1)+F(3,1)*NK(3)
    V(139)=F(2,3)*NK(1)+F(2,1)*NK(3)
    V(138)=F(1,3)*NK(1)+F(1,1)*NK(3)
    V(137)=F(3,2)*NK(1)+F(3,1)*NK(2)
    V(136)=F(2,2)*NK(1)+F(2,1)*NK(2)
    V(135)=F(1,2)*NK(1)+F(1,1)*NK(2)
    V(134)=F(3,3)*NK(3)
    V(133)=F(2,3)*NK(3)
    V(132)=F(1,3)*NK(3)
    V(131)=F(3,2)*NK(2)
    V(130)=F(2,2)*NK(2)
    V(129)=F(1,2)*NK(2)
    V(128)=F(3,1)*NK(1)
    V(127)=F(2,1)*NK(1)
    V(126)=F(1,1)*NK(1)
    V(125)=F(3,3)*NL(2)+F(3,2)*NL(3)
    V(124)=F(2,3)*NL(2)+F(2,2)*NL(3)
    V(123)=F(1,3)*NL(2)+F(1,2)*NL(3)
    V(122)=F(3,3)*NL(1)+F(3,1)*NL(3)
    V(121)=F(2,3)*NL(1)+F(2,1)*NL(3)
    V(120)=F(1,3)*NL(1)+F(1,1)*NL(3)
    V(119)=F(3,2)*NL(1)+F(3,1)*NL(2)
    V(118)=F(2,2)*NL(1)+F(2,1)*NL(2)
    V(117)=F(1,2)*NL(1)+F(1,1)*NL(2)
    V(116)=F(3,3)*NL(3)
    V(115)=F(2,3)*NL(3)
    V(168)=DS(3,3)*V(115)
    V(114)=F(1,3)*NL(3)
    V(147)=DS(3,3)*V(114)
    V(113)=F(3,2)*NL(2)
    V(160)=DS(1,2)*V(113)+DS(1,3)*V(116)+DS(1,4)*V(119)+DS(1,5)*V(122)+DS(1,6)*V(125)
    V(112)=F(2,2)*NL(2)
    V(167)=DS(2,2)*V(112)
    V(154)=DS(1,2)*V(112)+DS(1,3)*V(115)+DS(1,4)*V(118)+DS(1,5)*V(121)+DS(1,6)*V(124)
    V(111)=F(1,2)*NL(2)
    V(148)=DS(1,2)*V(111)+DS(1,3)*V(114)+DS(1,4)*V(117)+DS(1,5)*V(120)+DS(1,6)*V(123)
    V(146)=DS(2,2)*V(111)
    V(110)=F(3,1)*NL(1)
    V(165)=DS(6,1)*V(110)+DS(6,2)*V(113)+DS(6,3)*V(116)+DS(6,4)*V(119)+DS(6,5)*V(122)+DS(6,6)*V(125)
    V(164)=DS(5,1)*V(110)+DS(5,2)*V(113)+DS(5,3)*V(116)+DS(5,4)*V(119)+DS(5,5)*V(122)+DS(5,6)*V(125)
    V(163)=DS(4,1)*V(110)+DS(4,2)*V(113)+DS(4,3)*V(116)+DS(4,4)*V(119)+DS(4,5)*V(122)+DS(4,6)*V(125)
    V(162)=DS(3,1)*V(110)+DS(3,2)*V(113)+DS(3,4)*V(119)+DS(3,5)*V(122)+DS(3,6)*V(125)
    V(161)=DS(2,1)*V(110)+DS(2,3)*V(116)+DS(2,4)*V(119)+DS(2,5)*V(122)+DS(2,6)*V(125)
    V(109)=F(2,1)*NL(1)
    V(166)=DS(1,1)*V(109)
    V(159)=DS(6,1)*V(109)+DS(6,2)*V(112)+DS(6,3)*V(115)+DS(6,4)*V(118)+DS(6,5)*V(121)+DS(6,6)*V(124)
    V(158)=DS(5,1)*V(109)+DS(5,2)*V(112)+DS(5,3)*V(115)+DS(5,4)*V(118)+DS(5,5)*V(121)+DS(5,6)*V(124)
    V(157)=DS(4,1)*V(109)+DS(4,2)*V(112)+DS(4,3)*V(115)+DS(4,4)*V(118)+DS(4,5)*V(121)+DS(4,6)*V(124)
    V(156)=DS(3,1)*V(109)+DS(3,2)*V(112)+DS(3,4)*V(118)+DS(3,5)*V(121)+DS(3,6)*V(124)
    V(155)=DS(2,1)*V(109)+DS(2,3)*V(115)+DS(2,4)*V(118)+DS(2,5)*V(121)+DS(2,6)*V(124)
    V(108)=F(1,1)*NL(1)
    V(153)=DS(6,1)*V(108)+DS(6,2)*V(111)+DS(6,3)*V(114)+DS(6,4)*V(117)+DS(6,5)*V(120)+DS(6,6)*V(123)
    V(152)=DS(5,1)*V(108)+DS(5,2)*V(111)+DS(5,3)*V(114)+DS(5,4)*V(117)+DS(5,5)*V(120)+DS(5,6)*V(123)
    V(151)=DS(4,1)*V(108)+DS(4,2)*V(111)+DS(4,3)*V(114)+DS(4,4)*V(117)+DS(4,5)*V(120)+DS(4,6)*V(123)
    V(150)=DS(3,1)*V(108)+DS(3,2)*V(111)+DS(3,4)*V(117)+DS(3,5)*V(120)+DS(3,6)*V(123)
    V(149)=DS(2,1)*V(108)+DS(2,3)*V(114)+DS(2,4)*V(117)+DS(2,5)*V(120)+DS(2,6)*V(123)
    V(145)=DS(1,1)*V(108)
    V(104)=V(127)*V(145)+V(130)*V(146)+V(133)*V(147)
    V(105)=V(128)*V(145)+V(131)*V(146)+V(134)*V(147)
    V(106)=V(128)*V(166)+V(131)*V(167)+V(134)*V(168)
    KERNEL(1,1)=V(144)+V(126)*(V(145)+V(148))+V(129)*(V(146)+V(149))+V(132)*(V(147)+V(150))+V(135)*V(151)+V(138)*V(152)+V&
         &(141)*V(153)
    KERNEL(1,2)=V(104)+V(126)*V(154)+V(129)*V(155)+V(132)*V(156)+V(135)*V(157)+V(138)*V(158)+V(141)*V(159)
    KERNEL(1,3)=V(105)+V(126)*V(160)+V(129)*V(161)+V(132)*V(162)+V(135)*V(163)+V(138)*V(164)+V(141)*V(165)
    KERNEL(2,1)=V(104)+V(127)*V(148)+V(130)*V(149)+V(133)*V(150)+V(136)*V(151)+V(139)*V(152)+V(142)*V(153)
    KERNEL(2,2)=V(144)+V(136)*V(157)+V(139)*V(158)+V(142)*V(159)+V(127)*(V(154)+V(166))+V(130)*(V(155)+V(167))+V(133)*(V&
         &(156)+V(168))
    KERNEL(2,3)=V(106)+V(127)*V(160)+V(130)*V(161)+V(133)*V(162)+V(136)*V(163)+V(139)*V(164)+V(142)*V(165)
    KERNEL(3,1)=V(105)+V(128)*V(148)+V(131)*V(149)+V(134)*V(150)+V(137)*V(151)+V(140)*V(152)+V(143)*V(153)
    KERNEL(3,2)=V(106)+V(128)*V(154)+V(131)*V(155)+V(134)*V(156)+V(137)*V(157)+V(140)*V(158)+V(143)*V(159)
    KERNEL(3,3)=V(144)+V(128)*(DS(1,1)*V(110)+V(160))+V(131)*(DS(2,2)*V(113)+V(161))+V(134)*(DS(3,3)*V(116)+V(162))+V(137&
         &)*V(163)+V(140)*V(164)+V(143)*V(165)
  END SUBROUTINE MLS_3D

  SUBROUTINE mls_fbarforce(dck,s,nucleus)
    IMPLICIT NONE
    DOUBLE PRECISION v(54),dck(6,3),s(6),nucleus(3)
    v(49)=s(3)/2d0
    v(48)=s(2)/2d0
    v(47)=s(1)/2d0
    nucleus(1)=dck(4,1)*s(4)+dck(5,1)*s(5)+dck(6,1)*s(6)+dck(1,1)*v(47)+dck(2,1)*v(48)+dck(3,1)*v(49)
    nucleus(2)=dck(4,2)*s(4)+dck(5,2)*s(5)+dck(6,2)*s(6)+dck(1,2)*v(47)+dck(2,2)*v(48)+dck(3,2)*v(49)
    nucleus(3)=dck(4,3)*s(4)+dck(5,3)*s(5)+dck(6,3)*s(6)+dck(1,3)*v(47)+dck(2,3)*v(48)+dck(3,3)*v(49)
  END SUBROUTINE mls_fbarforce

  SUBROUTINE mls_fbarforce_2d(dck,s,nucleus)
    IMPLICIT NONE
    DOUBLE PRECISION v(28),dck(3,2),s(3),nucleus(2)
    v(23)=s(2)/2d0
    v(22)=s(1)/2d0
    nucleus(1)=dck(3,1)*s(3)+dck(1,1)*v(22)+dck(2,1)*v(23)
    nucleus(2)=dck(3,2)*s(3)+dck(1,2)*v(22)+dck(2,2)*v(23)
    nucleus(3)=nucleus(3)
  END SUBROUTINE mls_fbarforce_2d
  
  SUBROUTINE mls_fbar(dck,dcl,d2c,s,ds,kernel)
    IMPLICIT NONE
    DOUBLE PRECISION v(303),dck(6,3),dcl(6,3),d2c(6,3,3),s(6),ds(6,6),kernel(3,3)
    v(280)=s(3)/2d0
    v(279)=s(2)/2d0
    v(278)=s(1)/2d0
    v(277)=dcl(3,3)/2d0
    v(276)=dcl(2,3)/2d0
    v(275)=dcl(1,3)/2d0
    v(274)=ds(4,3)*v(277)
    v(273)=ds(4,2)*v(276)
    v(272)=ds(4,1)*v(275)
    v(298)=dcl(4,3)*ds(4,4)+dcl(5,3)*ds(4,5)+dcl(6,3)*ds(4,6)+v(272)+v(273)+v(274)
    v(271)=dcl(6,3)/2d0
    v(270)=dcl(5,3)/2d0
    v(269)=dcl(4,3)/2d0
    v(268)=dcl(3,3)/4d0
    v(267)=dcl(2,3)/4d0
    v(266)=dcl(1,3)/4d0
    v(265)=ds(1,6)*v(271)
    v(264)=ds(1,5)*v(270)
    v(263)=ds(1,4)*v(269)
    v(262)=ds(1,3)*v(268)
    v(261)=ds(1,2)*v(267)
    v(260)=ds(1,1)*v(266)
    v(297)=v(260)+v(261)+v(262)+v(263)+v(264)+v(265)
    v(259)=dcl(3,2)/2d0
    v(258)=dcl(2,2)/2d0
    v(257)=dcl(1,2)/2d0
    v(256)=ds(4,3)*v(259)
    v(255)=ds(4,2)*v(258)
    v(254)=ds(4,1)*v(257)
    v(292)=dcl(4,2)*ds(4,4)+dcl(5,2)*ds(4,5)+dcl(6,2)*ds(4,6)+v(254)+v(255)+v(256)
    v(253)=dcl(6,2)/2d0
    v(252)=dcl(5,2)/2d0
    v(251)=dcl(4,2)/2d0
    v(250)=dcl(3,2)/4d0
    v(249)=dcl(2,2)/4d0
    v(248)=dcl(1,2)/4d0
    v(247)=ds(1,6)*v(253)
    v(246)=ds(1,5)*v(252)
    v(245)=ds(1,4)*v(251)
    v(244)=ds(1,3)*v(250)
    v(243)=ds(1,2)*v(249)
    v(242)=ds(1,1)*v(248)
    v(291)=v(242)+v(243)+v(244)+v(245)+v(246)+v(247)
    v(241)=dcl(3,1)/2d0
    v(240)=dcl(2,1)/2d0
    v(239)=dcl(1,1)/2d0
    v(238)=ds(4,3)*v(241)
    v(237)=ds(4,2)*v(240)
    v(236)=ds(4,1)*v(239)
    v(286)=dcl(4,1)*ds(4,4)+dcl(5,1)*ds(4,5)+dcl(6,1)*ds(4,6)+v(236)+v(237)+v(238)
    v(235)=dcl(6,1)/2d0
    v(234)=dcl(5,1)/2d0
    v(233)=dcl(4,1)/2d0
    v(232)=dcl(3,1)/4d0
    v(231)=dcl(2,1)/4d0
    v(230)=dcl(1,1)/4d0
    v(229)=ds(1,6)*v(235)
    v(228)=ds(1,5)*v(234)
    v(227)=ds(1,4)*v(233)
    v(226)=ds(1,3)*v(232)
    v(225)=ds(1,2)*v(231)
    v(224)=ds(1,1)*v(230)
    v(285)=v(224)+v(225)+v(226)+v(227)+v(228)+v(229)
    v(148)=ds(2,1)*v(230)
    v(149)=ds(2,2)*v(231)
    v(150)=ds(2,3)*v(232)
    v(151)=ds(2,4)*v(233)
    v(152)=ds(2,5)*v(234)
    v(153)=ds(2,6)*v(235)
    v(281)=v(148)+v(149)+v(150)+v(151)+v(152)+v(153)
    v(154)=ds(3,1)*v(230)
    v(155)=ds(3,2)*v(231)
    v(156)=ds(3,3)*v(232)
    v(157)=ds(3,4)*v(233)
    v(158)=ds(3,5)*v(234)
    v(159)=ds(3,6)*v(235)
    v(282)=v(154)+v(155)+v(156)+v(157)+v(158)+v(159)
    v(163)=ds(5,1)*v(239)
    v(164)=ds(5,2)*v(240)
    v(165)=ds(5,3)*v(241)
    v(283)=dcl(4,1)*ds(5,4)+dcl(5,1)*ds(5,5)+dcl(6,1)*ds(5,6)+v(163)+v(164)+v(165)
    v(166)=ds(6,1)*v(239)
    v(167)=ds(6,2)*v(240)
    v(168)=ds(6,3)*v(241)
    v(284)=dcl(4,1)*ds(6,4)+dcl(5,1)*ds(6,5)+dcl(6,1)*ds(6,6)+v(166)+v(167)+v(168)
    v(175)=ds(2,1)*v(248)
    v(176)=ds(2,2)*v(249)
    v(177)=ds(2,3)*v(250)
    v(178)=ds(2,4)*v(251)
    v(179)=ds(2,5)*v(252)
    v(180)=ds(2,6)*v(253)
    v(287)=v(175)+v(176)+v(177)+v(178)+v(179)+v(180)
    v(181)=ds(3,1)*v(248)
    v(182)=ds(3,2)*v(249)
    v(183)=ds(3,3)*v(250)
    v(184)=ds(3,4)*v(251)
    v(185)=ds(3,5)*v(252)
    v(186)=ds(3,6)*v(253)
    v(288)=v(181)+v(182)+v(183)+v(184)+v(185)+v(186)
    v(190)=ds(5,1)*v(257)
    v(191)=ds(5,2)*v(258)
    v(192)=ds(5,3)*v(259)
    v(289)=dcl(4,2)*ds(5,4)+dcl(5,2)*ds(5,5)+dcl(6,2)*ds(5,6)+v(190)+v(191)+v(192)
    v(193)=ds(6,1)*v(257)
    v(194)=ds(6,2)*v(258)
    v(195)=ds(6,3)*v(259)
    v(290)=dcl(4,2)*ds(6,4)+dcl(5,2)*ds(6,5)+dcl(6,2)*ds(6,6)+v(193)+v(194)+v(195)
    v(202)=ds(2,1)*v(266)
    v(203)=ds(2,2)*v(267)
    v(204)=ds(2,3)*v(268)
    v(205)=ds(2,4)*v(269)
    v(206)=ds(2,5)*v(270)
    v(207)=ds(2,6)*v(271)
    v(293)=v(202)+v(203)+v(204)+v(205)+v(206)+v(207)
    v(208)=ds(3,1)*v(266)
    v(209)=ds(3,2)*v(267)
    v(210)=ds(3,3)*v(268)
    v(211)=ds(3,4)*v(269)
    v(212)=ds(3,5)*v(270)
    v(213)=ds(3,6)*v(271)
    v(294)=v(208)+v(209)+v(210)+v(211)+v(212)+v(213)
    v(217)=ds(5,1)*v(275)
    v(218)=ds(5,2)*v(276)
    v(219)=ds(5,3)*v(277)
    v(295)=dcl(4,3)*ds(5,4)+dcl(5,3)*ds(5,5)+dcl(6,3)*ds(5,6)+v(217)+v(218)+v(219)
    v(220)=ds(6,1)*v(275)
    v(221)=ds(6,2)*v(276)
    v(222)=ds(6,3)*v(277)
    v(296)=dcl(4,3)*ds(6,4)+dcl(5,3)*ds(6,5)+dcl(6,3)*ds(6,6)+v(220)+v(221)+v(222)
    kernel(1,1)=d2c(4,1,1)*s(4)+d2c(5,1,1)*s(5)+d2c(6,1,1)*s(6)+d2c(1,1,1)*v(278)+d2c(2,1,1)*v(279)+d2c(3,1,1)*v(280)+dck(2&
         &,1)*v(281)+dck(3,1)*v(282)+dck(5,1)*v(283)+dck(6,1)*v(284)+dck(1,1)*v(285)+dck(4,1)*v(286)
    kernel(1,2)=d2c(4,1,2)*s(4)+d2c(5,1,2)*s(5)+d2c(6,1,2)*s(6)+d2c(1,1,2)*v(278)+d2c(2,1,2)*v(279)+d2c(3,1,2)*v(280)+dck(2&
         &,1)*v(287)+dck(3,1)*v(288)+dck(5,1)*v(289)+dck(6,1)*v(290)+dck(1,1)*v(291)+dck(4,1)*v(292)
    kernel(1,3)=d2c(4,1,3)*s(4)+d2c(5,1,3)*s(5)+d2c(6,1,3)*s(6)+d2c(1,1,3)*v(278)+d2c(2,1,3)*v(279)+d2c(3,1,3)*v(280)+dck(2&
         &,1)*v(293)+dck(3,1)*v(294)+dck(5,1)*v(295)+dck(6,1)*v(296)+dck(1,1)*v(297)+dck(4,1)*v(298)
    kernel(2,1)=d2c(4,2,1)*s(4)+d2c(5,2,1)*s(5)+d2c(6,2,1)*s(6)+d2c(1,2,1)*v(278)+d2c(2,2,1)*v(279)+d2c(3,2,1)*v(280)+dck(2&
         &,2)*v(281)+dck(3,2)*v(282)+dck(5,2)*v(283)+dck(6,2)*v(284)+dck(1,2)*v(285)+dck(4,2)*v(286)
    kernel(2,2)=d2c(4,2,2)*s(4)+d2c(5,2,2)*s(5)+d2c(6,2,2)*s(6)+d2c(1,2,2)*v(278)+d2c(2,2,2)*v(279)+d2c(3,2,2)*v(280)+dck(2&
         &,2)*v(287)+dck(3,2)*v(288)+dck(5,2)*v(289)+dck(6,2)*v(290)+dck(1,2)*v(291)+dck(4,2)*v(292)
    kernel(2,3)=d2c(4,2,3)*s(4)+d2c(5,2,3)*s(5)+d2c(6,2,3)*s(6)+d2c(1,2,3)*v(278)+d2c(2,2,3)*v(279)+d2c(3,2,3)*v(280)+dck(2&
         &,2)*v(293)+dck(3,2)*v(294)+dck(5,2)*v(295)+dck(6,2)*v(296)+dck(1,2)*v(297)+dck(4,2)*v(298)
    kernel(3,1)=d2c(4,3,1)*s(4)+d2c(5,3,1)*s(5)+d2c(6,3,1)*s(6)+d2c(1,3,1)*v(278)+d2c(2,3,1)*v(279)+d2c(3,3,1)*v(280)+dck(2&
         &,3)*v(281)+dck(3,3)*v(282)+dck(5,3)*v(283)+dck(6,3)*v(284)+dck(1,3)*v(285)+dck(4,3)*v(286)
    kernel(3,2)=d2c(4,3,2)*s(4)+d2c(5,3,2)*s(5)+d2c(6,3,2)*s(6)+d2c(1,3,2)*v(278)+d2c(2,3,2)*v(279)+d2c(3,3,2)*v(280)+dck(2&
         &,3)*v(287)+dck(3,3)*v(288)+dck(5,3)*v(289)+dck(6,3)*v(290)+dck(1,3)*v(291)+dck(4,3)*v(292)
    kernel(3,3)=d2c(4,3,3)*s(4)+d2c(5,3,3)*s(5)+d2c(6,3,3)*s(6)+d2c(1,3,3)*v(278)+d2c(2,3,3)*v(279)+d2c(3,3,3)*v(280)+dck(2&
         &,3)*v(293)+dck(3,3)*v(294)+dck(5,3)*v(295)+dck(6,3)*v(296)+dck(1,3)*v(297)+dck(4,3)*v(298)
  END SUBROUTINE mls_fbar

SUBROUTINE mls_fbar_2d(dck,dcl,d2c,s,ds,kernel)
IMPLICIT NONE
DOUBLE PRECISION v(86),dck(3,2),dcl(3,2),d2c(3,2,2),s(3),ds(3,3),kernel(2,2)
v(75)=s(2)/2d0
v(74)=s(1)/2d0
v(73)=ds(3,2)/2d0
v(72)=ds(3,1)/2d0
v(71)=dcl(3,2)/2d0
v(70)=dcl(2,2)/4d0
v(69)=dcl(1,2)/4d0
v(68)=ds(1,3)*v(71)
v(67)=ds(1,2)*v(70)
v(66)=ds(1,1)*v(69)
v(81)=v(66)+v(67)+v(68)
v(65)=dcl(2,1)*v(73)
v(64)=dcl(1,1)*v(72)
v(78)=dcl(3,1)*ds(3,3)+v(64)+v(65)
v(63)=dcl(3,1)/2d0
v(62)=dcl(2,1)/4d0
v(61)=dcl(1,1)/4d0
v(60)=ds(1,3)*v(63)
v(59)=ds(1,2)*v(62)
v(58)=ds(1,1)*v(61)
v(77)=v(58)+v(59)+v(60)
v(44)=ds(2,1)*v(61)
v(45)=ds(2,2)*v(62)
v(46)=ds(2,3)*v(63)
v(76)=v(44)+v(45)+v(46)
v(52)=ds(2,1)*v(69)
v(53)=ds(2,2)*v(70)
v(54)=ds(2,3)*v(71)
v(79)=v(52)+v(53)+v(54)
v(55)=dcl(1,2)*v(72)
v(56)=dcl(2,2)*v(73)
v(80)=dcl(3,2)*ds(3,3)+v(55)+v(56)
kernel(1,1)=d2c(3,1,1)*s(3)+d2c(1,1,1)*v(74)+d2c(2,1,1)*v(75)+dck(2,1)*v(76)+dck(1,1)*v(77)+dck(3,1)*v(78)
kernel(1,2)=d2c(3,1,2)*s(3)+d2c(1,1,2)*v(74)+d2c(2,1,2)*v(75)+dck(2,1)*v(79)+dck(3,1)*v(80)+dck(1,1)*v(81)
kernel(2,1)=d2c(3,2,1)*s(3)+d2c(1,2,1)*v(74)+d2c(2,2,1)*v(75)+dck(2,2)*v(76)+dck(1,2)*v(77)+dck(3,2)*v(78)
kernel(2,2)=d2c(3,2,2)*s(3)+d2c(1,2,2)*v(74)+d2c(2,2,2)*v(75)+dck(2,2)*v(79)+dck(3,2)*v(80)+dck(1,2)*v(81)
END


  SUBROUTINE mls_classicaldc(nk,f,dc)
    IMPLICIT NONE
    DOUBLE PRECISION v(57),nk(3),f(3,3),dc(6,3)
    v(52)=2d0*nk(3)
    v(51)=2d0*nk(2)
    v(50)=2d0*nk(1)
    dc(1,1)=f(1,1)*v(50)
    dc(1,2)=f(2,1)*v(50)
    dc(1,3)=f(3,1)*v(50)
    dc(2,1)=f(1,2)*v(51)
    dc(2,2)=f(2,2)*v(51)
    dc(2,3)=f(3,2)*v(51)
    dc(3,1)=f(1,3)*v(52)
    dc(3,2)=f(2,3)*v(52)
    dc(3,3)=f(3,3)*v(52)
    dc(4,1)=f(1,2)*nk(1)+f(1,1)*nk(2)
    dc(4,2)=f(2,2)*nk(1)+f(2,1)*nk(2)
    dc(4,3)=f(3,2)*nk(1)+f(3,1)*nk(2)
    dc(5,1)=f(1,3)*nk(1)+f(1,1)*nk(3)
    dc(5,2)=f(2,3)*nk(1)+f(2,1)*nk(3)
    dc(5,3)=f(3,3)*nk(1)+f(3,1)*nk(3)
    dc(6,1)=f(1,3)*nk(2)+f(1,2)*nk(3)
    dc(6,2)=f(2,3)*nk(2)+f(2,2)*nk(3)
    dc(6,3)=f(3,3)*nk(2)+f(3,2)*nk(3)
  END SUBROUTINE mls_classicaldc


  SUBROUTINE mls_classicaldc_2d(nk,f,dc)
    IMPLICIT NONE
    DOUBLE PRECISION v(26),nk(2),f(2,2),dc(3,2)
    v(21)=2d0*nk(2)
    v(20)=2d0*nk(1)
    dc(1,1)=f(1,1)*v(20)
    dc(1,2)=f(2,1)*v(20)
    dc(2,1)=f(1,2)*v(21)
    dc(2,2)=f(2,2)*v(21)
    dc(3,1)=f(1,2)*nk(1)+f(1,1)*nk(2)
    dc(3,2)=f(2,2)*nk(1)+f(2,1)*nk(2)
  END SUBROUTINE mls_classicaldc_2d

  SUBROUTINE mls_classicaldc2(nk,nl,d2c)
    IMPLICIT NONE
    DOUBLE PRECISION v(78),nk(3),nl(3),d2c(6,3,3)
    v(73)=nk(3)*nl(2)+nk(2)*nl(3)
    v(72)=nk(3)*nl(1)+nk(1)*nl(3)
    v(71)=nk(2)*nl(1)+nk(1)*nl(2)
    v(70)=2d0*nk(3)*nl(3)
    v(69)=2d0*nk(2)*nl(2)
    v(68)=2d0*nk(1)*nl(1)
    d2c(1,1,1)=v(68)
    d2c(1,1,2)=0d0
    d2c(1,1,3)=0d0
    d2c(1,2,1)=0d0
    d2c(1,2,2)=v(68)
    d2c(1,2,3)=0d0
    d2c(1,3,1)=0d0
    d2c(1,3,2)=0d0
    d2c(1,3,3)=v(68)
    d2c(2,1,1)=v(69)
    d2c(2,1,2)=0d0
    d2c(2,1,3)=0d0
    d2c(2,2,1)=0d0
    d2c(2,2,2)=v(69)
    d2c(2,2,3)=0d0
    d2c(2,3,1)=0d0
    d2c(2,3,2)=0d0
    d2c(2,3,3)=v(69)
    d2c(3,1,1)=v(70)
    d2c(3,1,2)=0d0
    d2c(3,1,3)=0d0
    d2c(3,2,1)=0d0
    d2c(3,2,2)=v(70)
    d2c(3,2,3)=0d0
    d2c(3,3,1)=0d0
    d2c(3,3,2)=0d0
    d2c(3,3,3)=v(70)
    d2c(4,1,1)=v(71)
    d2c(4,1,2)=0d0
    d2c(4,1,3)=0d0
    d2c(4,2,1)=0d0
    d2c(4,2,2)=v(71)
    d2c(4,2,3)=0d0
    d2c(4,3,1)=0d0
    d2c(4,3,2)=0d0
    d2c(4,3,3)=v(71)
    d2c(5,1,1)=v(72)
    d2c(5,1,2)=0d0
    d2c(5,1,3)=0d0
    d2c(5,2,1)=0d0
    d2c(5,2,2)=v(72)
    d2c(5,2,3)=0d0
    d2c(5,3,1)=0d0
    d2c(5,3,2)=0d0
    d2c(5,3,3)=v(72)
    d2c(6,1,1)=v(73)
    d2c(6,1,2)=0d0
    d2c(6,1,3)=0d0
    d2c(6,2,1)=0d0
    d2c(6,2,2)=v(73)
    d2c(6,2,3)=0d0
    d2c(6,3,1)=0d0
    d2c(6,3,2)=0d0
    d2c(6,3,3)=v(73)
  END SUBROUTINE mls_classicaldc2


  SUBROUTINE mls_classicaldc2_2d(nk,nl,d2c)
    IMPLICIT NONE
    DOUBLE PRECISION v(28),nk(2),nl(2),d2c(3,2,2)
    v(23)=nk(2)*nl(1)+nk(1)*nl(2)
    v(22)=2d0*nk(2)*nl(2)
    v(21)=2d0*nk(1)*nl(1)
    d2c(1,1,1)=v(21)
    d2c(1,1,2)=0d0
    d2c(1,2,1)=0d0
    d2c(1,2,2)=v(21)
    d2c(2,1,1)=v(22)
    d2c(2,1,2)=0d0
    d2c(2,2,1)=0d0
    d2c(2,2,2)=v(22)
    d2c(3,1,1)=v(23)
    d2c(3,1,2)=0d0
    d2c(3,2,1)=0d0
    d2c(3,2,2)=v(23)
  END SUBROUTINE mls_classicaldc2_2d

  
!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           5 Feb 22 20:02:45  *
!**************************************************************
! User     : Full professional version
! Notebook : mls_combinedfbar0
! Evaluation time                 : 3 s     Mode  : Optimal
! Number of formulae              : 2       Method: Automatic
! Subroutine                      : mls_combinedfbar0 size: 203
! Total size of Mathematica  code : 203 subexpressions
! Total size of Fortran code      : 472 bytes

!******************* S U B R O U T I N E **********************
  SUBROUTINE mls_combinedfbar0(c,cb,cs)
    IMPLICIT NONE
    DOUBLE PRECISION v(32),c(6),cb(6),cs(6)
    v(27)=((cb(3)*cb(4)**2+cb(5)*(cb(2)*cb(5)-2d0*cb(4)*cb(6))+cb(1)*(-(cb(2)*cb(3))+cb(6)**2))/(c(3)*c(4)**2+c(5)*(c(2)*c&
         &(5)-2d0*c(4)*c(6))+c(1)*(-(c(2)*c(3))+c(6)**2)))**0.3333333333333333d0
    cs(1)=c(1)*v(27)
    cs(2)=c(2)*v(27)
    cs(3)=c(3)*v(27)
    cs(4)=c(4)*v(27)
    cs(5)=c(5)*v(27)
    cs(6)=c(6)*v(27)
  END SUBROUTINE mls_combinedfbar0

  SUBROUTINE mls_combinedfbar0_2d(c,cb,cs)
    IMPLICIT NONE
    DOUBLE PRECISION v(20),c(3),cb(3),cs(3)
    v(15)=SQRT((cb(1)*cb(2)-cb(3)**2)/(c(1)*c(2)-c(3)**2))
    cs(1)=c(1)*v(15)
    cs(2)=c(2)*v(15)
    cs(3)=c(3)*v(15)
  END SUBROUTINE mls_combinedfbar0_2d

  SUBROUTINE mls_combinedfbar1(c,cb,dc,dcb,dcs)
    IMPLICIT NONE
    DOUBLE PRECISION v(166),c(6),cb(6),dc(6,3),dcb(6,3),dcs(6,3)
    v(161)=cb(2)*dcb(3,3)
    v(160)=cb(3)*dcb(2,3)
    v(159)=-(cb(2)*dcb(5,3))
    v(158)=cb(6)*dcb(4,3)
    v(157)=c(2)*dc(3,3)
    v(156)=c(3)*dc(2,3)
    v(155)=-(c(2)*dc(5,3))
    v(154)=c(6)*dc(4,3)
    v(153)=2d0*dcb(6,3)
    v(152)=2d0*dc(6,3)
    v(151)=cb(2)*dcb(3,2)
    v(150)=cb(3)*dcb(2,2)
    v(149)=-(cb(2)*dcb(5,2))
    v(148)=cb(6)*dcb(4,2)
    v(147)=c(2)*dc(3,2)
    v(146)=c(3)*dc(2,2)
    v(145)=-(c(2)*dc(5,2))
    v(144)=c(6)*dc(4,2)
    v(143)=2d0*dcb(6,2)
    v(142)=cb(4)**2
    v(141)=cb(5)**2
    v(140)=2d0*dc(6,2)
    v(139)=c(4)**2
    v(138)=c(5)**2
    v(135)=cb(2)*dcb(3,1)
    v(134)=cb(3)*dcb(2,1)
    v(133)=-(cb(2)*dcb(5,1))
    v(132)=cb(6)*dcb(4,1)
    v(131)=-(dcb(3,1)*v(142))
    v(130)=-(dcb(2,1)*v(141))
    v(129)=c(2)*dc(3,1)
    v(128)=c(3)*dc(2,1)
    v(127)=-(c(2)*dc(5,1))
    v(126)=c(6)*dc(4,1)
    v(125)=-(dc(3,1)*v(139))
    v(124)=-(dc(2,1)*v(138))
    v(123)=2d0*dcb(6,1)
    v(122)=cb(4)*cb(6)
    v(121)=cb(3)*cb(4)
    v(120)=2d0*dc(6,1)
    v(119)=c(4)*c(6)
    v(118)=c(3)*c(4)
    v(117)=cb(2)*cb(3)-cb(6)**2
    v(116)=-(cb(5)*cb(6))+v(121)
    v(115)=-(cb(2)*cb(5))+v(122)
    v(114)=c(2)*c(3)-c(6)**2
    v(113)=-(c(5)*c(6))+v(118)
    v(112)=-(c(2)*c(5))+v(119)
    v(78)=c(5)*v(112)-c(4)*v(113)+c(1)*v(114)
    v(88)=1d0/v(78)**2
    v(79)=cb(5)*v(115)-cb(4)*v(116)+cb(1)*v(117)
    v(137)=-(v(79)*v(88))
    v(77)=v(79)/v(78)
    v(87)=1d0/v(77)**0.6666666666666666d0
    v(136)=v(87)/3d0
    v(71)=v(77)**0.3333333333333333d0
    v(86)=v(136)*((dc(5,1)*v(112)-dc(4,1)*v(113)+dc(1,1)*v(114)-dc(4,1)*v(118)+dc(5,1)*v(119)+v(124)+v(125)+c(5)*(c(4)*v&
         &(120)+v(126)+v(127))+c(1)*(-(c(6)*v(120))+v(128)+v(129)))*v(137)+(dcb(5,1)*v(115)-dcb(4,1)*v(116)+dcb(1,1)*v(117)-dcb(4&
         &,1)*v(121)+dcb(5,1)*v(122)+v(130)+v(131)+cb(5)*(cb(4)*v(123)+v(132)+v(133))+cb(1)*(-(cb(6)*v(123))+v(134)+v(135)))/v(78&
         &))
    v(90)=v(136)*(v(137)*(dc(5,2)*v(112)-dc(4,2)*v(113)+dc(1,2)*v(114)-dc(4,2)*v(118)+dc(5,2)*v(119)-dc(2,2)*v(138)-dc(3,2&
         &)*v(139)+c(5)*(c(4)*v(140)+v(144)+v(145))+c(1)*(-(c(6)*v(140))+v(146)+v(147)))+(dcb(5,2)*v(115)-dcb(4,2)*v(116)+dcb(1,2&
         &)*v(117)-dcb(4,2)*v(121)+dcb(5,2)*v(122)-dcb(2,2)*v(141)-dcb(3,2)*v(142)+cb(5)*(cb(4)*v(143)+v(148)+v(149))+cb(1)*(-(cb&
         &(6)*v(143))+v(150)+v(151)))/v(78))
    v(92)=v(136)*(v(137)*(dc(5,3)*v(112)-dc(4,3)*v(113)+dc(1,3)*v(114)-dc(4,3)*v(118)+dc(5,3)*v(119)-dc(2,3)*v(138)-dc(3,3&
         &)*v(139)+c(5)*(c(4)*v(152)+v(154)+v(155))+c(1)*(-(c(6)*v(152))+v(156)+v(157)))+(dcb(5,3)*v(115)-dcb(4,3)*v(116)+dcb(1,3&
         &)*v(117)-dcb(4,3)*v(121)+dcb(5,3)*v(122)-dcb(2,3)*v(141)-dcb(3,3)*v(142)+cb(5)*(cb(4)*v(153)+v(158)+v(159))+cb(1)*(-(cb&
         &(6)*v(153))+v(160)+v(161)))/v(78))
    dcs(1,1)=dc(1,1)*v(71)+c(1)*v(86)
    dcs(1,2)=dc(1,2)*v(71)+c(1)*v(90)
    dcs(1,3)=dc(1,3)*v(71)+c(1)*v(92)
    dcs(2,1)=dc(2,1)*v(71)+c(2)*v(86)
    dcs(2,2)=dc(2,2)*v(71)+c(2)*v(90)
    dcs(2,3)=dc(2,3)*v(71)+c(2)*v(92)
    dcs(3,1)=dc(3,1)*v(71)+c(3)*v(86)
    dcs(3,2)=dc(3,2)*v(71)+c(3)*v(90)
    dcs(3,3)=dc(3,3)*v(71)+c(3)*v(92)
    dcs(4,1)=dc(4,1)*v(71)+c(4)*v(86)
    dcs(4,2)=dc(4,2)*v(71)+c(4)*v(90)
    dcs(4,3)=dc(4,3)*v(71)+c(4)*v(92)
    dcs(5,1)=dc(5,1)*v(71)+c(5)*v(86)
    dcs(5,2)=dc(5,2)*v(71)+c(5)*v(90)
    dcs(5,3)=dc(5,3)*v(71)+c(5)*v(92)
    dcs(6,1)=dc(6,1)*v(71)+c(6)*v(86)
    dcs(6,2)=dc(6,2)*v(71)+c(6)*v(90)
    dcs(6,3)=dc(6,3)*v(71)+c(6)*v(92)
  END SUBROUTINE mls_combinedfbar1

  SUBROUTINE mls_combinedfbar1_2d(c,cb,dc,dcb,dcs)
    IMPLICIT NONE
    DOUBLE PRECISION v(62),c(3),cb(3),dc(3,2),dcb(3,2),dcs(3,2)
    v(57)=c(1)*dc(2,2)
    v(56)=c(2)*dc(1,2)
    v(55)=cb(1)*dcb(2,2)
    v(54)=cb(2)*dcb(1,2)
    v(53)=(-2d0)*c(3)
    v(51)=(-2d0)*cb(3)
    v(49)=c(2)*dc(1,1)+c(1)*dc(2,1)+dc(3,1)*v(53)
    v(48)=cb(2)*dcb(1,1)+cb(1)*dcb(2,1)+dcb(3,1)*v(51)
    v(47)=cb(1)*cb(2)-cb(3)**2
    v(46)=c(1)*c(2)-c(3)**2
    v(36)=1d0/v(46)**2
    v(52)=-(v(36)*v(47))
    v(31)=v(47)/v(46)
    v(35)=1d0/sqrt(v(31))
    v(50)=v(35)/2d0
    v(28)=sqrt(v(31))
    v(34)=v(50)*(v(48)/v(46)+v(49)*v(52))
    v(38)=v(50)*((dcb(3,2)*v(51)+v(54)+v(55))/v(46)+v(52)*(dc(3,2)*v(53)+v(56)+v(57)))
    dcs(1,1)=dc(1,1)*v(28)+c(1)*v(34)
    dcs(1,2)=dc(1,2)*v(28)+c(1)*v(38)
    dcs(2,1)=dc(2,1)*v(28)+c(2)*v(34)
    dcs(2,2)=dc(2,2)*v(28)+c(2)*v(38)
    dcs(3,1)=dc(3,1)*v(28)+c(3)*v(34)
    dcs(3,2)=dc(3,2)*v(28)+c(3)*v(38)
  END SUBROUTINE mls_combinedfbar1_2d

  SUBROUTINE mls_combinedfbar2_2d(c1,c0,dck1,dcl1,dck0,dcl0,d2c1,d2c0,d2c)
    IMPLICIT NONE
    DOUBLE PRECISION v(202),c1(3),c0(3),dck1(3,2),dcl1(3,2),dck0(3,2),dcl0(3,2),d2c1(3,2,2),d2c0(3,2,2),d2c(3,2,2&
         &)
    v(197)=dck1(1,1)*dcl1(2,2)
    v(196)=dck1(2,1)*dcl1(1,2)
    v(195)=dck0(1,1)*dcl0(2,2)
    v(194)=dck0(2,1)*dcl0(1,2)
    v(193)=c0(1)*d2c0(2,1,2)
    v(192)=c0(2)*d2c0(1,1,2)
    v(191)=(-2d0)*dcl1(3,2)
    v(188)=(-2d0)*dcl0(3,2)
    v(187)=dck1(3,2)*v(191)
    v(186)=dck1(1,2)*dcl1(2,2)
    v(185)=dck1(2,2)*dcl1(1,2)
    v(184)=dck0(3,2)*v(188)
    v(183)=dck0(1,2)*dcl0(2,2)
    v(182)=dck0(2,2)*dcl0(1,2)
    v(181)=c0(1)*d2c0(2,2,2)
    v(180)=c0(2)*d2c0(1,2,2)
    v(178)=dck1(1,1)*dcl1(2,1)
    v(177)=dck1(2,1)*dcl1(1,1)
    v(176)=dck0(1,1)*dcl0(2,1)
    v(175)=dck0(2,1)*dcl0(1,1)
    v(174)=c0(1)*d2c0(2,1,1)
    v(173)=c0(2)*d2c0(1,1,1)
    v(172)=(-2d0)*dcl1(3,1)
    v(169)=(-2d0)*dcl0(3,1)
    v(168)=dck1(3,2)*v(172)
    v(167)=dck1(1,2)*dcl1(2,1)
    v(166)=dck1(2,2)*dcl1(1,1)
    v(165)=dck0(3,2)*v(169)
    v(164)=dck0(1,2)*dcl0(2,1)
    v(163)=dck0(2,2)*dcl0(1,1)
    v(162)=c0(1)*d2c0(2,2,1)
    v(161)=c0(2)*d2c0(1,2,1)
    v(160)=c0(1)*dcl0(2,2)
    v(159)=c0(2)*dcl0(1,2)
    v(158)=c0(1)*dcl0(2,1)
    v(157)=c0(2)*dcl0(1,1)
    v(155)=c1(1)*dcl1(2,2)
    v(154)=c1(2)*dcl1(1,2)
    v(153)=c1(1)*dcl1(2,1)
    v(152)=c1(2)*dcl1(1,1)
    v(150)=d2c1(3,2,2)
    v(149)=d2c1(3,2,1)
    v(148)=d2c1(3,1,2)
    v(147)=d2c1(3,1,1)
    v(146)=d2c1(2,2,2)
    v(145)=d2c1(2,2,1)
    v(144)=d2c1(2,1,2)
    v(143)=d2c1(2,1,1)
    v(142)=d2c1(1,2,2)
    v(141)=d2c1(1,2,1)
    v(140)=d2c1(1,1,2)
    v(139)=d2c1(1,1,1)
    v(137)=c0(1)*dck0(2,2)
    v(136)=c0(2)*dck0(1,2)
    v(135)=(-2d0)*c0(3)
    v(134)=c0(2)*dck0(1,1)+c0(1)*dck0(2,1)+dck0(3,1)*v(135)
    v(133)=c1(1)*dck1(2,2)
    v(132)=c1(2)*dck1(1,2)
    v(131)=(-2d0)*c1(3)
    v(130)=c1(2)*dck1(1,1)+c1(1)*dck1(2,1)+dck1(3,1)*v(131)
    v(129)=c0(1)*c0(2)-c0(3)**2
    v(128)=c1(1)*c1(2)-c1(3)**2
    v(94)=1d0/v(128)**3
    v(156)=(-2d0)*v(94)
    v(80)=1d0/v(128)**2
    v(138)=-(v(129)*v(80))
    v(75)=v(129)/v(128)
    v(101)=1d0/v(75)**0.15d1
    v(179)=(-0.5d0)*v(101)
    v(79)=1d0/sqrt(v(75))
    v(151)=v(79)/2d0
    v(72)=sqrt(v(75))
    v(111)=dck1(3,2)*v(131)+v(132)+v(133)
    v(107)=v(134)/v(128)+v(130)*v(138)
    v(110)=dck0(3,2)*v(135)+v(136)+v(137)
    v(112)=v(110)/v(128)+v(111)*v(138)
    v(78)=v(107)*v(151)
    v(82)=v(112)*v(151)
    v(89)=dcl1(3,1)*v(131)+v(152)+v(153)
    v(90)=dcl1(3,2)*v(131)+v(154)+v(155)
    v(91)=-(v(80)*v(89))
    v(92)=-(v(80)*v(90))
    v(93)=v(156)*v(89)
    v(171)=-(v(129)*v(93))
    v(95)=v(156)*v(90)
    v(190)=-(v(129)*v(95))
    v(96)=dcl0(3,1)*v(135)+v(157)+v(158)
    v(170)=-(v(80)*v(96))
    v(97)=dcl0(3,2)*v(135)+v(159)+v(160)
    v(189)=-(v(80)*v(97))
    v(98)=v(129)*v(91)+v(96)/v(128)
    v(99)=v(129)*v(92)+v(97)/v(128)
    v(100)=v(179)*v(98)
    v(113)=(v(100)*v(112)+v(79)*((d2c0(3,2,1)*v(135)+v(161)+v(162)+v(163)+v(164)+v(165))/v(128)+v(138)*(c1(2)*v(141)+c1(1&
         &)*v(145)+v(131)*v(149)+v(166)+v(167)+v(168))+v(111)*v(170)+v(111)*v(171)+v(110)*v(91)))/2d0
    v(108)=(v(100)*v(107)+v(79)*(v(130)*v(170)+v(130)*v(171)+(d2c0(3,1,1)*v(135)+dck0(3,1)*v(169)+v(173)+v(174)+v(175)+v&
         &(176))/v(128)+v(138)*(c1(2)*v(139)+c1(1)*v(143)+v(131)*v(147)+dck1(3,1)*v(172)+v(177)+v(178))+v(134)*v(91)))/2d0
    v(102)=v(179)*v(99)
    v(114)=(v(102)*v(112)+v(79)*((d2c0(3,2,2)*v(135)+v(180)+v(181)+v(182)+v(183)+v(184))/v(128)+v(138)*(c1(2)*v(142)+c1(1&
         &)*v(146)+v(131)*v(150)+v(185)+v(186)+v(187))+v(111)*v(189)+v(111)*v(190)+v(110)*v(92)))/2d0
    v(109)=(v(102)*v(107)+v(79)*(v(130)*v(189)+v(130)*v(190)+(d2c0(3,1,2)*v(135)+dck0(3,1)*v(188)+v(192)+v(193)+v(194)+v&
         &(195))/v(128)+v(138)*(c1(2)*v(140)+c1(1)*v(144)+v(131)*v(148)+dck1(3,1)*v(191)+v(196)+v(197))+v(134)*v(92)))/2d0
    v(103)=v(151)*v(98)
    v(104)=v(151)*v(99)
    d2c(1,1,1)=dck1(1,1)*v(103)+c1(1)*v(108)+v(139)*v(72)+dcl1(1,1)*v(78)
    d2c(1,1,2)=dck1(1,1)*v(104)+c1(1)*v(109)+v(140)*v(72)+dcl1(1,2)*v(78)
    d2c(1,2,1)=dck1(1,2)*v(103)+c1(1)*v(113)+v(141)*v(72)+dcl1(1,1)*v(82)
    d2c(1,2,2)=dck1(1,2)*v(104)+c1(1)*v(114)+v(142)*v(72)+dcl1(1,2)*v(82)
    d2c(2,1,1)=dck1(2,1)*v(103)+c1(2)*v(108)+v(143)*v(72)+dcl1(2,1)*v(78)
    d2c(2,1,2)=dck1(2,1)*v(104)+c1(2)*v(109)+v(144)*v(72)+dcl1(2,2)*v(78)
    d2c(2,2,1)=dck1(2,2)*v(103)+c1(2)*v(113)+v(145)*v(72)+dcl1(2,1)*v(82)
    d2c(2,2,2)=dck1(2,2)*v(104)+c1(2)*v(114)+v(146)*v(72)+dcl1(2,2)*v(82)
    d2c(3,1,1)=dck1(3,1)*v(103)+c1(3)*v(108)+v(147)*v(72)+dcl1(3,1)*v(78)
    d2c(3,1,2)=dck1(3,1)*v(104)+c1(3)*v(109)+v(148)*v(72)+dcl1(3,2)*v(78)
    d2c(3,2,1)=dck1(3,2)*v(103)+c1(3)*v(113)+v(149)*v(72)+dcl1(3,1)*v(82)
    d2c(3,2,2)=dck1(3,2)*v(104)+c1(3)*v(114)+v(150)*v(72)+dcl1(3,2)*v(82)
  END SUBROUTINE mls_combinedfbar2_2d

!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           9 Feb 22 10:52:10  *
!**************************************************************
! User     : Full professional version
! Notebook : mls_combinedfbar2
! Evaluation time                 : 25 s    Mode  : Optimal
! Number of formulae              : 416     Method: Automatic
! Subroutine                      : mls_combinedfbar2 size: 13062
! Total size of Mathematica  code : 13062 subexpressions
! Total size of Fortran code      : 24811 bytes

!******************* S U B R O U T I N E **********************
  SUBROUTINE mls_combinedfbar2(c1,c0,dck1,dcl1,dck0,dcl0,d2c1,d2c0,d2c)
    IMPLICIT NONE
    DOUBLE PRECISION v(785),c1(6),c0(6),dck1(6,3),dcl1(6,3),dck0(6,3),dcl0(6,3),d2c1(6,3,3),d2c0(6,3,3),d2c(6,3,3&
         &)
    v(780)=dck0(2,1)*dcl0(3,3)
    v(779)=dck0(3,1)*dcl0(2,3)
    v(778)=-(dck0(5,1)*dcl0(6,3))
    v(777)=-(dck0(6,1)*dcl0(5,3))
    v(776)=dck0(3,1)*dcl0(4,3)
    v(775)=dck0(4,1)*dcl0(3,3)
    v(774)=dck0(4,1)*dcl0(6,3)
    v(773)=-(dck0(2,1)*dcl0(5,3))
    v(772)=dck0(6,1)*dcl0(4,3)
    v(771)=-(dck0(5,1)*dcl0(2,3))
    v(770)=dck1(2,1)*dcl1(3,3)
    v(769)=dck1(3,1)*dcl1(2,3)
    v(768)=-(dck1(5,1)*dcl1(6,3))
    v(767)=-(dck1(6,1)*dcl1(5,3))
    v(766)=dck1(3,1)*dcl1(4,3)
    v(765)=dck1(4,1)*dcl1(3,3)
    v(764)=dck1(4,1)*dcl1(6,3)
    v(763)=-(dck1(2,1)*dcl1(5,3))
    v(762)=dck1(6,1)*dcl1(4,3)
    v(761)=-(dck1(5,1)*dcl1(2,3))
    v(760)=dck0(2,2)*dcl0(3,3)
    v(759)=dck0(3,2)*dcl0(2,3)
    v(758)=-(dck0(5,2)*dcl0(6,3))
    v(757)=-(dck0(6,2)*dcl0(5,3))
    v(756)=dck0(3,2)*dcl0(4,3)
    v(755)=dck0(4,2)*dcl0(3,3)
    v(754)=dck0(4,2)*dcl0(6,3)
    v(753)=-(dck0(2,2)*dcl0(5,3))
    v(752)=dck0(6,2)*dcl0(4,3)
    v(751)=-(dck0(5,2)*dcl0(2,3))
    v(750)=dck1(2,2)*dcl1(3,3)
    v(749)=dck1(3,2)*dcl1(2,3)
    v(748)=-(dck1(5,2)*dcl1(6,3))
    v(747)=-(dck1(6,2)*dcl1(5,3))
    v(746)=dck1(3,2)*dcl1(4,3)
    v(745)=dck1(4,2)*dcl1(3,3)
    v(744)=dck1(4,2)*dcl1(6,3)
    v(743)=-(dck1(2,2)*dcl1(5,3))
    v(742)=dck1(6,2)*dcl1(4,3)
    v(741)=-(dck1(5,2)*dcl1(2,3))
    v(740)=(-2d0)*dcl0(6,3)
    v(739)=(-2d0)*dcl1(6,3)
    v(736)=dck0(6,3)*v(740)
    v(735)=dck0(2,3)*dcl0(3,3)
    v(734)=dck0(3,3)*dcl0(2,3)
    v(733)=-(dck0(5,3)*dcl0(6,3))
    v(732)=-(dck0(6,3)*dcl0(5,3))
    v(731)=dck0(3,3)*dcl0(4,3)
    v(730)=dck0(4,3)*dcl0(3,3)
    v(729)=dck0(4,3)*dcl0(6,3)
    v(728)=-(dck0(2,3)*dcl0(5,3))
    v(727)=dck0(6,3)*dcl0(4,3)
    v(726)=-(dck0(5,3)*dcl0(2,3))
    v(725)=dck1(6,3)*v(739)
    v(724)=dck1(2,3)*dcl1(3,3)
    v(723)=dck1(3,3)*dcl1(2,3)
    v(722)=-(dck1(5,3)*dcl1(6,3))
    v(721)=-(dck1(6,3)*dcl1(5,3))
    v(720)=dck1(3,3)*dcl1(4,3)
    v(719)=dck1(4,3)*dcl1(3,3)
    v(718)=dck1(4,3)*dcl1(6,3)
    v(717)=-(dck1(2,3)*dcl1(5,3))
    v(716)=dck1(6,3)*dcl1(4,3)
    v(715)=-(dck1(5,3)*dcl1(2,3))
    v(714)=dck0(2,1)*dcl0(3,2)
    v(713)=dck0(3,1)*dcl0(2,2)
    v(712)=-(dck0(5,1)*dcl0(6,2))
    v(711)=-(dck0(6,1)*dcl0(5,2))
    v(710)=dck0(3,1)*dcl0(4,2)
    v(709)=dck0(4,1)*dcl0(3,2)
    v(708)=dck0(4,1)*dcl0(6,2)
    v(707)=-(dck0(2,1)*dcl0(5,2))
    v(706)=dck0(6,1)*dcl0(4,2)
    v(705)=-(dck0(5,1)*dcl0(2,2))
    v(704)=dck1(2,1)*dcl1(3,2)
    v(703)=dck1(3,1)*dcl1(2,2)
    v(702)=-(dck1(5,1)*dcl1(6,2))
    v(701)=-(dck1(6,1)*dcl1(5,2))
    v(700)=dck1(3,1)*dcl1(4,2)
    v(699)=dck1(4,1)*dcl1(3,2)
    v(698)=dck1(4,1)*dcl1(6,2)
    v(697)=-(dck1(2,1)*dcl1(5,2))
    v(696)=dck1(6,1)*dcl1(4,2)
    v(695)=-(dck1(5,1)*dcl1(2,2))
    v(694)=dck0(2,2)*dcl0(3,2)
    v(693)=dck0(3,2)*dcl0(2,2)
    v(692)=-(dck0(5,2)*dcl0(6,2))
    v(691)=-(dck0(6,2)*dcl0(5,2))
    v(690)=dck0(3,2)*dcl0(4,2)
    v(689)=dck0(4,2)*dcl0(3,2)
    v(688)=dck0(4,2)*dcl0(6,2)
    v(687)=-(dck0(2,2)*dcl0(5,2))
    v(686)=dck0(6,2)*dcl0(4,2)
    v(685)=-(dck0(5,2)*dcl0(2,2))
    v(684)=dck1(2,2)*dcl1(3,2)
    v(683)=dck1(3,2)*dcl1(2,2)
    v(682)=-(dck1(5,2)*dcl1(6,2))
    v(681)=-(dck1(6,2)*dcl1(5,2))
    v(680)=dck1(3,2)*dcl1(4,2)
    v(679)=dck1(4,2)*dcl1(3,2)
    v(678)=dck1(4,2)*dcl1(6,2)
    v(677)=-(dck1(2,2)*dcl1(5,2))
    v(676)=dck1(6,2)*dcl1(4,2)
    v(675)=-(dck1(5,2)*dcl1(2,2))
    v(674)=(-2d0)*dcl0(6,2)
    v(673)=(-2d0)*dcl1(6,2)
    v(670)=dck0(6,3)*v(674)
    v(669)=dck0(2,3)*dcl0(3,2)
    v(668)=dck0(3,3)*dcl0(2,2)
    v(667)=-(dck0(5,3)*dcl0(6,2))
    v(666)=-(dck0(6,3)*dcl0(5,2))
    v(665)=dck0(3,3)*dcl0(4,2)
    v(664)=dck0(4,3)*dcl0(3,2)
    v(663)=dck0(4,3)*dcl0(6,2)
    v(662)=-(dck0(2,3)*dcl0(5,2))
    v(661)=dck0(6,3)*dcl0(4,2)
    v(660)=-(dck0(5,3)*dcl0(2,2))
    v(659)=dck1(6,3)*v(673)
    v(658)=dck1(2,3)*dcl1(3,2)
    v(657)=dck1(3,3)*dcl1(2,2)
    v(656)=-(dck1(5,3)*dcl1(6,2))
    v(655)=-(dck1(6,3)*dcl1(5,2))
    v(654)=dck1(3,3)*dcl1(4,2)
    v(653)=dck1(4,3)*dcl1(3,2)
    v(652)=dck1(4,3)*dcl1(6,2)
    v(651)=-(dck1(2,3)*dcl1(5,2))
    v(650)=dck1(6,3)*dcl1(4,2)
    v(649)=-(dck1(5,3)*dcl1(2,2))
    v(647)=dck0(2,1)*dcl0(3,1)
    v(646)=dck0(3,1)*dcl0(2,1)
    v(645)=-(dck0(5,1)*dcl0(6,1))
    v(644)=-(dck0(6,1)*dcl0(5,1))
    v(643)=dck0(3,1)*dcl0(4,1)
    v(642)=dck0(4,1)*dcl0(3,1)
    v(641)=dck0(4,1)*dcl0(6,1)
    v(640)=-(dck0(2,1)*dcl0(5,1))
    v(639)=dck0(6,1)*dcl0(4,1)
    v(638)=-(dck0(5,1)*dcl0(2,1))
    v(637)=dck1(2,1)*dcl1(3,1)
    v(636)=dck1(3,1)*dcl1(2,1)
    v(635)=-(dck1(5,1)*dcl1(6,1))
    v(634)=-(dck1(6,1)*dcl1(5,1))
    v(633)=dck1(3,1)*dcl1(4,1)
    v(632)=dck1(4,1)*dcl1(3,1)
    v(631)=dck1(4,1)*dcl1(6,1)
    v(630)=-(dck1(2,1)*dcl1(5,1))
    v(629)=dck1(6,1)*dcl1(4,1)
    v(628)=-(dck1(5,1)*dcl1(2,1))
    v(627)=dck0(2,2)*dcl0(3,1)
    v(626)=dck0(3,2)*dcl0(2,1)
    v(625)=-(dck0(5,2)*dcl0(6,1))
    v(624)=-(dck0(6,2)*dcl0(5,1))
    v(623)=dck0(3,2)*dcl0(4,1)
    v(622)=dck0(4,2)*dcl0(3,1)
    v(621)=dck0(4,2)*dcl0(6,1)
    v(620)=-(dck0(2,2)*dcl0(5,1))
    v(619)=dck0(6,2)*dcl0(4,1)
    v(618)=-(dck0(5,2)*dcl0(2,1))
    v(617)=dck1(2,2)*dcl1(3,1)
    v(616)=dck1(3,2)*dcl1(2,1)
    v(615)=-(dck1(5,2)*dcl1(6,1))
    v(614)=-(dck1(6,2)*dcl1(5,1))
    v(613)=dck1(3,2)*dcl1(4,1)
    v(612)=dck1(4,2)*dcl1(3,1)
    v(611)=dck1(4,2)*dcl1(6,1)
    v(610)=-(dck1(2,2)*dcl1(5,1))
    v(609)=dck1(6,2)*dcl1(4,1)
    v(608)=-(dck1(5,2)*dcl1(2,1))
    v(607)=(-2d0)*dcl0(6,1)
    v(606)=(-2d0)*dcl1(6,1)
    v(603)=dck0(6,3)*v(607)
    v(602)=dck0(2,3)*dcl0(3,1)
    v(601)=dck0(3,3)*dcl0(2,1)
    v(600)=-(dck0(5,3)*dcl0(6,1))
    v(599)=-(dck0(6,3)*dcl0(5,1))
    v(598)=dck0(3,3)*dcl0(4,1)
    v(597)=dck0(4,3)*dcl0(3,1)
    v(596)=dck0(4,3)*dcl0(6,1)
    v(595)=-(dck0(2,3)*dcl0(5,1))
    v(594)=dck0(6,3)*dcl0(4,1)
    v(593)=-(dck0(5,3)*dcl0(2,1))
    v(592)=dck1(6,3)*v(606)
    v(591)=dck1(2,3)*dcl1(3,1)
    v(590)=dck1(3,3)*dcl1(2,1)
    v(589)=-(dck1(5,3)*dcl1(6,1))
    v(588)=-(dck1(6,3)*dcl1(5,1))
    v(587)=dck1(3,3)*dcl1(4,1)
    v(586)=dck1(4,3)*dcl1(3,1)
    v(585)=dck1(4,3)*dcl1(6,1)
    v(584)=-(dck1(2,3)*dcl1(5,1))
    v(583)=dck1(6,3)*dcl1(4,1)
    v(582)=-(dck1(5,3)*dcl1(2,1))
    v(581)=c0(2)*dcl0(3,3)
    v(580)=c0(3)*dcl0(2,3)
    v(579)=c0(2)*dcl0(3,2)
    v(578)=c0(3)*dcl0(2,2)
    v(577)=c0(2)*dcl0(3,1)
    v(576)=c0(3)*dcl0(2,1)
    v(575)=c0(4)*dcl0(3,3)+c0(3)*dcl0(4,3)-c0(6)*dcl0(5,3)-c0(5)*dcl0(6,3)
    v(574)=c0(4)*dcl0(3,2)+c0(3)*dcl0(4,2)-c0(6)*dcl0(5,2)-c0(5)*dcl0(6,2)
    v(573)=c0(4)*dcl0(3,1)+c0(3)*dcl0(4,1)-c0(6)*dcl0(5,1)-c0(5)*dcl0(6,1)
    v(572)=-(c0(5)*dcl0(2,3))+c0(6)*dcl0(4,3)-c0(2)*dcl0(5,3)+c0(4)*dcl0(6,3)
    v(571)=-(c0(5)*dcl0(2,2))+c0(6)*dcl0(4,2)-c0(2)*dcl0(5,2)+c0(4)*dcl0(6,2)
    v(570)=-(c0(5)*dcl0(2,1))+c0(6)*dcl0(4,1)-c0(2)*dcl0(5,1)+c0(4)*dcl0(6,1)
    v(568)=c1(2)*dcl1(3,3)
    v(567)=c1(3)*dcl1(2,3)
    v(566)=c1(2)*dcl1(3,2)
    v(565)=c1(3)*dcl1(2,2)
    v(564)=c1(2)*dcl1(3,1)
    v(563)=c1(3)*dcl1(2,1)
    v(562)=c1(4)*dcl1(3,3)+c1(3)*dcl1(4,3)-c1(6)*dcl1(5,3)-c1(5)*dcl1(6,3)
    v(561)=c1(4)*dcl1(3,2)+c1(3)*dcl1(4,2)-c1(6)*dcl1(5,2)-c1(5)*dcl1(6,2)
    v(560)=c1(4)*dcl1(3,1)+c1(3)*dcl1(4,1)-c1(6)*dcl1(5,1)-c1(5)*dcl1(6,1)
    v(559)=-(c1(5)*dcl1(2,3))+c1(6)*dcl1(4,3)-c1(2)*dcl1(5,3)+c1(4)*dcl1(6,3)
    v(558)=-(c1(5)*dcl1(2,2))+c1(6)*dcl1(4,2)-c1(2)*dcl1(5,2)+c1(4)*dcl1(6,2)
    v(557)=-(c1(5)*dcl1(2,1))+c1(6)*dcl1(4,1)-c1(2)*dcl1(5,1)+c1(4)*dcl1(6,1)
    v(555)=d2c0(6,3,3)
    v(554)=d2c0(6,3,2)
    v(553)=d2c0(6,3,1)
    v(552)=d2c0(6,2,3)
    v(551)=d2c0(6,2,2)
    v(550)=d2c0(6,2,1)
    v(549)=d2c0(6,1,3)
    v(548)=d2c0(6,1,2)
    v(547)=d2c0(6,1,1)
    v(546)=d2c0(5,3,3)
    v(545)=d2c0(5,3,2)
    v(544)=d2c0(5,3,1)
    v(543)=d2c0(5,2,3)
    v(542)=d2c0(5,2,2)
    v(541)=d2c0(5,2,1)
    v(540)=d2c0(5,1,3)
    v(539)=d2c0(5,1,2)
    v(538)=d2c0(5,1,1)
    v(537)=d2c0(4,3,3)
    v(536)=d2c0(4,3,2)
    v(535)=d2c0(4,3,1)
    v(534)=d2c0(4,2,3)
    v(533)=d2c0(4,2,2)
    v(532)=d2c0(4,2,1)
    v(531)=d2c0(4,1,3)
    v(530)=d2c0(4,1,2)
    v(529)=d2c0(4,1,1)
    v(528)=d2c0(3,3,3)
    v(527)=d2c0(3,3,2)
    v(526)=d2c0(3,3,1)
    v(525)=d2c0(3,2,3)
    v(524)=d2c0(3,2,2)
    v(523)=d2c0(3,2,1)
    v(522)=d2c0(3,1,3)
    v(521)=d2c0(3,1,2)
    v(520)=d2c0(3,1,1)
    v(519)=d2c0(2,3,3)
    v(518)=d2c0(2,3,2)
    v(517)=d2c0(2,3,1)
    v(516)=d2c0(2,2,3)
    v(515)=d2c0(2,2,2)
    v(514)=d2c0(2,2,1)
    v(513)=d2c0(2,1,3)
    v(512)=d2c0(2,1,2)
    v(511)=d2c0(2,1,1)
    v(510)=d2c1(6,3,3)
    v(509)=d2c1(6,3,2)
    v(508)=d2c1(6,3,1)
    v(507)=d2c1(6,2,3)
    v(506)=d2c1(6,2,2)
    v(505)=d2c1(6,2,1)
    v(504)=d2c1(6,1,3)
    v(503)=d2c1(6,1,2)
    v(502)=d2c1(6,1,1)
    v(501)=d2c1(5,3,3)
    v(500)=d2c1(5,3,2)
    v(499)=d2c1(5,3,1)
    v(498)=d2c1(5,2,3)
    v(497)=d2c1(5,2,2)
    v(496)=d2c1(5,2,1)
    v(495)=d2c1(5,1,3)
    v(494)=d2c1(5,1,2)
    v(493)=d2c1(5,1,1)
    v(492)=d2c1(4,3,3)
    v(491)=d2c1(4,3,2)
    v(490)=d2c1(4,3,1)
    v(489)=d2c1(4,2,3)
    v(488)=d2c1(4,2,2)
    v(487)=d2c1(4,2,1)
    v(486)=d2c1(4,1,3)
    v(485)=d2c1(4,1,2)
    v(484)=d2c1(4,1,1)
    v(483)=d2c1(3,3,3)
    v(482)=d2c1(3,3,2)
    v(481)=d2c1(3,3,1)
    v(480)=d2c1(3,2,3)
    v(479)=d2c1(3,2,2)
    v(478)=d2c1(3,2,1)
    v(477)=d2c1(3,1,3)
    v(476)=d2c1(3,1,2)
    v(475)=d2c1(3,1,1)
    v(474)=d2c1(2,3,3)
    v(473)=d2c1(2,3,2)
    v(472)=d2c1(2,3,1)
    v(471)=d2c1(2,2,3)
    v(470)=d2c1(2,2,2)
    v(469)=d2c1(2,2,1)
    v(468)=d2c1(2,1,3)
    v(467)=d2c1(2,1,2)
    v(466)=d2c1(2,1,1)
    v(465)=d2c1(1,3,3)
    v(464)=d2c1(1,3,2)
    v(463)=d2c1(1,3,1)
    v(462)=d2c1(1,2,3)
    v(461)=d2c1(1,2,2)
    v(460)=d2c1(1,2,1)
    v(459)=d2c1(1,1,3)
    v(458)=d2c1(1,1,2)
    v(457)=d2c1(1,1,1)
    v(456)=c0(2)*dck0(3,3)
    v(455)=c0(3)*dck0(2,3)
    v(454)=c0(4)*dck0(3,3)+c0(3)*dck0(4,3)-c0(6)*dck0(5,3)-c0(5)*dck0(6,3)
    v(453)=-(c0(5)*dck0(2,3))+c0(6)*dck0(4,3)-c0(2)*dck0(5,3)+c0(4)*dck0(6,3)
    v(451)=c0(2)*dck0(3,2)
    v(450)=c0(3)*dck0(2,2)
    v(449)=(-2d0)*c0(6)
    v(448)=c0(4)*dck0(3,2)+c0(3)*dck0(4,2)-c0(6)*dck0(5,2)-c0(5)*dck0(6,2)
    v(447)=-(c0(5)*dck0(2,2))+c0(6)*dck0(4,2)-c0(2)*dck0(5,2)+c0(4)*dck0(6,2)
    v(446)=c0(3)*dck0(2,1)+c0(2)*dck0(3,1)+dck0(6,1)*v(449)
    v(445)=c0(4)*dck0(3,1)+c0(3)*dck0(4,1)-c0(6)*dck0(5,1)-c0(5)*dck0(6,1)
    v(444)=-(c0(5)*dck0(2,1))+c0(6)*dck0(4,1)-c0(2)*dck0(5,1)+c0(4)*dck0(6,1)
    v(443)=c1(2)*dck1(3,3)
    v(442)=c1(3)*dck1(2,3)
    v(441)=c1(4)*dck1(3,3)+c1(3)*dck1(4,3)-c1(6)*dck1(5,3)-c1(5)*dck1(6,3)
    v(440)=-(c1(5)*dck1(2,3))+c1(6)*dck1(4,3)-c1(2)*dck1(5,3)+c1(4)*dck1(6,3)
    v(439)=c1(2)*dck1(3,2)
    v(438)=c1(3)*dck1(2,2)
    v(437)=(-2d0)*c1(6)
    v(436)=c1(4)*dck1(3,2)+c1(3)*dck1(4,2)-c1(6)*dck1(5,2)-c1(5)*dck1(6,2)
    v(435)=-(c1(5)*dck1(2,2))+c1(6)*dck1(4,2)-c1(2)*dck1(5,2)+c1(4)*dck1(6,2)
    v(434)=c1(3)*dck1(2,1)+c1(2)*dck1(3,1)+dck1(6,1)*v(437)
    v(433)=c1(4)*dck1(3,1)+c1(3)*dck1(4,1)-c1(6)*dck1(5,1)-c1(5)*dck1(6,1)
    v(432)=-(c1(5)*dck1(2,1))+c1(6)*dck1(4,1)-c1(2)*dck1(5,1)+c1(4)*dck1(6,1)
    v(431)=c0(2)*c0(3)-c0(6)**2
    v(430)=c0(3)*c0(4)-c0(5)*c0(6)
    v(429)=-(c0(2)*c0(5))+c0(4)*c0(6)
    v(428)=c1(2)*c1(3)-c1(6)**2
    v(427)=c1(3)*c1(4)-c1(5)*c1(6)
    v(426)=-(c1(2)*c1(5))+c1(4)*c1(6)
    v(261)=c1(5)*v(426)-c1(4)*v(427)+c1(1)*v(428)
    v(310)=1d0/v(261)**3
    v(569)=(-2d0)*v(310)
    v(271)=1d0/v(261)**2
    v(262)=c0(5)*v(429)-c0(4)*v(430)+c0(1)*v(431)
    v(452)=-(v(262)*v(271))
    v(260)=v(262)/v(261)
    v(329)=1d0/v(260)**0.16666666666666669d1
    v(648)=(-2d0/3d0)*v(329)
    v(270)=1d0/v(260)**0.6666666666666666d0
    v(556)=v(270)/3d0
    v(254)=v(260)**0.3333333333333333d0
    v(339)=dck1(5,1)*v(426)-dck1(4,1)*v(427)+dck1(1,1)*v(428)+c1(5)*v(432)-c1(4)*v(433)+c1(1)*v(434)
    v(347)=dck1(6,2)*v(437)+v(438)+v(439)
    v(351)=c1(1)*v(347)+dck1(5,2)*v(426)-dck1(4,2)*v(427)+dck1(1,2)*v(428)+c1(5)*v(435)-c1(4)*v(436)
    v(359)=dck1(6,3)*v(437)+v(442)+v(443)
    v(363)=c1(1)*v(359)+dck1(5,3)*v(426)-dck1(4,3)*v(427)+dck1(1,3)*v(428)+c1(5)*v(440)-c1(4)*v(441)
    v(338)=dck0(5,1)*v(429)-dck0(4,1)*v(430)+dck0(1,1)*v(431)+c0(5)*v(444)-c0(4)*v(445)+c0(1)*v(446)
    v(343)=v(338)/v(261)+v(339)*v(452)
    v(352)=dck0(6,2)*v(449)+v(450)+v(451)
    v(350)=c0(1)*v(352)+dck0(5,2)*v(429)-dck0(4,2)*v(430)+dck0(1,2)*v(431)+c0(5)*v(447)-c0(4)*v(448)
    v(355)=v(350)/v(261)+v(351)*v(452)
    v(364)=dck0(6,3)*v(449)+v(455)+v(456)
    v(362)=c0(1)*v(364)+dck0(5,3)*v(429)-dck0(4,3)*v(430)+dck0(1,3)*v(431)+c0(5)*v(453)-c0(4)*v(454)
    v(367)=v(362)/v(261)+v(363)*v(452)
    v(269)=v(343)*v(556)
    v(273)=v(355)*v(556)
    v(275)=v(367)*v(556)
    v(300)=dcl1(6,1)*v(437)+v(563)+v(564)
    v(301)=dcl1(6,2)*v(437)+v(565)+v(566)
    v(302)=dcl1(6,3)*v(437)+v(567)+v(568)
    v(303)=c1(1)*v(300)+dcl1(5,1)*v(426)-dcl1(4,1)*v(427)+dcl1(1,1)*v(428)+c1(5)*v(557)-c1(4)*v(560)
    v(304)=c1(1)*v(301)+dcl1(5,2)*v(426)-dcl1(4,2)*v(427)+dcl1(1,2)*v(428)+c1(5)*v(558)-c1(4)*v(561)
    v(305)=c1(1)*v(302)+dcl1(5,3)*v(426)-dcl1(4,3)*v(427)+dcl1(1,3)*v(428)+c1(5)*v(559)-c1(4)*v(562)
    v(306)=-(v(271)*v(303))
    v(307)=-(v(271)*v(304))
    v(308)=-(v(271)*v(305))
    v(309)=v(303)*v(569)
    v(604)=-(v(262)*v(309))
    v(311)=v(304)*v(569)
    v(671)=-(v(262)*v(311))
    v(312)=v(305)*v(569)
    v(737)=-(v(262)*v(312))
    v(319)=dcl0(6,1)*v(449)+v(576)+v(577)
    v(320)=dcl0(6,2)*v(449)+v(578)+v(579)
    v(321)=dcl0(6,3)*v(449)+v(580)+v(581)
    v(322)=c0(1)*v(319)+dcl0(5,1)*v(429)-dcl0(4,1)*v(430)+dcl0(1,1)*v(431)+c0(5)*v(570)-c0(4)*v(573)
    v(605)=-(v(271)*v(322))
    v(323)=c0(1)*v(320)+dcl0(5,2)*v(429)-dcl0(4,2)*v(430)+dcl0(1,2)*v(431)+c0(5)*v(571)-c0(4)*v(574)
    v(672)=-(v(271)*v(323))
    v(324)=c0(1)*v(321)+dcl0(5,3)*v(429)-dcl0(4,3)*v(430)+dcl0(1,3)*v(431)+c0(5)*v(572)-c0(4)*v(575)
    v(738)=-(v(271)*v(324))
    v(325)=v(262)*v(306)+v(322)/v(261)
    v(326)=v(262)*v(307)+v(323)/v(261)
    v(327)=v(262)*v(308)+v(324)/v(261)
    v(328)=v(325)*v(648)
    v(368)=(v(328)*v(367)+v(270)*(v(306)*v(362)+v(452)*(dck1(1,3)*v(300)+dcl1(1,1)*v(359)+dcl1(5,1)*v(440)-dcl1(4,1)*v(441)&
         &+v(428)*v(463)-v(427)*v(490)+v(426)*v(499)+dck1(5,3)*v(557)-dck1(4,3)*v(560)+c1(5)*(-(c1(5)*v(472))+c1(6)*v(490)-c1(2&
         &)*v(499)+c1(4)*v(508)+v(582)+v(583)+v(584)+v(585))-c1(4)*(c1(4)*v(481)+c1(3)*v(490)-c1(6)*v(499)-c1(5)*v(508)+v(586)+v&
         &(587)+v(588)+v(589))+c1(1)*(c1(3)*v(472)+c1(2)*v(481)+v(437)*v(508)+v(590)+v(591)+v(592)))+(dck0(1,3)*v(319)+dcl0(1,1&
         &)*v(364)+d2c0(1,3,1)*v(431)+dcl0(5,1)*v(453)-dcl0(4,1)*v(454)-v(430)*v(535)+v(429)*v(544)+dck0(5,3)*v(570)-dck0(4,3)*v&
         &(573)+c0(5)*(-(c0(5)*v(517))+c0(6)*v(535)-c0(2)*v(544)+c0(4)*v(553)+v(593)+v(594)+v(595)+v(596))-c0(4)*(c0(4)*v(526)+c0&
         &(3)*v(535)-c0(6)*v(544)-c0(5)*v(553)+v(597)+v(598)+v(599)+v(600))+c0(1)*(c0(3)*v(517)+c0(2)*v(526)+v(449)*v(553)+v(601)&
         &+v(602)+v(603)))/v(261)+v(363)*v(604)+v(363)*v(605)))/3d0
    v(356)=(v(328)*v(355)+v(270)*(v(306)*v(350)+v(351)*v(604)+v(351)*v(605)+v(452)*(dck1(1,2)*v(300)+dcl1(1,1)*v(347)+dcl1&
         &(5,1)*v(435)-dcl1(4,1)*v(436)+v(428)*v(460)-v(427)*v(487)+v(426)*v(496)+dck1(5,2)*v(557)-dck1(4,2)*v(560)+c1(5)*(-(c1(5&
         &)*v(469))+c1(6)*v(487)-c1(2)*v(496)+c1(4)*v(505)+v(608)+v(609)+v(610)+v(611))-c1(4)*(c1(4)*v(478)+c1(3)*v(487)-c1(6)*v&
         &(496)-c1(5)*v(505)+v(612)+v(613)+v(614)+v(615))+c1(1)*(c1(3)*v(469)+c1(2)*v(478)+v(437)*v(505)+dck1(6,2)*v(606)+v(616)&
         &+v(617)))+(dck0(1,2)*v(319)+dcl0(1,1)*v(352)+d2c0(1,2,1)*v(431)+dcl0(5,1)*v(447)-dcl0(4,1)*v(448)-v(430)*v(532)+v(429&
         &)*v(541)+dck0(5,2)*v(570)-dck0(4,2)*v(573)+c0(5)*(-(c0(5)*v(514))+c0(6)*v(532)-c0(2)*v(541)+c0(4)*v(550)+v(618)+v(619)&
         &+v(620)+v(621))-c0(4)*(c0(4)*v(523)+c0(3)*v(532)-c0(6)*v(541)-c0(5)*v(550)+v(622)+v(623)+v(624)+v(625))+c0(1)*(c0(3)*v&
         &(514)+c0(2)*v(523)+v(449)*v(550)+dck0(6,2)*v(607)+v(626)+v(627)))/v(261)))/3d0
    v(344)=(v(328)*v(343)+v(270)*(v(306)*v(338)+v(339)*v(604)+v(339)*v(605)+v(452)*(dck1(1,1)*v(300)+dcl1(5,1)*v(432)-dcl1&
         &(4,1)*v(433)+dcl1(1,1)*v(434)+v(428)*v(457)-v(427)*v(484)+v(426)*v(493)+dck1(5,1)*v(557)-dck1(4,1)*v(560)+c1(5)*(-(c1(5&
         &)*v(466))+c1(6)*v(484)-c1(2)*v(493)+c1(4)*v(502)+v(628)+v(629)+v(630)+v(631))-c1(4)*(c1(4)*v(475)+c1(3)*v(484)-c1(6)*v&
         &(493)-c1(5)*v(502)+v(632)+v(633)+v(634)+v(635))+c1(1)*(c1(3)*v(466)+c1(2)*v(475)+v(437)*v(502)+dck1(6,1)*v(606)+v(636)&
         &+v(637)))+(dck0(1,1)*v(319)+d2c0(1,1,1)*v(431)+dcl0(5,1)*v(444)-dcl0(4,1)*v(445)+dcl0(1,1)*v(446)-v(430)*v(529)+v(429&
         &)*v(538)+dck0(5,1)*v(570)-dck0(4,1)*v(573)+c0(5)*(-(c0(5)*v(511))+c0(6)*v(529)-c0(2)*v(538)+c0(4)*v(547)+v(638)+v(639)&
         &+v(640)+v(641))-c0(4)*(c0(4)*v(520)+c0(3)*v(529)-c0(6)*v(538)-c0(5)*v(547)+v(642)+v(643)+v(644)+v(645))+c0(1)*(c0(3)*v&
         &(511)+c0(2)*v(520)+v(449)*v(547)+dck0(6,1)*v(607)+v(646)+v(647)))/v(261)))/3d0
    v(330)=v(326)*v(648)
    v(369)=(v(330)*v(367)+v(270)*(v(307)*v(362)+v(452)*(dck1(1,3)*v(301)+dcl1(1,2)*v(359)+dcl1(5,2)*v(440)-dcl1(4,2)*v(441)&
         &+v(428)*v(464)-v(427)*v(491)+v(426)*v(500)+dck1(5,3)*v(558)-dck1(4,3)*v(561)+c1(5)*(-(c1(5)*v(473))+c1(6)*v(491)-c1(2&
         &)*v(500)+c1(4)*v(509)+v(649)+v(650)+v(651)+v(652))-c1(4)*(c1(4)*v(482)+c1(3)*v(491)-c1(6)*v(500)-c1(5)*v(509)+v(653)+v&
         &(654)+v(655)+v(656))+c1(1)*(c1(3)*v(473)+c1(2)*v(482)+v(437)*v(509)+v(657)+v(658)+v(659)))+(dck0(1,3)*v(320)+dcl0(1,2&
         &)*v(364)+d2c0(1,3,2)*v(431)+dcl0(5,2)*v(453)-dcl0(4,2)*v(454)-v(430)*v(536)+v(429)*v(545)+dck0(5,3)*v(571)-dck0(4,3)*v&
         &(574)+c0(5)*(-(c0(5)*v(518))+c0(6)*v(536)-c0(2)*v(545)+c0(4)*v(554)+v(660)+v(661)+v(662)+v(663))-c0(4)*(c0(4)*v(527)+c0&
         &(3)*v(536)-c0(6)*v(545)-c0(5)*v(554)+v(664)+v(665)+v(666)+v(667))+c0(1)*(c0(3)*v(518)+c0(2)*v(527)+v(449)*v(554)+v(668)&
         &+v(669)+v(670)))/v(261)+v(363)*v(671)+v(363)*v(672)))/3d0
    v(357)=(v(330)*v(355)+v(270)*(v(307)*v(350)+v(351)*v(671)+v(351)*v(672)+v(452)*(dck1(1,2)*v(301)+dcl1(1,2)*v(347)+dcl1&
         &(5,2)*v(435)-dcl1(4,2)*v(436)+v(428)*v(461)-v(427)*v(488)+v(426)*v(497)+dck1(5,2)*v(558)-dck1(4,2)*v(561)+c1(5)*(-(c1(5&
         &)*v(470))+c1(6)*v(488)-c1(2)*v(497)+c1(4)*v(506)+v(675)+v(676)+v(677)+v(678))-c1(4)*(c1(4)*v(479)+c1(3)*v(488)-c1(6)*v&
         &(497)-c1(5)*v(506)+v(679)+v(680)+v(681)+v(682))+c1(1)*(c1(3)*v(470)+c1(2)*v(479)+v(437)*v(506)+dck1(6,2)*v(673)+v(683)&
         &+v(684)))+(dck0(1,2)*v(320)+dcl0(1,2)*v(352)+d2c0(1,2,2)*v(431)+dcl0(5,2)*v(447)-dcl0(4,2)*v(448)-v(430)*v(533)+v(429&
         &)*v(542)+dck0(5,2)*v(571)-dck0(4,2)*v(574)+c0(5)*(-(c0(5)*v(515))+c0(6)*v(533)-c0(2)*v(542)+c0(4)*v(551)+v(685)+v(686)&
         &+v(687)+v(688))-c0(4)*(c0(4)*v(524)+c0(3)*v(533)-c0(6)*v(542)-c0(5)*v(551)+v(689)+v(690)+v(691)+v(692))+c0(1)*(c0(3)*v&
         &(515)+c0(2)*v(524)+v(449)*v(551)+dck0(6,2)*v(674)+v(693)+v(694)))/v(261)))/3d0
    v(345)=(v(330)*v(343)+v(270)*(v(307)*v(338)+v(339)*v(671)+v(339)*v(672)+v(452)*(dck1(1,1)*v(301)+dcl1(5,2)*v(432)-dcl1&
         &(4,2)*v(433)+dcl1(1,2)*v(434)+v(428)*v(458)-v(427)*v(485)+v(426)*v(494)+dck1(5,1)*v(558)-dck1(4,1)*v(561)+c1(5)*(-(c1(5&
         &)*v(467))+c1(6)*v(485)-c1(2)*v(494)+c1(4)*v(503)+v(695)+v(696)+v(697)+v(698))-c1(4)*(c1(4)*v(476)+c1(3)*v(485)-c1(6)*v&
         &(494)-c1(5)*v(503)+v(699)+v(700)+v(701)+v(702))+c1(1)*(c1(3)*v(467)+c1(2)*v(476)+v(437)*v(503)+dck1(6,1)*v(673)+v(703)&
         &+v(704)))+(dck0(1,1)*v(320)+d2c0(1,1,2)*v(431)+dcl0(5,2)*v(444)-dcl0(4,2)*v(445)+dcl0(1,2)*v(446)-v(430)*v(530)+v(429&
         &)*v(539)+dck0(5,1)*v(571)-dck0(4,1)*v(574)+c0(5)*(-(c0(5)*v(512))+c0(6)*v(530)-c0(2)*v(539)+c0(4)*v(548)+v(705)+v(706)&
         &+v(707)+v(708))-c0(4)*(c0(4)*v(521)+c0(3)*v(530)-c0(6)*v(539)-c0(5)*v(548)+v(709)+v(710)+v(711)+v(712))+c0(1)*(c0(3)*v&
         &(512)+c0(2)*v(521)+v(449)*v(548)+dck0(6,1)*v(674)+v(713)+v(714)))/v(261)))/3d0
    v(331)=v(327)*v(648)
    v(370)=(v(331)*v(367)+v(270)*(v(308)*v(362)+v(452)*(dck1(1,3)*v(302)+dcl1(1,3)*v(359)+dcl1(5,3)*v(440)-dcl1(4,3)*v(441)&
         &+v(428)*v(465)-v(427)*v(492)+v(426)*v(501)+dck1(5,3)*v(559)-dck1(4,3)*v(562)+c1(5)*(-(c1(5)*v(474))+c1(6)*v(492)-c1(2&
         &)*v(501)+c1(4)*v(510)+v(715)+v(716)+v(717)+v(718))-c1(4)*(c1(4)*v(483)+c1(3)*v(492)-c1(6)*v(501)-c1(5)*v(510)+v(719)+v&
         &(720)+v(721)+v(722))+c1(1)*(c1(3)*v(474)+c1(2)*v(483)+v(437)*v(510)+v(723)+v(724)+v(725)))+(dck0(1,3)*v(321)+dcl0(1,3&
         &)*v(364)+d2c0(1,3,3)*v(431)+dcl0(5,3)*v(453)-dcl0(4,3)*v(454)-v(430)*v(537)+v(429)*v(546)+dck0(5,3)*v(572)-dck0(4,3)*v&
         &(575)+c0(5)*(-(c0(5)*v(519))+c0(6)*v(537)-c0(2)*v(546)+c0(4)*v(555)+v(726)+v(727)+v(728)+v(729))-c0(4)*(c0(4)*v(528)+c0&
         &(3)*v(537)-c0(6)*v(546)-c0(5)*v(555)+v(730)+v(731)+v(732)+v(733))+c0(1)*(c0(3)*v(519)+c0(2)*v(528)+v(449)*v(555)+v(734)&
         &+v(735)+v(736)))/v(261)+v(363)*v(737)+v(363)*v(738)))/3d0
    v(358)=(v(331)*v(355)+v(270)*(v(308)*v(350)+v(351)*v(737)+v(351)*v(738)+v(452)*(dck1(1,2)*v(302)+dcl1(1,3)*v(347)+dcl1&
         &(5,3)*v(435)-dcl1(4,3)*v(436)+v(428)*v(462)-v(427)*v(489)+v(426)*v(498)+dck1(5,2)*v(559)-dck1(4,2)*v(562)+c1(5)*(-(c1(5&
         &)*v(471))+c1(6)*v(489)-c1(2)*v(498)+c1(4)*v(507)+v(741)+v(742)+v(743)+v(744))-c1(4)*(c1(4)*v(480)+c1(3)*v(489)-c1(6)*v&
         &(498)-c1(5)*v(507)+v(745)+v(746)+v(747)+v(748))+c1(1)*(c1(3)*v(471)+c1(2)*v(480)+v(437)*v(507)+dck1(6,2)*v(739)+v(749)&
         &+v(750)))+(dck0(1,2)*v(321)+dcl0(1,3)*v(352)+d2c0(1,2,3)*v(431)+dcl0(5,3)*v(447)-dcl0(4,3)*v(448)-v(430)*v(534)+v(429&
         &)*v(543)+dck0(5,2)*v(572)-dck0(4,2)*v(575)+c0(5)*(-(c0(5)*v(516))+c0(6)*v(534)-c0(2)*v(543)+c0(4)*v(552)+v(751)+v(752)&
         &+v(753)+v(754))-c0(4)*(c0(4)*v(525)+c0(3)*v(534)-c0(6)*v(543)-c0(5)*v(552)+v(755)+v(756)+v(757)+v(758))+c0(1)*(c0(3)*v&
         &(516)+c0(2)*v(525)+v(449)*v(552)+dck0(6,2)*v(740)+v(759)+v(760)))/v(261)))/3d0
    v(346)=(v(331)*v(343)+v(270)*(v(308)*v(338)+v(339)*v(737)+v(339)*v(738)+v(452)*(dck1(1,1)*v(302)+dcl1(5,3)*v(432)-dcl1&
         &(4,3)*v(433)+dcl1(1,3)*v(434)+v(428)*v(459)-v(427)*v(486)+v(426)*v(495)+dck1(5,1)*v(559)-dck1(4,1)*v(562)+c1(5)*(-(c1(5&
         &)*v(468))+c1(6)*v(486)-c1(2)*v(495)+c1(4)*v(504)+v(761)+v(762)+v(763)+v(764))-c1(4)*(c1(4)*v(477)+c1(3)*v(486)-c1(6)*v&
         &(495)-c1(5)*v(504)+v(765)+v(766)+v(767)+v(768))+c1(1)*(c1(3)*v(468)+c1(2)*v(477)+v(437)*v(504)+dck1(6,1)*v(739)+v(769)&
         &+v(770)))+(dck0(1,1)*v(321)+d2c0(1,1,3)*v(431)+dcl0(5,3)*v(444)-dcl0(4,3)*v(445)+dcl0(1,3)*v(446)-v(430)*v(531)+v(429&
         &)*v(540)+dck0(5,1)*v(572)-dck0(4,1)*v(575)+c0(5)*(-(c0(5)*v(513))+c0(6)*v(531)-c0(2)*v(540)+c0(4)*v(549)+v(771)+v(772)&
         &+v(773)+v(774))-c0(4)*(c0(4)*v(522)+c0(3)*v(531)-c0(6)*v(540)-c0(5)*v(549)+v(775)+v(776)+v(777)+v(778))+c0(1)*(c0(3)*v&
         &(513)+c0(2)*v(522)+v(449)*v(549)+dck0(6,1)*v(740)+v(779)+v(780)))/v(261)))/3d0
    v(332)=v(325)*v(556)
    v(333)=v(326)*v(556)
    v(334)=v(327)*v(556)
    d2c(1,1,1)=dcl1(1,1)*v(269)+dck1(1,1)*v(332)+c1(1)*v(344)+v(254)*v(457)
    d2c(1,1,2)=dcl1(1,2)*v(269)+dck1(1,1)*v(333)+c1(1)*v(345)+v(254)*v(458)
    d2c(1,1,3)=dcl1(1,3)*v(269)+dck1(1,1)*v(334)+c1(1)*v(346)+v(254)*v(459)
    d2c(1,2,1)=dcl1(1,1)*v(273)+dck1(1,2)*v(332)+c1(1)*v(356)+v(254)*v(460)
    d2c(1,2,2)=dcl1(1,2)*v(273)+dck1(1,2)*v(333)+c1(1)*v(357)+v(254)*v(461)
    d2c(1,2,3)=dcl1(1,3)*v(273)+dck1(1,2)*v(334)+c1(1)*v(358)+v(254)*v(462)
    d2c(1,3,1)=dcl1(1,1)*v(275)+dck1(1,3)*v(332)+c1(1)*v(368)+v(254)*v(463)
    d2c(1,3,2)=dcl1(1,2)*v(275)+dck1(1,3)*v(333)+c1(1)*v(369)+v(254)*v(464)
    d2c(1,3,3)=dcl1(1,3)*v(275)+dck1(1,3)*v(334)+c1(1)*v(370)+v(254)*v(465)
    d2c(2,1,1)=dcl1(2,1)*v(269)+dck1(2,1)*v(332)+c1(2)*v(344)+v(254)*v(466)
    d2c(2,1,2)=dcl1(2,2)*v(269)+dck1(2,1)*v(333)+c1(2)*v(345)+v(254)*v(467)
    d2c(2,1,3)=dcl1(2,3)*v(269)+dck1(2,1)*v(334)+c1(2)*v(346)+v(254)*v(468)
    d2c(2,2,1)=dcl1(2,1)*v(273)+dck1(2,2)*v(332)+c1(2)*v(356)+v(254)*v(469)
    d2c(2,2,2)=dcl1(2,2)*v(273)+dck1(2,2)*v(333)+c1(2)*v(357)+v(254)*v(470)
    d2c(2,2,3)=dcl1(2,3)*v(273)+dck1(2,2)*v(334)+c1(2)*v(358)+v(254)*v(471)
    d2c(2,3,1)=dcl1(2,1)*v(275)+dck1(2,3)*v(332)+c1(2)*v(368)+v(254)*v(472)
    d2c(2,3,2)=dcl1(2,2)*v(275)+dck1(2,3)*v(333)+c1(2)*v(369)+v(254)*v(473)
    d2c(2,3,3)=dcl1(2,3)*v(275)+dck1(2,3)*v(334)+c1(2)*v(370)+v(254)*v(474)
    d2c(3,1,1)=dcl1(3,1)*v(269)+dck1(3,1)*v(332)+c1(3)*v(344)+v(254)*v(475)
    d2c(3,1,2)=dcl1(3,2)*v(269)+dck1(3,1)*v(333)+c1(3)*v(345)+v(254)*v(476)
    d2c(3,1,3)=dcl1(3,3)*v(269)+dck1(3,1)*v(334)+c1(3)*v(346)+v(254)*v(477)
    d2c(3,2,1)=dcl1(3,1)*v(273)+dck1(3,2)*v(332)+c1(3)*v(356)+v(254)*v(478)
    d2c(3,2,2)=dcl1(3,2)*v(273)+dck1(3,2)*v(333)+c1(3)*v(357)+v(254)*v(479)
    d2c(3,2,3)=dcl1(3,3)*v(273)+dck1(3,2)*v(334)+c1(3)*v(358)+v(254)*v(480)
    d2c(3,3,1)=dcl1(3,1)*v(275)+dck1(3,3)*v(332)+c1(3)*v(368)+v(254)*v(481)
    d2c(3,3,2)=dcl1(3,2)*v(275)+dck1(3,3)*v(333)+c1(3)*v(369)+v(254)*v(482)
    d2c(3,3,3)=dcl1(3,3)*v(275)+dck1(3,3)*v(334)+c1(3)*v(370)+v(254)*v(483)
    d2c(4,1,1)=dcl1(4,1)*v(269)+dck1(4,1)*v(332)+c1(4)*v(344)+v(254)*v(484)
    d2c(4,1,2)=dcl1(4,2)*v(269)+dck1(4,1)*v(333)+c1(4)*v(345)+v(254)*v(485)
    d2c(4,1,3)=dcl1(4,3)*v(269)+dck1(4,1)*v(334)+c1(4)*v(346)+v(254)*v(486)
    d2c(4,2,1)=dcl1(4,1)*v(273)+dck1(4,2)*v(332)+c1(4)*v(356)+v(254)*v(487)
    d2c(4,2,2)=dcl1(4,2)*v(273)+dck1(4,2)*v(333)+c1(4)*v(357)+v(254)*v(488)
    d2c(4,2,3)=dcl1(4,3)*v(273)+dck1(4,2)*v(334)+c1(4)*v(358)+v(254)*v(489)
    d2c(4,3,1)=dcl1(4,1)*v(275)+dck1(4,3)*v(332)+c1(4)*v(368)+v(254)*v(490)
    d2c(4,3,2)=dcl1(4,2)*v(275)+dck1(4,3)*v(333)+c1(4)*v(369)+v(254)*v(491)
    d2c(4,3,3)=dcl1(4,3)*v(275)+dck1(4,3)*v(334)+c1(4)*v(370)+v(254)*v(492)
    d2c(5,1,1)=dcl1(5,1)*v(269)+dck1(5,1)*v(332)+c1(5)*v(344)+v(254)*v(493)
    d2c(5,1,2)=dcl1(5,2)*v(269)+dck1(5,1)*v(333)+c1(5)*v(345)+v(254)*v(494)
    d2c(5,1,3)=dcl1(5,3)*v(269)+dck1(5,1)*v(334)+c1(5)*v(346)+v(254)*v(495)
    d2c(5,2,1)=dcl1(5,1)*v(273)+dck1(5,2)*v(332)+c1(5)*v(356)+v(254)*v(496)
    d2c(5,2,2)=dcl1(5,2)*v(273)+dck1(5,2)*v(333)+c1(5)*v(357)+v(254)*v(497)
    d2c(5,2,3)=dcl1(5,3)*v(273)+dck1(5,2)*v(334)+c1(5)*v(358)+v(254)*v(498)
    d2c(5,3,1)=dcl1(5,1)*v(275)+dck1(5,3)*v(332)+c1(5)*v(368)+v(254)*v(499)
    d2c(5,3,2)=dcl1(5,2)*v(275)+dck1(5,3)*v(333)+c1(5)*v(369)+v(254)*v(500)
    d2c(5,3,3)=dcl1(5,3)*v(275)+dck1(5,3)*v(334)+c1(5)*v(370)+v(254)*v(501)
    d2c(6,1,1)=dcl1(6,1)*v(269)+dck1(6,1)*v(332)+c1(6)*v(344)+v(254)*v(502)
    d2c(6,1,2)=dcl1(6,2)*v(269)+dck1(6,1)*v(333)+c1(6)*v(345)+v(254)*v(503)
    d2c(6,1,3)=dcl1(6,3)*v(269)+dck1(6,1)*v(334)+c1(6)*v(346)+v(254)*v(504)
    d2c(6,2,1)=dcl1(6,1)*v(273)+dck1(6,2)*v(332)+c1(6)*v(356)+v(254)*v(505)
    d2c(6,2,2)=dcl1(6,2)*v(273)+dck1(6,2)*v(333)+c1(6)*v(357)+v(254)*v(506)
    d2c(6,2,3)=dcl1(6,3)*v(273)+dck1(6,2)*v(334)+c1(6)*v(358)+v(254)*v(507)
    d2c(6,3,1)=dcl1(6,1)*v(275)+dck1(6,3)*v(332)+c1(6)*v(368)+v(254)*v(508)
    d2c(6,3,2)=dcl1(6,2)*v(275)+dck1(6,3)*v(333)+c1(6)*v(369)+v(254)*v(509)
    d2c(6,3,3)=dcl1(6,3)*v(275)+dck1(6,3)*v(334)+c1(6)*v(370)+v(254)*v(510)
  END SUBROUTINE mls_combinedfbar2


!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           8 Feb 22 12:29:31  *
!**************************************************************
! User     : Full professional version
! Notebook : mls_combinedfbar2
! Evaluation time                 : 19 s    Mode  : Optimal
! Number of formulae              : 401     Method: Automatic
! Subroutine                      : mls_combinedfbar2 size: 12767
! Total size of Mathematica  code : 12767 subexpressions
! Total size of Fortran code      : 26108 bytes

!******************* S U B R O U T I N E **********************
  SUBROUTINE mls_combinedfbar2old(c,cb,dcleft,dcright,dcbleft,dcbright,dc2,dcb2,dcs2)
    IMPLICIT NONE
    DOUBLE PRECISION v(770),c(6),cb(6),dcleft(6,3),dcright(6,3),dcbleft(6,3),dcbright(6,3),dc2(6,3,3),dcb2(6,3,3)&
         &,dcs2(6,3,3)
    v(765)=dcbleft(2,1)*dcbright(3,3)
    v(764)=dcbleft(3,1)*dcbright(2,3)
    v(763)=-(dcbleft(5,1)*dcbright(6,3))
    v(762)=-(dcbleft(6,1)*dcbright(5,3))
    v(761)=dcbleft(3,1)*dcbright(4,3)
    v(760)=dcbleft(4,1)*dcbright(3,3)
    v(759)=dcbleft(4,1)*dcbright(6,3)
    v(758)=-(dcbleft(2,1)*dcbright(5,3))
    v(757)=dcbleft(6,1)*dcbright(4,3)
    v(756)=-(dcbleft(5,1)*dcbright(2,3))
    v(755)=dcleft(2,1)*dcright(3,3)
    v(754)=dcleft(3,1)*dcright(2,3)
    v(753)=-(dcleft(5,1)*dcright(6,3))
    v(752)=-(dcleft(6,1)*dcright(5,3))
    v(751)=dcleft(3,1)*dcright(4,3)
    v(750)=dcleft(4,1)*dcright(3,3)
    v(749)=dcleft(4,1)*dcright(6,3)
    v(748)=-(dcleft(2,1)*dcright(5,3))
    v(747)=dcleft(6,1)*dcright(4,3)
    v(746)=-(dcleft(5,1)*dcright(2,3))
    v(745)=dcbleft(2,2)*dcbright(3,3)
    v(744)=dcbleft(3,2)*dcbright(2,3)
    v(743)=-(dcbleft(5,2)*dcbright(6,3))
    v(742)=-(dcbleft(6,2)*dcbright(5,3))
    v(741)=dcbleft(3,2)*dcbright(4,3)
    v(740)=dcbleft(4,2)*dcbright(3,3)
    v(739)=dcbleft(4,2)*dcbright(6,3)
    v(738)=-(dcbleft(2,2)*dcbright(5,3))
    v(737)=dcbleft(6,2)*dcbright(4,3)
    v(736)=-(dcbleft(5,2)*dcbright(2,3))
    v(735)=dcleft(2,2)*dcright(3,3)
    v(734)=dcleft(3,2)*dcright(2,3)
    v(733)=-(dcleft(5,2)*dcright(6,3))
    v(732)=-(dcleft(6,2)*dcright(5,3))
    v(731)=dcleft(3,2)*dcright(4,3)
    v(730)=dcleft(4,2)*dcright(3,3)
    v(729)=dcleft(4,2)*dcright(6,3)
    v(728)=-(dcleft(2,2)*dcright(5,3))
    v(727)=dcleft(6,2)*dcright(4,3)
    v(726)=-(dcleft(5,2)*dcright(2,3))
    v(725)=(-2d0)*dcbright(6,3)
    v(724)=(-2d0)*dcright(6,3)
    v(722)=dcbleft(6,3)*v(725)
    v(721)=dcbleft(2,3)*dcbright(3,3)
    v(720)=dcbleft(3,3)*dcbright(2,3)
    v(719)=-(dcbleft(5,3)*dcbright(6,3))
    v(718)=-(dcbleft(6,3)*dcbright(5,3))
    v(717)=dcbleft(3,3)*dcbright(4,3)
    v(716)=dcbleft(4,3)*dcbright(3,3)
    v(715)=dcbleft(4,3)*dcbright(6,3)
    v(714)=-(dcbleft(2,3)*dcbright(5,3))
    v(713)=dcbleft(6,3)*dcbright(4,3)
    v(712)=-(dcbleft(5,3)*dcbright(2,3))
    v(711)=dcleft(6,3)*v(724)
    v(710)=dcleft(2,3)*dcright(3,3)
    v(709)=dcleft(3,3)*dcright(2,3)
    v(708)=-(dcleft(5,3)*dcright(6,3))
    v(707)=-(dcleft(6,3)*dcright(5,3))
    v(706)=dcleft(3,3)*dcright(4,3)
    v(705)=dcleft(4,3)*dcright(3,3)
    v(704)=dcleft(4,3)*dcright(6,3)
    v(703)=-(dcleft(2,3)*dcright(5,3))
    v(702)=dcleft(6,3)*dcright(4,3)
    v(701)=-(dcleft(5,3)*dcright(2,3))
    v(700)=dcbleft(2,1)*dcbright(3,2)
    v(699)=dcbleft(3,1)*dcbright(2,2)
    v(698)=-(dcbleft(5,1)*dcbright(6,2))
    v(697)=-(dcbleft(6,1)*dcbright(5,2))
    v(696)=dcbleft(3,1)*dcbright(4,2)
    v(695)=dcbleft(4,1)*dcbright(3,2)
    v(694)=dcbleft(4,1)*dcbright(6,2)
    v(693)=-(dcbleft(2,1)*dcbright(5,2))
    v(692)=dcbleft(6,1)*dcbright(4,2)
    v(691)=-(dcbleft(5,1)*dcbright(2,2))
    v(690)=dcleft(2,1)*dcright(3,2)
    v(689)=dcleft(3,1)*dcright(2,2)
    v(688)=-(dcleft(5,1)*dcright(6,2))
    v(687)=-(dcleft(6,1)*dcright(5,2))
    v(686)=dcleft(3,1)*dcright(4,2)
    v(685)=dcleft(4,1)*dcright(3,2)
    v(684)=dcleft(4,1)*dcright(6,2)
    v(683)=-(dcleft(2,1)*dcright(5,2))
    v(682)=dcleft(6,1)*dcright(4,2)
    v(681)=-(dcleft(5,1)*dcright(2,2))
    v(680)=dcbleft(2,2)*dcbright(3,2)
    v(679)=dcbleft(3,2)*dcbright(2,2)
    v(678)=-(dcbleft(5,2)*dcbright(6,2))
    v(677)=-(dcbleft(6,2)*dcbright(5,2))
    v(676)=dcbleft(3,2)*dcbright(4,2)
    v(675)=dcbleft(4,2)*dcbright(3,2)
    v(674)=dcbleft(4,2)*dcbright(6,2)
    v(673)=-(dcbleft(2,2)*dcbright(5,2))
    v(672)=dcbleft(6,2)*dcbright(4,2)
    v(671)=-(dcbleft(5,2)*dcbright(2,2))
    v(670)=dcleft(2,2)*dcright(3,2)
    v(669)=dcleft(3,2)*dcright(2,2)
    v(668)=-(dcleft(5,2)*dcright(6,2))
    v(667)=-(dcleft(6,2)*dcright(5,2))
    v(666)=dcleft(3,2)*dcright(4,2)
    v(665)=dcleft(4,2)*dcright(3,2)
    v(664)=dcleft(4,2)*dcright(6,2)
    v(663)=-(dcleft(2,2)*dcright(5,2))
    v(662)=dcleft(6,2)*dcright(4,2)
    v(661)=-(dcleft(5,2)*dcright(2,2))
    v(660)=(-2d0)*dcbright(6,2)
    v(659)=(-2d0)*dcright(6,2)
    v(657)=dcbleft(6,3)*v(660)
    v(656)=dcbleft(2,3)*dcbright(3,2)
    v(655)=dcbleft(3,3)*dcbright(2,2)
    v(654)=-(dcbleft(5,3)*dcbright(6,2))
    v(653)=-(dcbleft(6,3)*dcbright(5,2))
    v(652)=dcbleft(3,3)*dcbright(4,2)
    v(651)=dcbleft(4,3)*dcbright(3,2)
    v(650)=dcbleft(4,3)*dcbright(6,2)
    v(649)=-(dcbleft(2,3)*dcbright(5,2))
    v(648)=dcbleft(6,3)*dcbright(4,2)
    v(647)=-(dcbleft(5,3)*dcbright(2,2))
    v(646)=dcleft(6,3)*v(659)
    v(645)=dcleft(2,3)*dcright(3,2)
    v(644)=dcleft(3,3)*dcright(2,2)
    v(643)=-(dcleft(5,3)*dcright(6,2))
    v(642)=-(dcleft(6,3)*dcright(5,2))
    v(641)=dcleft(3,3)*dcright(4,2)
    v(640)=dcleft(4,3)*dcright(3,2)
    v(639)=dcleft(4,3)*dcright(6,2)
    v(638)=-(dcleft(2,3)*dcright(5,2))
    v(637)=dcleft(6,3)*dcright(4,2)
    v(636)=-(dcleft(5,3)*dcright(2,2))
    v(634)=dcbleft(2,1)*dcbright(3,1)
    v(633)=dcbleft(3,1)*dcbright(2,1)
    v(632)=-(dcbleft(5,1)*dcbright(6,1))
    v(631)=-(dcbleft(6,1)*dcbright(5,1))
    v(630)=dcbleft(3,1)*dcbright(4,1)
    v(629)=dcbleft(4,1)*dcbright(3,1)
    v(628)=dcbleft(4,1)*dcbright(6,1)
    v(627)=-(dcbleft(2,1)*dcbright(5,1))
    v(626)=dcbleft(6,1)*dcbright(4,1)
    v(625)=-(dcbleft(5,1)*dcbright(2,1))
    v(624)=dcleft(2,1)*dcright(3,1)
    v(623)=dcleft(3,1)*dcright(2,1)
    v(622)=-(dcleft(5,1)*dcright(6,1))
    v(621)=-(dcleft(6,1)*dcright(5,1))
    v(620)=dcleft(3,1)*dcright(4,1)
    v(619)=dcleft(4,1)*dcright(3,1)
    v(618)=dcleft(4,1)*dcright(6,1)
    v(617)=-(dcleft(2,1)*dcright(5,1))
    v(616)=dcleft(6,1)*dcright(4,1)
    v(615)=-(dcleft(5,1)*dcright(2,1))
    v(614)=dcbleft(2,2)*dcbright(3,1)
    v(613)=dcbleft(3,2)*dcbright(2,1)
    v(612)=-(dcbleft(5,2)*dcbright(6,1))
    v(611)=-(dcbleft(6,2)*dcbright(5,1))
    v(610)=dcbleft(3,2)*dcbright(4,1)
    v(609)=dcbleft(4,2)*dcbright(3,1)
    v(608)=dcbleft(4,2)*dcbright(6,1)
    v(607)=-(dcbleft(2,2)*dcbright(5,1))
    v(606)=dcbleft(6,2)*dcbright(4,1)
    v(605)=-(dcbleft(5,2)*dcbright(2,1))
    v(604)=dcleft(2,2)*dcright(3,1)
    v(603)=dcleft(3,2)*dcright(2,1)
    v(602)=-(dcleft(5,2)*dcright(6,1))
    v(601)=-(dcleft(6,2)*dcright(5,1))
    v(600)=dcleft(3,2)*dcright(4,1)
    v(599)=dcleft(4,2)*dcright(3,1)
    v(598)=dcleft(4,2)*dcright(6,1)
    v(597)=-(dcleft(2,2)*dcright(5,1))
    v(596)=dcleft(6,2)*dcright(4,1)
    v(595)=-(dcleft(5,2)*dcright(2,1))
    v(594)=(-2d0)*dcbright(6,1)
    v(593)=(-2d0)*dcright(6,1)
    v(591)=dcbleft(6,3)*v(594)
    v(590)=dcbleft(2,3)*dcbright(3,1)
    v(589)=dcbleft(3,3)*dcbright(2,1)
    v(588)=-(dcbleft(5,3)*dcbright(6,1))
    v(587)=-(dcbleft(6,3)*dcbright(5,1))
    v(586)=dcbleft(3,3)*dcbright(4,1)
    v(585)=dcbleft(4,3)*dcbright(3,1)
    v(584)=dcbleft(4,3)*dcbright(6,1)
    v(583)=-(dcbleft(2,3)*dcbright(5,1))
    v(582)=dcbleft(6,3)*dcbright(4,1)
    v(581)=-(dcbleft(5,3)*dcbright(2,1))
    v(580)=dcleft(6,3)*v(593)
    v(579)=dcleft(2,3)*dcright(3,1)
    v(578)=dcleft(3,3)*dcright(2,1)
    v(577)=-(dcleft(5,3)*dcright(6,1))
    v(576)=-(dcleft(6,3)*dcright(5,1))
    v(575)=dcleft(3,3)*dcright(4,1)
    v(574)=dcleft(4,3)*dcright(3,1)
    v(573)=dcleft(4,3)*dcright(6,1)
    v(572)=-(dcleft(2,3)*dcright(5,1))
    v(571)=dcleft(6,3)*dcright(4,1)
    v(570)=-(dcleft(5,3)*dcright(2,1))
    v(569)=cb(2)*dcbright(3,3)
    v(568)=cb(3)*dcbright(2,3)
    v(567)=cb(2)*dcbright(3,2)
    v(566)=cb(3)*dcbright(2,2)
    v(565)=cb(2)*dcbright(3,1)
    v(564)=cb(3)*dcbright(2,1)
    v(563)=cb(4)*dcbright(3,3)+cb(3)*dcbright(4,3)-cb(6)*dcbright(5,3)-cb(5)*dcbright(6,3)
    v(562)=cb(4)*dcbright(3,2)+cb(3)*dcbright(4,2)-cb(6)*dcbright(5,2)-cb(5)*dcbright(6,2)
    v(561)=cb(4)*dcbright(3,1)+cb(3)*dcbright(4,1)-cb(6)*dcbright(5,1)-cb(5)*dcbright(6,1)
    v(560)=-(cb(5)*dcbright(2,3))+cb(6)*dcbright(4,3)-cb(2)*dcbright(5,3)+cb(4)*dcbright(6,3)
    v(559)=-(cb(5)*dcbright(2,2))+cb(6)*dcbright(4,2)-cb(2)*dcbright(5,2)+cb(4)*dcbright(6,2)
    v(558)=-(cb(5)*dcbright(2,1))+cb(6)*dcbright(4,1)-cb(2)*dcbright(5,1)+cb(4)*dcbright(6,1)
    v(557)=c(2)*dcright(3,3)
    v(556)=c(3)*dcright(2,3)
    v(555)=c(2)*dcright(3,2)
    v(554)=c(3)*dcright(2,2)
    v(553)=c(2)*dcright(3,1)
    v(552)=c(3)*dcright(2,1)
    v(551)=c(4)*dcright(3,3)+c(3)*dcright(4,3)-c(6)*dcright(5,3)-c(5)*dcright(6,3)
    v(550)=c(4)*dcright(3,2)+c(3)*dcright(4,2)-c(6)*dcright(5,2)-c(5)*dcright(6,2)
    v(549)=c(4)*dcright(3,1)+c(3)*dcright(4,1)-c(6)*dcright(5,1)-c(5)*dcright(6,1)
    v(548)=-(c(5)*dcright(2,3))+c(6)*dcright(4,3)-c(2)*dcright(5,3)+c(4)*dcright(6,3)
    v(547)=-(c(5)*dcright(2,2))+c(6)*dcright(4,2)-c(2)*dcright(5,2)+c(4)*dcright(6,2)
    v(546)=-(c(5)*dcright(2,1))+c(6)*dcright(4,1)-c(2)*dcright(5,1)+c(4)*dcright(6,1)
    v(544)=dcb2(6,3,3)
    v(543)=dcb2(6,3,2)
    v(542)=dcb2(6,3,1)
    v(541)=dcb2(6,2,3)
    v(540)=dcb2(6,2,2)
    v(539)=dcb2(6,2,1)
    v(538)=dcb2(6,1,3)
    v(537)=dcb2(6,1,2)
    v(536)=dcb2(6,1,1)
    v(535)=dcb2(5,3,3)
    v(534)=dcb2(5,3,2)
    v(533)=dcb2(5,3,1)
    v(532)=dcb2(5,2,3)
    v(531)=dcb2(5,2,2)
    v(530)=dcb2(5,2,1)
    v(529)=dcb2(5,1,3)
    v(528)=dcb2(5,1,2)
    v(527)=dcb2(5,1,1)
    v(526)=dcb2(4,3,3)
    v(525)=dcb2(4,3,2)
    v(524)=dcb2(4,3,1)
    v(523)=dcb2(4,2,3)
    v(522)=dcb2(4,2,2)
    v(521)=dcb2(4,2,1)
    v(520)=dcb2(4,1,3)
    v(519)=dcb2(4,1,2)
    v(518)=dcb2(4,1,1)
    v(517)=dcb2(3,3,3)
    v(516)=dcb2(3,3,2)
    v(515)=dcb2(3,3,1)
    v(514)=dcb2(3,2,3)
    v(513)=dcb2(3,2,2)
    v(512)=dcb2(3,2,1)
    v(511)=dcb2(3,1,3)
    v(510)=dcb2(3,1,2)
    v(509)=dcb2(3,1,1)
    v(508)=dcb2(2,3,3)
    v(507)=dcb2(2,3,2)
    v(506)=dcb2(2,3,1)
    v(505)=dcb2(2,2,3)
    v(504)=dcb2(2,2,2)
    v(503)=dcb2(2,2,1)
    v(502)=dcb2(2,1,3)
    v(501)=dcb2(2,1,2)
    v(500)=dcb2(2,1,1)
    v(499)=dc2(6,3,3)
    v(498)=dc2(6,3,2)
    v(497)=dc2(6,3,1)
    v(496)=dc2(6,2,3)
    v(495)=dc2(6,2,2)
    v(494)=dc2(6,2,1)
    v(493)=dc2(6,1,3)
    v(492)=dc2(6,1,2)
    v(491)=dc2(6,1,1)
    v(490)=dc2(5,3,3)
    v(489)=dc2(5,3,2)
    v(488)=dc2(5,3,1)
    v(487)=dc2(5,2,3)
    v(486)=dc2(5,2,2)
    v(485)=dc2(5,2,1)
    v(484)=dc2(5,1,3)
    v(483)=dc2(5,1,2)
    v(482)=dc2(5,1,1)
    v(481)=dc2(4,3,3)
    v(480)=dc2(4,3,2)
    v(479)=dc2(4,3,1)
    v(478)=dc2(4,2,3)
    v(477)=dc2(4,2,2)
    v(476)=dc2(4,2,1)
    v(475)=dc2(4,1,3)
    v(474)=dc2(4,1,2)
    v(473)=dc2(4,1,1)
    v(472)=dc2(3,3,3)
    v(471)=dc2(3,3,2)
    v(470)=dc2(3,3,1)
    v(469)=dc2(3,2,3)
    v(468)=dc2(3,2,2)
    v(467)=dc2(3,2,1)
    v(466)=dc2(3,1,3)
    v(465)=dc2(3,1,2)
    v(464)=dc2(3,1,1)
    v(463)=dc2(2,3,3)
    v(462)=dc2(2,3,2)
    v(461)=dc2(2,3,1)
    v(460)=dc2(2,2,3)
    v(459)=dc2(2,2,2)
    v(458)=dc2(2,2,1)
    v(457)=dc2(2,1,3)
    v(456)=dc2(2,1,2)
    v(455)=dc2(2,1,1)
    v(454)=dc2(1,3,3)
    v(453)=dc2(1,3,2)
    v(452)=dc2(1,3,1)
    v(451)=dc2(1,2,3)
    v(450)=dc2(1,2,2)
    v(449)=dc2(1,2,1)
    v(448)=dc2(1,1,3)
    v(447)=dc2(1,1,2)
    v(446)=dc2(1,1,1)
    v(445)=cb(2)*dcbleft(3,3)
    v(444)=cb(3)*dcbleft(2,3)
    v(443)=cb(4)*dcbleft(3,3)+cb(3)*dcbleft(4,3)-cb(6)*dcbleft(5,3)-cb(5)*dcbleft(6,3)
    v(442)=-(cb(5)*dcbleft(2,3))+cb(6)*dcbleft(4,3)-cb(2)*dcbleft(5,3)+cb(4)*dcbleft(6,3)
    v(441)=cb(2)*dcbleft(3,2)
    v(440)=cb(3)*dcbleft(2,2)
    v(439)=(-2d0)*cb(6)
    v(438)=cb(4)*dcbleft(3,2)+cb(3)*dcbleft(4,2)-cb(6)*dcbleft(5,2)-cb(5)*dcbleft(6,2)
    v(437)=-(cb(5)*dcbleft(2,2))+cb(6)*dcbleft(4,2)-cb(2)*dcbleft(5,2)+cb(4)*dcbleft(6,2)
    v(436)=cb(3)*dcbleft(2,1)+cb(2)*dcbleft(3,1)+dcbleft(6,1)*v(439)
    v(435)=cb(4)*dcbleft(3,1)+cb(3)*dcbleft(4,1)-cb(6)*dcbleft(5,1)-cb(5)*dcbleft(6,1)
    v(434)=-(cb(5)*dcbleft(2,1))+cb(6)*dcbleft(4,1)-cb(2)*dcbleft(5,1)+cb(4)*dcbleft(6,1)
    v(433)=c(2)*dcleft(3,3)
    v(432)=c(3)*dcleft(2,3)
    v(431)=c(4)*dcleft(3,3)+c(3)*dcleft(4,3)-c(6)*dcleft(5,3)-c(5)*dcleft(6,3)
    v(430)=-(c(5)*dcleft(2,3))+c(6)*dcleft(4,3)-c(2)*dcleft(5,3)+c(4)*dcleft(6,3)
    v(429)=c(2)*dcleft(3,2)
    v(428)=c(3)*dcleft(2,2)
    v(427)=(-2d0)*c(6)
    v(426)=c(4)*dcleft(3,2)+c(3)*dcleft(4,2)-c(6)*dcleft(5,2)-c(5)*dcleft(6,2)
    v(425)=-(c(5)*dcleft(2,2))+c(6)*dcleft(4,2)-c(2)*dcleft(5,2)+c(4)*dcleft(6,2)
    v(424)=c(3)*dcleft(2,1)+c(2)*dcleft(3,1)+dcleft(6,1)*v(427)
    v(423)=c(4)*dcleft(3,1)+c(3)*dcleft(4,1)-c(6)*dcleft(5,1)-c(5)*dcleft(6,1)
    v(422)=-(c(5)*dcleft(2,1))+c(6)*dcleft(4,1)-c(2)*dcleft(5,1)+c(4)*dcleft(6,1)
    v(420)=cb(2)*cb(3)-cb(6)**2
    v(419)=cb(3)*cb(4)-cb(5)*cb(6)
    v(418)=-(cb(2)*cb(5))+cb(4)*cb(6)
    v(417)=c(2)*c(3)-c(6)**2
    v(416)=c(3)*c(4)-c(5)*c(6)
    v(415)=-(c(2)*c(5))+c(4)*c(6)
    v(260)=c(5)*v(415)-c(4)*v(416)+c(1)*v(417)
    v(545)=1d0/(3d0*v(260))
    v(304)=1d0/v(260)**2
    v(635)=(-1d0/3d0)*v(304)
    v(261)=cb(5)*v(418)-cb(4)*v(419)+cb(1)*v(420)
    v(421)=v(261)/3d0
    v(321)=((2d0/3d0)*v(261))/v(260)**3
    v(269)=-(v(304)*v(421))
    v(254)=v(421)/v(260)
    v(334)=dcleft(5,1)*v(415)-dcleft(4,1)*v(416)+dcleft(1,1)*v(417)+c(5)*v(422)-c(4)*v(423)+c(1)*v(424)
    v(338)=dcleft(6,2)*v(427)+v(428)+v(429)
    v(345)=c(1)*v(338)+dcleft(5,2)*v(415)-dcleft(4,2)*v(416)+dcleft(1,2)*v(417)+c(5)*v(425)-c(4)*v(426)
    v(349)=dcleft(6,3)*v(427)+v(432)+v(433)
    v(356)=c(1)*v(349)+dcleft(5,3)*v(415)-dcleft(4,3)*v(416)+dcleft(1,3)*v(417)+c(5)*v(430)-c(4)*v(431)
    v(330)=dcbleft(5,1)*v(418)-dcbleft(4,1)*v(419)+dcbleft(1,1)*v(420)+cb(5)*v(434)-cb(4)*v(435)+cb(1)*v(436)
    v(342)=dcbleft(6,2)*v(439)+v(440)+v(441)
    v(341)=cb(1)*v(342)+dcbleft(5,2)*v(418)-dcbleft(4,2)*v(419)+dcbleft(1,2)*v(420)+cb(5)*v(437)-cb(4)*v(438)
    v(353)=dcbleft(6,3)*v(439)+v(444)+v(445)
    v(352)=cb(1)*v(353)+dcbleft(5,3)*v(418)-dcbleft(4,3)*v(419)+dcbleft(1,3)*v(420)+cb(5)*v(442)-cb(4)*v(443)
    v(268)=v(269)*v(334)+v(330)*v(545)
    v(271)=v(269)*v(345)+v(341)*v(545)
    v(273)=v(269)*v(356)+v(352)*v(545)
    v(298)=dcright(6,1)*v(427)+v(552)+v(553)
    v(299)=dcright(6,2)*v(427)+v(554)+v(555)
    v(300)=dcright(6,3)*v(427)+v(556)+v(557)
    v(301)=c(1)*v(298)+dcright(5,1)*v(415)-dcright(4,1)*v(416)+dcright(1,1)*v(417)+c(5)*v(546)-c(4)*v(549)
    v(302)=c(1)*v(299)+dcright(5,2)*v(415)-dcright(4,2)*v(416)+dcright(1,2)*v(417)+c(5)*v(547)-c(4)*v(550)
    v(303)=c(1)*v(300)+dcright(5,3)*v(415)-dcright(4,3)*v(416)+dcright(1,3)*v(417)+c(5)*v(548)-c(4)*v(551)
    v(305)=-(v(301)*v(304))
    v(592)=v(305)/3d0
    v(306)=-(v(302)*v(304))
    v(658)=v(306)/3d0
    v(307)=-(v(303)*v(304))
    v(723)=v(307)/3d0
    v(314)=dcbright(6,1)*v(439)+v(564)+v(565)
    v(315)=dcbright(6,2)*v(439)+v(566)+v(567)
    v(316)=dcbright(6,3)*v(439)+v(568)+v(569)
    v(317)=cb(1)*v(314)+dcbright(5,1)*v(418)-dcbright(4,1)*v(419)+dcbright(1,1)*v(420)+cb(5)*v(558)-cb(4)*v(561)
    v(318)=cb(1)*v(315)+dcbright(5,2)*v(418)-dcbright(4,2)*v(419)+dcbright(1,2)*v(420)+cb(5)*v(559)-cb(4)*v(562)
    v(319)=cb(1)*v(316)+dcbright(5,3)*v(418)-dcbright(4,3)*v(419)+dcbright(1,3)*v(420)+cb(5)*v(560)-cb(4)*v(563)
    v(320)=v(301)*v(321)+v(317)*v(635)
    v(357)=v(320)*v(356)+v(269)*(dcleft(1,3)*v(298)+dcright(1,1)*v(349)+dcright(5,1)*v(430)-dcright(4,1)*v(431)+v(417)*v&
         &(452)-v(416)*v(479)+v(415)*v(488)+dcleft(5,3)*v(546)-dcleft(4,3)*v(549)+c(5)*(-(c(5)*v(461))+c(6)*v(479)-c(2)*v(488)+c&
         &(4)*v(497)+v(570)+v(571)+v(572)+v(573))-c(4)*(c(4)*v(470)+c(3)*v(479)-c(6)*v(488)-c(5)*v(497)+v(574)+v(575)+v(576)+v&
         &(577))+c(1)*(c(3)*v(461)+c(2)*v(470)+v(427)*v(497)+v(578)+v(579)+v(580)))+v(545)*(dcbleft(1,3)*v(314)+dcbright(1,1)*v&
         &(353)+dcb2(1,3,1)*v(420)+dcbright(5,1)*v(442)-dcbright(4,1)*v(443)-v(419)*v(524)+v(418)*v(533)+dcbleft(5,3)*v(558)&
         &-dcbleft(4,3)*v(561)+cb(5)*(-(cb(5)*v(506))+cb(6)*v(524)-cb(2)*v(533)+cb(4)*v(542)+v(581)+v(582)+v(583)+v(584))-cb(4)*&
         &(cb(4)*v(515)+cb(3)*v(524)-cb(6)*v(533)-cb(5)*v(542)+v(585)+v(586)+v(587)+v(588))+cb(1)*(cb(3)*v(506)+cb(2)*v(515)+v&
         &(439)*v(542)+v(589)+v(590)+v(591)))+v(352)*v(592)
    v(346)=v(320)*v(345)+v(341)*v(592)+v(269)*(dcleft(1,2)*v(298)+dcright(1,1)*v(338)+dcright(5,1)*v(425)-dcright(4,1)*v&
         &(426)+v(417)*v(449)-v(416)*v(476)+v(415)*v(485)+dcleft(5,2)*v(546)-dcleft(4,2)*v(549)+c(5)*(-(c(5)*v(458))+c(6)*v(476)&
         &-c(2)*v(485)+c(4)*v(494)+v(595)+v(596)+v(597)+v(598))-c(4)*(c(4)*v(467)+c(3)*v(476)-c(6)*v(485)-c(5)*v(494)+v(599)+v&
         &(600)+v(601)+v(602))+c(1)*(c(3)*v(458)+c(2)*v(467)+v(427)*v(494)+dcleft(6,2)*v(593)+v(603)+v(604)))+v(545)*(dcbleft(1,2&
         &)*v(314)+dcbright(1,1)*v(342)+dcb2(1,2,1)*v(420)+dcbright(5,1)*v(437)-dcbright(4,1)*v(438)-v(419)*v(521)+v(418)*v(530)&
         &+dcbleft(5,2)*v(558)-dcbleft(4,2)*v(561)+cb(5)*(-(cb(5)*v(503))+cb(6)*v(521)-cb(2)*v(530)+cb(4)*v(539)+v(605)+v(606)+v&
         &(607)+v(608))-cb(4)*(cb(4)*v(512)+cb(3)*v(521)-cb(6)*v(530)-cb(5)*v(539)+v(609)+v(610)+v(611)+v(612))+cb(1)*(cb(3)*v&
         &(503)+cb(2)*v(512)+v(439)*v(539)+dcbleft(6,2)*v(594)+v(613)+v(614)))
    v(335)=v(320)*v(334)+v(330)*v(592)+v(269)*(dcleft(1,1)*v(298)+dcright(5,1)*v(422)-dcright(4,1)*v(423)+dcright(1,1)*v&
         &(424)+v(417)*v(446)-v(416)*v(473)+v(415)*v(482)+dcleft(5,1)*v(546)-dcleft(4,1)*v(549)+c(5)*(-(c(5)*v(455))+c(6)*v(473)&
         &-c(2)*v(482)+c(4)*v(491)+v(615)+v(616)+v(617)+v(618))-c(4)*(c(4)*v(464)+c(3)*v(473)-c(6)*v(482)-c(5)*v(491)+v(619)+v&
         &(620)+v(621)+v(622))+c(1)*(c(3)*v(455)+c(2)*v(464)+v(427)*v(491)+dcleft(6,1)*v(593)+v(623)+v(624)))+v(545)*(dcbleft(1,1&
         &)*v(314)+dcb2(1,1,1)*v(420)+dcbright(5,1)*v(434)-dcbright(4,1)*v(435)+dcbright(1,1)*v(436)-v(419)*v(518)+v(418)*v(527)&
         &+dcbleft(5,1)*v(558)-dcbleft(4,1)*v(561)+cb(5)*(-(cb(5)*v(500))+cb(6)*v(518)-cb(2)*v(527)+cb(4)*v(536)+v(625)+v(626)+v&
         &(627)+v(628))-cb(4)*(cb(4)*v(509)+cb(3)*v(518)-cb(6)*v(527)-cb(5)*v(536)+v(629)+v(630)+v(631)+v(632))+cb(1)*(cb(3)*v&
         &(500)+cb(2)*v(509)+v(439)*v(536)+dcbleft(6,1)*v(594)+v(633)+v(634)))
    v(322)=v(302)*v(321)+v(318)*v(635)
    v(358)=v(322)*v(356)+v(269)*(dcleft(1,3)*v(299)+dcright(1,2)*v(349)+dcright(5,2)*v(430)-dcright(4,2)*v(431)+v(417)*v&
         &(453)-v(416)*v(480)+v(415)*v(489)+dcleft(5,3)*v(547)-dcleft(4,3)*v(550)+c(5)*(-(c(5)*v(462))+c(6)*v(480)-c(2)*v(489)+c&
         &(4)*v(498)+v(636)+v(637)+v(638)+v(639))-c(4)*(c(4)*v(471)+c(3)*v(480)-c(6)*v(489)-c(5)*v(498)+v(640)+v(641)+v(642)+v&
         &(643))+c(1)*(c(3)*v(462)+c(2)*v(471)+v(427)*v(498)+v(644)+v(645)+v(646)))+v(545)*(dcbleft(1,3)*v(315)+dcbright(1,2)*v&
         &(353)+dcb2(1,3,2)*v(420)+dcbright(5,2)*v(442)-dcbright(4,2)*v(443)-v(419)*v(525)+v(418)*v(534)+dcbleft(5,3)*v(559)&
         &-dcbleft(4,3)*v(562)+cb(5)*(-(cb(5)*v(507))+cb(6)*v(525)-cb(2)*v(534)+cb(4)*v(543)+v(647)+v(648)+v(649)+v(650))-cb(4)*&
         &(cb(4)*v(516)+cb(3)*v(525)-cb(6)*v(534)-cb(5)*v(543)+v(651)+v(652)+v(653)+v(654))+cb(1)*(cb(3)*v(507)+cb(2)*v(516)+v&
         &(439)*v(543)+v(655)+v(656)+v(657)))+v(352)*v(658)
    v(347)=v(322)*v(345)+v(341)*v(658)+v(269)*(dcleft(1,2)*v(299)+dcright(1,2)*v(338)+dcright(5,2)*v(425)-dcright(4,2)*v&
         &(426)+v(417)*v(450)-v(416)*v(477)+v(415)*v(486)+dcleft(5,2)*v(547)-dcleft(4,2)*v(550)+c(5)*(-(c(5)*v(459))+c(6)*v(477)&
         &-c(2)*v(486)+c(4)*v(495)+v(661)+v(662)+v(663)+v(664))-c(4)*(c(4)*v(468)+c(3)*v(477)-c(6)*v(486)-c(5)*v(495)+v(665)+v&
         &(666)+v(667)+v(668))+c(1)*(c(3)*v(459)+c(2)*v(468)+v(427)*v(495)+dcleft(6,2)*v(659)+v(669)+v(670)))+v(545)*(dcbleft(1,2&
         &)*v(315)+dcbright(1,2)*v(342)+dcb2(1,2,2)*v(420)+dcbright(5,2)*v(437)-dcbright(4,2)*v(438)-v(419)*v(522)+v(418)*v(531)&
         &+dcbleft(5,2)*v(559)-dcbleft(4,2)*v(562)+cb(5)*(-(cb(5)*v(504))+cb(6)*v(522)-cb(2)*v(531)+cb(4)*v(540)+v(671)+v(672)+v&
         &(673)+v(674))-cb(4)*(cb(4)*v(513)+cb(3)*v(522)-cb(6)*v(531)-cb(5)*v(540)+v(675)+v(676)+v(677)+v(678))+cb(1)*(cb(3)*v&
         &(504)+cb(2)*v(513)+v(439)*v(540)+dcbleft(6,2)*v(660)+v(679)+v(680)))
    v(336)=v(322)*v(334)+v(330)*v(658)+v(269)*(dcleft(1,1)*v(299)+dcright(5,2)*v(422)-dcright(4,2)*v(423)+dcright(1,2)*v&
         &(424)+v(417)*v(447)-v(416)*v(474)+v(415)*v(483)+dcleft(5,1)*v(547)-dcleft(4,1)*v(550)+c(5)*(-(c(5)*v(456))+c(6)*v(474)&
         &-c(2)*v(483)+c(4)*v(492)+v(681)+v(682)+v(683)+v(684))-c(4)*(c(4)*v(465)+c(3)*v(474)-c(6)*v(483)-c(5)*v(492)+v(685)+v&
         &(686)+v(687)+v(688))+c(1)*(c(3)*v(456)+c(2)*v(465)+v(427)*v(492)+dcleft(6,1)*v(659)+v(689)+v(690)))+v(545)*(dcbleft(1,1&
         &)*v(315)+dcb2(1,1,2)*v(420)+dcbright(5,2)*v(434)-dcbright(4,2)*v(435)+dcbright(1,2)*v(436)-v(419)*v(519)+v(418)*v(528)&
         &+dcbleft(5,1)*v(559)-dcbleft(4,1)*v(562)+cb(5)*(-(cb(5)*v(501))+cb(6)*v(519)-cb(2)*v(528)+cb(4)*v(537)+v(691)+v(692)+v&
         &(693)+v(694))-cb(4)*(cb(4)*v(510)+cb(3)*v(519)-cb(6)*v(528)-cb(5)*v(537)+v(695)+v(696)+v(697)+v(698))+cb(1)*(cb(3)*v&
         &(501)+cb(2)*v(510)+v(439)*v(537)+dcbleft(6,1)*v(660)+v(699)+v(700)))
    v(323)=v(303)*v(321)+v(319)*v(635)
    v(359)=v(323)*v(356)+v(269)*(dcleft(1,3)*v(300)+dcright(1,3)*v(349)+dcright(5,3)*v(430)-dcright(4,3)*v(431)+v(417)*v&
         &(454)-v(416)*v(481)+v(415)*v(490)+dcleft(5,3)*v(548)-dcleft(4,3)*v(551)+c(5)*(-(c(5)*v(463))+c(6)*v(481)-c(2)*v(490)+c&
         &(4)*v(499)+v(701)+v(702)+v(703)+v(704))-c(4)*(c(4)*v(472)+c(3)*v(481)-c(6)*v(490)-c(5)*v(499)+v(705)+v(706)+v(707)+v&
         &(708))+c(1)*(c(3)*v(463)+c(2)*v(472)+v(427)*v(499)+v(709)+v(710)+v(711)))+v(545)*(dcbleft(1,3)*v(316)+dcbright(1,3)*v&
         &(353)+dcb2(1,3,3)*v(420)+dcbright(5,3)*v(442)-dcbright(4,3)*v(443)-v(419)*v(526)+v(418)*v(535)+dcbleft(5,3)*v(560)&
         &-dcbleft(4,3)*v(563)+cb(5)*(-(cb(5)*v(508))+cb(6)*v(526)-cb(2)*v(535)+cb(4)*v(544)+v(712)+v(713)+v(714)+v(715))-cb(4)*&
         &(cb(4)*v(517)+cb(3)*v(526)-cb(6)*v(535)-cb(5)*v(544)+v(716)+v(717)+v(718)+v(719))+cb(1)*(cb(3)*v(508)+cb(2)*v(517)+v&
         &(439)*v(544)+v(720)+v(721)+v(722)))+v(352)*v(723)
    v(348)=v(323)*v(345)+v(341)*v(723)+v(269)*(dcleft(1,2)*v(300)+dcright(1,3)*v(338)+dcright(5,3)*v(425)-dcright(4,3)*v&
         &(426)+v(417)*v(451)-v(416)*v(478)+v(415)*v(487)+dcleft(5,2)*v(548)-dcleft(4,2)*v(551)+c(5)*(-(c(5)*v(460))+c(6)*v(478)&
         &-c(2)*v(487)+c(4)*v(496)+v(726)+v(727)+v(728)+v(729))-c(4)*(c(4)*v(469)+c(3)*v(478)-c(6)*v(487)-c(5)*v(496)+v(730)+v&
         &(731)+v(732)+v(733))+c(1)*(c(3)*v(460)+c(2)*v(469)+v(427)*v(496)+dcleft(6,2)*v(724)+v(734)+v(735)))+v(545)*(dcbleft(1,2&
         &)*v(316)+dcbright(1,3)*v(342)+dcb2(1,2,3)*v(420)+dcbright(5,3)*v(437)-dcbright(4,3)*v(438)-v(419)*v(523)+v(418)*v(532)&
         &+dcbleft(5,2)*v(560)-dcbleft(4,2)*v(563)+cb(5)*(-(cb(5)*v(505))+cb(6)*v(523)-cb(2)*v(532)+cb(4)*v(541)+v(736)+v(737)+v&
         &(738)+v(739))-cb(4)*(cb(4)*v(514)+cb(3)*v(523)-cb(6)*v(532)-cb(5)*v(541)+v(740)+v(741)+v(742)+v(743))+cb(1)*(cb(3)*v&
         &(505)+cb(2)*v(514)+v(439)*v(541)+dcbleft(6,2)*v(725)+v(744)+v(745)))
    v(337)=v(323)*v(334)+v(330)*v(723)+v(269)*(dcleft(1,1)*v(300)+dcright(5,3)*v(422)-dcright(4,3)*v(423)+dcright(1,3)*v&
         &(424)+v(417)*v(448)-v(416)*v(475)+v(415)*v(484)+dcleft(5,1)*v(548)-dcleft(4,1)*v(551)+c(5)*(-(c(5)*v(457))+c(6)*v(475)&
         &-c(2)*v(484)+c(4)*v(493)+v(746)+v(747)+v(748)+v(749))-c(4)*(c(4)*v(466)+c(3)*v(475)-c(6)*v(484)-c(5)*v(493)+v(750)+v&
         &(751)+v(752)+v(753))+c(1)*(c(3)*v(457)+c(2)*v(466)+v(427)*v(493)+dcleft(6,1)*v(724)+v(754)+v(755)))+v(545)*(dcbleft(1,1&
         &)*v(316)+dcb2(1,1,3)*v(420)+dcbright(5,3)*v(434)-dcbright(4,3)*v(435)+dcbright(1,3)*v(436)-v(419)*v(520)+v(418)*v(529)&
         &+dcbleft(5,1)*v(560)-dcbleft(4,1)*v(563)+cb(5)*(-(cb(5)*v(502))+cb(6)*v(520)-cb(2)*v(529)+cb(4)*v(538)+v(756)+v(757)+v&
         &(758)+v(759))-cb(4)*(cb(4)*v(511)+cb(3)*v(520)-cb(6)*v(529)-cb(5)*v(538)+v(760)+v(761)+v(762)+v(763))+cb(1)*(cb(3)*v&
         &(502)+cb(2)*v(511)+v(439)*v(538)+dcbleft(6,1)*v(725)+v(764)+v(765)))
    v(324)=(v(261)*v(305)+v(317)/v(260))/3d0
    v(325)=(v(261)*v(306)+v(318)/v(260))/3d0
    v(326)=(v(261)*v(307)+v(319)/v(260))/3d0
    dcs2(1,1,1)=dcright(1,1)*v(268)+dcleft(1,1)*v(324)+c(1)*v(335)+v(254)*v(446)
    dcs2(1,1,2)=dcright(1,2)*v(268)+dcleft(1,1)*v(325)+c(1)*v(336)+v(254)*v(447)
    dcs2(1,1,3)=dcright(1,3)*v(268)+dcleft(1,1)*v(326)+c(1)*v(337)+v(254)*v(448)
    dcs2(1,2,1)=dcright(1,1)*v(271)+dcleft(1,2)*v(324)+c(1)*v(346)+v(254)*v(449)
    dcs2(1,2,2)=dcright(1,2)*v(271)+dcleft(1,2)*v(325)+c(1)*v(347)+v(254)*v(450)
    dcs2(1,2,3)=dcright(1,3)*v(271)+dcleft(1,2)*v(326)+c(1)*v(348)+v(254)*v(451)
    dcs2(1,3,1)=dcright(1,1)*v(273)+dcleft(1,3)*v(324)+c(1)*v(357)+v(254)*v(452)
    dcs2(1,3,2)=dcright(1,2)*v(273)+dcleft(1,3)*v(325)+c(1)*v(358)+v(254)*v(453)
    dcs2(1,3,3)=dcright(1,3)*v(273)+dcleft(1,3)*v(326)+c(1)*v(359)+v(254)*v(454)
    dcs2(2,1,1)=dcright(2,1)*v(268)+dcleft(2,1)*v(324)+c(2)*v(335)+v(254)*v(455)
    dcs2(2,1,2)=dcright(2,2)*v(268)+dcleft(2,1)*v(325)+c(2)*v(336)+v(254)*v(456)
    dcs2(2,1,3)=dcright(2,3)*v(268)+dcleft(2,1)*v(326)+c(2)*v(337)+v(254)*v(457)
    dcs2(2,2,1)=dcright(2,1)*v(271)+dcleft(2,2)*v(324)+c(2)*v(346)+v(254)*v(458)
    dcs2(2,2,2)=dcright(2,2)*v(271)+dcleft(2,2)*v(325)+c(2)*v(347)+v(254)*v(459)
    dcs2(2,2,3)=dcright(2,3)*v(271)+dcleft(2,2)*v(326)+c(2)*v(348)+v(254)*v(460)
    dcs2(2,3,1)=dcright(2,1)*v(273)+dcleft(2,3)*v(324)+c(2)*v(357)+v(254)*v(461)
    dcs2(2,3,2)=dcright(2,2)*v(273)+dcleft(2,3)*v(325)+c(2)*v(358)+v(254)*v(462)
    dcs2(2,3,3)=dcright(2,3)*v(273)+dcleft(2,3)*v(326)+c(2)*v(359)+v(254)*v(463)
    dcs2(3,1,1)=dcright(3,1)*v(268)+dcleft(3,1)*v(324)+c(3)*v(335)+v(254)*v(464)
    dcs2(3,1,2)=dcright(3,2)*v(268)+dcleft(3,1)*v(325)+c(3)*v(336)+v(254)*v(465)
    dcs2(3,1,3)=dcright(3,3)*v(268)+dcleft(3,1)*v(326)+c(3)*v(337)+v(254)*v(466)
    dcs2(3,2,1)=dcright(3,1)*v(271)+dcleft(3,2)*v(324)+c(3)*v(346)+v(254)*v(467)
    dcs2(3,2,2)=dcright(3,2)*v(271)+dcleft(3,2)*v(325)+c(3)*v(347)+v(254)*v(468)
    dcs2(3,2,3)=dcright(3,3)*v(271)+dcleft(3,2)*v(326)+c(3)*v(348)+v(254)*v(469)
    dcs2(3,3,1)=dcright(3,1)*v(273)+dcleft(3,3)*v(324)+c(3)*v(357)+v(254)*v(470)
    dcs2(3,3,2)=dcright(3,2)*v(273)+dcleft(3,3)*v(325)+c(3)*v(358)+v(254)*v(471)
    dcs2(3,3,3)=dcright(3,3)*v(273)+dcleft(3,3)*v(326)+c(3)*v(359)+v(254)*v(472)
    dcs2(4,1,1)=dcright(4,1)*v(268)+dcleft(4,1)*v(324)+c(4)*v(335)+v(254)*v(473)
    dcs2(4,1,2)=dcright(4,2)*v(268)+dcleft(4,1)*v(325)+c(4)*v(336)+v(254)*v(474)
    dcs2(4,1,3)=dcright(4,3)*v(268)+dcleft(4,1)*v(326)+c(4)*v(337)+v(254)*v(475)
    dcs2(4,2,1)=dcright(4,1)*v(271)+dcleft(4,2)*v(324)+c(4)*v(346)+v(254)*v(476)
    dcs2(4,2,2)=dcright(4,2)*v(271)+dcleft(4,2)*v(325)+c(4)*v(347)+v(254)*v(477)
    dcs2(4,2,3)=dcright(4,3)*v(271)+dcleft(4,2)*v(326)+c(4)*v(348)+v(254)*v(478)
    dcs2(4,3,1)=dcright(4,1)*v(273)+dcleft(4,3)*v(324)+c(4)*v(357)+v(254)*v(479)
    dcs2(4,3,2)=dcright(4,2)*v(273)+dcleft(4,3)*v(325)+c(4)*v(358)+v(254)*v(480)
    dcs2(4,3,3)=dcright(4,3)*v(273)+dcleft(4,3)*v(326)+c(4)*v(359)+v(254)*v(481)
    dcs2(5,1,1)=dcright(5,1)*v(268)+dcleft(5,1)*v(324)+c(5)*v(335)+v(254)*v(482)
    dcs2(5,1,2)=dcright(5,2)*v(268)+dcleft(5,1)*v(325)+c(5)*v(336)+v(254)*v(483)
    dcs2(5,1,3)=dcright(5,3)*v(268)+dcleft(5,1)*v(326)+c(5)*v(337)+v(254)*v(484)
    dcs2(5,2,1)=dcright(5,1)*v(271)+dcleft(5,2)*v(324)+c(5)*v(346)+v(254)*v(485)
    dcs2(5,2,2)=dcright(5,2)*v(271)+dcleft(5,2)*v(325)+c(5)*v(347)+v(254)*v(486)
    dcs2(5,2,3)=dcright(5,3)*v(271)+dcleft(5,2)*v(326)+c(5)*v(348)+v(254)*v(487)
    dcs2(5,3,1)=dcright(5,1)*v(273)+dcleft(5,3)*v(324)+c(5)*v(357)+v(254)*v(488)
    dcs2(5,3,2)=dcright(5,2)*v(273)+dcleft(5,3)*v(325)+c(5)*v(358)+v(254)*v(489)
    dcs2(5,3,3)=dcright(5,3)*v(273)+dcleft(5,3)*v(326)+c(5)*v(359)+v(254)*v(490)
    dcs2(6,1,1)=dcright(6,1)*v(268)+dcleft(6,1)*v(324)+c(6)*v(335)+v(254)*v(491)
    dcs2(6,1,2)=dcright(6,2)*v(268)+dcleft(6,1)*v(325)+c(6)*v(336)+v(254)*v(492)
    dcs2(6,1,3)=dcright(6,3)*v(268)+dcleft(6,1)*v(326)+c(6)*v(337)+v(254)*v(493)
    dcs2(6,2,1)=dcright(6,1)*v(271)+dcleft(6,2)*v(324)+c(6)*v(346)+v(254)*v(494)
    dcs2(6,2,2)=dcright(6,2)*v(271)+dcleft(6,2)*v(325)+c(6)*v(347)+v(254)*v(495)
    dcs2(6,2,3)=dcright(6,3)*v(271)+dcleft(6,2)*v(326)+c(6)*v(348)+v(254)*v(496)
    dcs2(6,3,1)=dcright(6,1)*v(273)+dcleft(6,3)*v(324)+c(6)*v(357)+v(254)*v(497)
    dcs2(6,3,2)=dcright(6,2)*v(273)+dcleft(6,3)*v(325)+c(6)*v(358)+v(254)*v(498)
    dcs2(6,3,3)=dcright(6,3)*v(273)+dcleft(6,3)*v(326)+c(6)*v(359)+v(254)*v(499)
  END SUBROUTINE mls_combinedfbar2old



!----------------------
!*** END MLS
!*** START RBF MESHLESS
!----------------------

!---------------------------
!*** RADIUS AND DERIVATIVES
!*** 
!*** DETERMINED
!*** CHK0
!---------------------------
  SUBROUTINE RBF_RADIUS(NDI,X,XI,R,DRDX,D2RDX2)
    IMPLICIT REAL(8) (A-H,O-Z)
    REAL(8),DIMENSION(NDI)::X,XI
    REAL(8),DIMENSION(NDI)::DRDX
    REAL(8),DIMENSION(NDI,NDI)::D2RDX2
    R=RNORM2(NDI,X-XI)
    DRDX=0.0D00
    D2RDX2=0.0D00
    IF(R.GT.1.0D-20)THEN
       DO ID=1,NDI
          DRDX(ID)=(X(ID)-XI(ID))/R
          DO JD=1,NDI
             D2RDX2(ID,JD)=DELTAK(ID,JD)/R-(1.0D00/(R**3.0D00))*(X(ID)-XI(ID))*(X(JD)-XI(JD))
          END DO
       END DO
    END IF
  END SUBROUTINE RBF_RADIUS

!--------------------------
!*** RADIAL BASIS FUNCTION
!*** BASIS FUNCTION
!*** chk0
!--------------------------
  SUBROUTINE RBF_C2(RMAX,R,A,DADR,D2ADR2)
    IMPLICIT REAL(8) (A-H,O-Z)
    REAL(8)::RMAX,R,A,DADR,D2ADR2
    A=0.0D00
    DADR=0.0D00
    D2ADR2=0.0D00
    IF(R.LE.RMAX)THEN
       A=((1.0d00-(R/RMAX))**4)*(1.0D00+4.0d00*(R/RMAX))
       DADR=(20.0d00*R*(R-RMAX)**3.0D00)/(RMAX**5.0D00)
       D2ADR2=20.0d00*(R-RMAX)**2.0D00*(4.0D00*R-RMAX)/(RMAX**5.0D00)
    END IF
  END SUBROUTINE RBF_C2

!------------------------------------
!*** RADIAL AND DERIVATIVES
!*** XI IS THE "NODE" COORDINATE
!*** X IS THE GAUSS POINT COORDINATE
!*** chk0
!------------------------------------
  SUBROUTINE RBF_ONE(RMAX,NDI,X,XI,A,DADX,D2ADX2)
    IMPLICIT REAL(8) (A-H,O-Z)
    REAL(8)::RMAX,R,A,DADR,D2ADR2
    REAL(8),DIMENSION(NDI)::X,XI
    REAL(8),DIMENSION(NDI)::DADX,DRDX
    REAL(8),DIMENSION(NDI,NDI)::D2ADX2,D2RDX2
    CALL RBF_RADIUS(NDI,X,XI,R,DRDX,D2RDX2)
    CALL RBF_C2(RMAX,R,A,DADR,D2ADR2)
    DO ID=1,NDI
       DADX(ID)=DADR*DRDX(ID)
       DO JD=1,NDI
          D2ADX2(ID,JD)=D2ADR2*DRDX(ID)*DRDX(JD)+DADR*D2RDX2(ID,JD)
       END DO
    END DO
  END SUBROUTINE RBF_ONE

!-------------------------------------------
!*** OBTAIN SHAPE FUNCTIONS AND DERIVATIVES
!*** FROM A GIVEN LIST OF NODES
!*** AT ONE POINT WITH COORDINATES X
!*** CHK0
!-------------------------------------------
  SUBROUTINE RBF_SHAPEFUNCTIONS(IMESHLESS,RMAX,NDI,N,XN,X,XC,FF,DFF,DFF2)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER,PARAMETER::MPOL=20
    REAL(8)::RMAX
    REAL(8),DIMENSION(N,N)::AMATRIX,INVMATRIX,GMATRIX
    REAL(8),DIMENSION(N,MPOL)::BMATRIX
    REAL(8),DIMENSION(MPOL,N)::HMATRIX
    REAL(8),DIMENSION(NDI,NDI,N)::DFF2
    REAL(8),DIMENSION(NDI,N)::DFF
    REAL(8),DIMENSION(N)::FF,A
    REAL(8),DIMENSION(NDI,N)::XN
    REAL(8),DIMENSION(NDI)::X,XC
    REAL(8),DIMENSION(MPOL)::B
    REAL(8),DIMENSION(NDI,MPOL)::DBDX
    REAL(8),DIMENSION(NDI,NDI,MPOL)::DBDX2
    REAL(8),DIMENSION(NDI,NDI,MPOL)::D2BDX2
    REAL(8),DIMENSION(MPOL,MPOL)::B1B,INVB1B
    REAL(8),DIMENSION(NDI,N)::DADX
    REAL(8),DIMENSION(NDI,NDI,N)::D2ADX2
    REAL(8),DIMENSION(MPOL)::POLYN
    REAL(8),DIMENSION(NDI,MPOL)::DPOLYN
    REAL(8),DIMENSION(NDI,NDI,MPOL)::DPOLYN2
!------------------
!*** FORMS AMATRIX
!------------------
    DO JN=1,N
       DO IN=1,N
!-----------------------------------
!*** DADX AND D2ADX2 ARE TRASH here
!-----------------------------------
          CALL RBF_ONE(RMAX,NDI,XN(1:NDI,IN),XN(1:NDI,JN),AMATRIX(IN,JN),DADX(1:NDI,IN),D2ADX2(1:NDI,1:NDI,IN))
       END DO
    END DO
!-------------------------
!*** INVMATRIX=AMATRIX^-1
!-------------------------
    CALL INVMAT(N,DETA,AMATRIX,INVMATRIX)
!------------------
!*** FORMS BMATRIX
!------------------
    DO IN=1,N
       MP=0
       SELECT CASE(IMESHLESS)
       CASE(1)
          CALL MLS_POLYNLINBASE(NDI,1.0d00,XN(1:NDI,IN),XC(1:NDI),POLYN(1:MPOL),DPOLYN(1:NDI,1:MPOL),DPOLYN2(1:NDI,1:NDI,1:MPOL),MP)
       CASE(2)
          CALL MLS_POLYNQBASE(NDI,1.0d00,XN(1:NDI,IN),XC(1:NDI),POLYN(1:MPOL),DPOLYN(1:NDI,1:MPOL),DPOLYN2(1:NDI,1:NDI,1:MPOL),MP)
       CASE(3)
          CALL MLS_POLYNCUBBASE(NDI,1.0d00,XN(1:NDI,IN),XC(1:NDI),POLYN(1:MPOL),DPOLYN(1:NDI,1:MPOL),DPOLYN2(1:NDI,1:NDI,1:MPOL),MP)
       END SELECT
       DO IP=1,MP
          BMATRIX(IN,IP)=POLYN(IP)
       END DO
    END DO
!-----------------------
!*** FORMS BIG MATRICES
!-----------------------
    DO JD=1,MP
       DO ID=1,MP
          B1B(ID,JD)=0.0D00
          DO IN=1,N
             DO JN=1,N
                B1B(ID,JD)=B1B(ID,JD)+BMATRIX(IN,ID)*INVMATRIX(IN,JN)*BMATRIX(JN,JD)
             END DO
          END DO
       END DO
    END DO
!----------------
!*** INVERT B1B
!*** INTO INVB1B
!----------------
    CALL INVMAT(MP,DET,B1B(1:MP,1:MP),INVB1B(1:MP,1:MP))
!------------
!*** HMATRIX
!------------
    DO IN=1,N
       DO ID=1,MP
          HMATRIX(ID,IN)=0.0D00          
          DO JN=1,N
             DO JD=1,MP
                HMATRIX(ID,IN)=HMATRIX(ID,IN)+INVB1B(ID,JD)*BMATRIX(JN,JD)*INVMATRIX(JN,IN)
             END DO
          END DO
       END DO
    END DO
!------------
!*** GMATRIX
!------------
    DO JN=1,N
       DO IN=1,N
          GMATRIX(IN,JN)=INVMATRIX(IN,JN)
          DO ID=1,MP
             DO KN=1,N
                GMATRIX(IN,JN)=GMATRIX(IN,JN)-INVMATRIX(IN,KN)*BMATRIX(KN,ID)*HMATRIX(ID,JN)
             END DO
          END DO
       END DO
    END DO
!--------------------------
!*** FORMS VECTORS A AND B
!--------------------------
    DO IN=1,N
       CALL RBF_ONE(RMAX,NDI,X(1:NDI),XN(1:NDI,IN),A(IN),DADX(1:NDI,IN),D2ADX2(1:NDI,1:NDI,IN))
       SELECT CASE(IMESHLESS)
       CASE(1)
          CALL MLS_POLYNLINBASE(NDI,1.0d00,X(1:NDI),XC(1:NDI),B,DBDX(1:NDI,1:MP),DBDX2(1:NDI,1:NDI,1:MP),MP)
       CASE(2)
          CALL MLS_POLYNQBASE(NDI,1.0d00,X(1:NDI),XC(1:NDI),B,DBDX(1:NDI,1:MP),DBDX2(1:NDI,1:NDI,1:MP),MP)
       CASE(3)
          CALL MLS_POLYNCUBBASE(NDI,1.0d00,X(1:NDI),XC(1:NDI),B,DBDX(1:NDI,1:MP),DBDX2(1:NDI,1:NDI,1:MP),MP)
       END SELECT
    END DO
    D2BDX2=0.0D00   
!------------------------------------------------------
!*** NOW DETERMINES RB SHAPE FUNCTIONS AND DERIVATIVES
!------------------------------------------------------
    IF(IMESHLESS.EQ.0)THEN
       DO IN=1,N
          FF(IN)=DOTPROD(N,A(1:N),INVMATRIX(1:N,IN))
          DO ID=1,NDI
             DFF(ID,IN)=DOTPROD(N,DADX(ID,1:N),INVMATRIX(1:N,IN))
             DO JD=1,NDI
                DFF2(ID,JD,IN)=DOTPROD(N,D2ADX2(ID,JD,1:N),INVMATRIX(1:N,IN))
             ENDDO
          ENDDO
       ENDDO
    ELSE
       DO IN=1,N
          FF(IN)=DOTPROD(N,A(1:N),GMATRIX(1:N,IN))+DOTPROD(MP,B(1:MP),HMATRIX(1:MP,IN))
          DO ID=1,NDI
             DFF(ID,IN)=DOTPROD(N,DADX(ID,1:N),GMATRIX(1:N,IN))+DOTPROD(MP,DBDX(ID,1:MP),HMATRIX(1:MP,IN))
             DO JD=1,NDI
                DFF2(ID,JD,IN)=DOTPROD(N,D2ADX2(ID,JD,1:N),GMATRIX(1:N,IN))+DOTPROD(MP,D2BDX2(ID,JD,1:MP),HMATRIX(1:MP,IN))
             ENDDO
          ENDDO
       ENDDO
    END IF
  END SUBROUTINE RBF_SHAPEFUNCTIONS

!-------------------------------------------------
!*** FROM AN ELEMENT SUPPORT ESTIMATE, DETERMINES
!*** GAUSS POINT SUPPORT RADIUS 
!-------------------------------------------------
  SUBROUTINE RBF_GAUSSLEVELSHAPEFUNCTIONS(IMESHLESS,RCENTROID,NDI,XG,XC,NTOT,XTOT,FF,DFF,DFF2)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER::IMESHLESS
    REAL(8)::RCENTROID,RGAUSS
    REAL(8),DIMENSION(NDI)::XG,XC
    REAL(8),DIMENSION(NTOT)::FF
    REAL(8),DIMENSION(NDI,NTOT)::DFF
    REAL(8),DIMENSION(NDI,NDI,NTOT)::DFF2
    REAL(8),DIMENSION(NDI,NTOT)::XTOT
    RGAUSS=RCENTROID!-RNORM2(NDI,XG-XC)
    CALL RBF_TESTINGALLREQUIRED(IMESHLESS,RGAUSS,NDI,NTOT,XTOT,XG,XC,FF,DFF,DFF2)    
  END SUBROUTINE RBF_GAUSSLEVELSHAPEFUNCTIONS

!--------------------------------------------------------------------
!*** FOR REPRESENTATION PURPOSES ONLY:
!*** GET ALL NODES, MAKES A SELECTION
!*** AND RETURNS ALL QUANTITIES FOR ALL NODES, NOT ONLY THE SELECTION
!--------------------------------------------------------------------
  SUBROUTINE RBF_TESTINGALLREQUIRED(IMESHLESS,RMAX,NDI,NTOT,XTOT,X,XC,FF,DFF,DFF2)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER::IMESHLESS
    INTEGER,PARAMETER::MPOL=20
    INTEGER,DIMENSION(:),ALLOCATABLE::LISTN
    REAL(8)::RMAX
    REAL(8),DIMENSION(NTOT)::FF,FFLOC
    REAL(8),DIMENSION(NDI,NTOT)::DFF,DFFLOC
    REAL(8),DIMENSION(NDI,NDI,NTOT)::DFF2,DFF2LOC
    REAL(8),DIMENSION(NDI,NTOT)::XTOT,XN
    REAL(8),DIMENSION(NDI)::X,XC
    REAL(8),DIMENSION(MPOL)::POLYN
    REAL(8),DIMENSION(NDI,MPOL)::DPOLYN
    REAL(8),DIMENSION(NDI,NDI,MPOL)::DPOLYN2
    REAL(8),DIMENSION(MPOL,NTOT)::U2
!-----------------------
!*** GETS CLOSEST NODES
!-----------------------
    CALL MLS_GETSNODES(RMAX,X,NDI,NTOT,XTOT,N,LISTN)
    DO I=1,N
       DO ID=1,NDI
          XN(ID,I)=XTOT(ID,LISTN(I))
       END DO
    END DO
!-----------------------------------
!*** CHECK FOR REPEATED COORDINATES
!-----------------------------------
    GOTO 313
    DO IN=1,N
       DO JN=IN+1,N
          IF(RNORM2(NDI,XN(1:NDI,IN)-XN(1:NDI,JN)).LE.1.0D-20)THEN
             WRITE(*,*)"NODE IN AND JN",IN,JN," ARE THE SAME"
             WRITE(*,*)"IN,X=",IN,XN(1:NDI,IN)
             WRITE(*,*)"JN,X=",JN,XN(1:NDI,JN)
          END IF
       END DO
    END DO
313 CONTINUE
!-------------------------------
!*** SHAPE FUNCTION DERIVATIVES
!-------------------------------
    IF(IMESHLESS.GE.0)THEN
       CALL RBF_SHAPEFUNCTIONS(IMESHLESS,RMAX,NDI,N,XN,X,XC,FFLOC,DFFLOC,DFF2LOC)
    ELSE
!---------  
!**** MLS
!---------
       TOL=1.0D-2
       CALL MLS_U2ATACOORDINATE(IMESHLESS,RMAX,TOL,NDI,N,XC,XN,X,U2(1:MPOL,1:N))
       SELECT CASE(IMESHLESS)
       CASE(-1)
          CALL MLS_POLYNLINBASE(NDI,RMAX,X(1:NDI),XC(1:NDI),POLYN(1:MPOL),DPOLYN(1:NDI,1:MPOL),DPOLYN2(1:NDI,1:NDI,1:MPOL),M)
       CASE(-2)
          CALL MLS_POLYNQBASE(NDI,RMAX,X(1:NDI),XC(1:NDI),POLYN(1:MPOL),DPOLYN(1:NDI,1:MPOL),DPOLYN2(1:NDI,1:NDI,1:MPOL),M)
       CASE(-3)
          CALL MLS_POLYNCUBBASE(NDI,RMAX,X(1:NDI),XC(1:NDI),POLYN(1:MPOL),DPOLYN(1:NDI,1:MPOL),DPOLYN2(1:NDI,1:NDI,1:MPOL),M)
       END SELECT
       CALL MLS_SFDER(NDI,N,M,POLYN(1:M),DPOLYN(1:NDI,1:M),DPOLYN2(1:NDI,1:NDI,1:M),U2(1:M,1:N),FFLOC,DFFLOC,DFF2LOC)
    END IF
    FF=0.0D00
    DFF=0.0D00
    DFF2=0.0D00
    DO I=1,N
       FF(LISTN(I))=FFLOC(I)
       DFF(1:NDI,LISTN(I))=DFFLOC(1:NDI,I)
       DFF2(1:NDI,1:NDI,LISTN(I))=DFF2LOC(1:NDI,1:NDI,I)
    END DO
    DEALLOCATE(LISTN)
  END SUBROUTINE RBF_TESTINGALLREQUIRED

!--------------------------
!*** END RBF, END MESHLESS
!--------------------------
