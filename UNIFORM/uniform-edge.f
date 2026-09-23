C=======================================================================
C  UNIFORM-MAIN.F
C
C  Heavily loaded actuator disk with UNIFORM load (h - h0 = b = const
C  inside the slipstream, J = 0, i.e. no swirl) in the framework of
C  Conway's exact actuator-disk theory (Conway, J.Fluid Mech. 297 (1995)
C  and 365 (1998)); companion program to Conway-main.f (parabolic load).
C
C  For h(Psi) = h0 + b*Theta(Psi_e - Psi) the vorticity   omega/r = -dh/dPsi
C  degenerates to a vortex sheet on the slipstream boundary r = R(z).
C  Its strength per unit axial length (ring circulation d(Gamma)/dz) is
C
C        gamma(z) = b / Vzm(z),     Vzm = mean axial velocity on the sheet
C
C  (Bernoulli jump condition), and the polynomial radial expansion of
C  Conway (eq. 6.59/6.60 of the book) is replaced by the ring-vortex
C  kernels  I(0,1,1), I(1,1,0), I(1,1,1)  (closed form, complete
C  elliptic integrals):
C
C    Psi(r,z) = r^2/2 + r/2 * int_0^inf gamma R' I(0,1,1)(R',r,z-z') dz'
C    Vz (r,z) =   1   + 1/2 * int_0^inf gamma R' I(1,1,0)(R',r,z-z') dz'
C    Vr (r,z) =         1/2 * int_0^inf sgn(z-z') gamma R' I(1,1,1)   dz'
C
C  R' = R(z'). The numerical integration runs over 0 <= z' <= X2; the
C  remainder X2..infinity is a cylinder of radius R(X2) and strength
C  gamma(X2) whose influence is given analytically by
C  I(-1,1,1), I(0,1,0), I(0,1,1) (same integrals as in GENINT of math.f,
C  here in closed form, see CYLINT).
C
C  Iteration (as in Conway-main.f):
C     Psi_e = Psi(1,0) ; Psi(R(z),z) = Psi_e  -->  R(z)  (quadratic)
C     gamma(z) = b / Vzm(z)
C  Slipstream radius R(z): natural cubic spline through the nodes,
C  gamma(z): piecewise linear.  On the sheet the log- and 1/(z-z')-
C  singularities of the kernels are subtracted analytically, so that
C  Vzm is the Cauchy principal value = mean of both sides.
C
C  Non-dimensional: U_inf = 1, R_disk = 1, rho = 1.
C     cT = 2 b            (propeller sign; wind turbine: b < 0, cT = -2b)
C     cP = 4 b Psi_e      (power = b * mass flux 2 pi Psi_e)
C     far wake: (1+w)^2 = 1 + 2b,  gamma_inf = w,  R_inf^2 (1+w) = 2 Psi_e
C
C  Input (screen):  b    (e.g. -0.44444 = Betz optimum cT = 8/9,
C                         b > -0.5 required; propeller b > 0)
C                   X2, NN, NITER  (end of numerical domain, number of
C                         nodes, max. iterations; 0 = default 40 40 30;
C                         for cT -> 1 (b -> -0.5) the wake develops
C                         slowly: use X2 = 100 or more)
C  Output files: Stream.DAT, Solution.DAT, Gamma.DAT, Psi.DAT,
C                Vz1..Vz10.DAT, Vr1..Vr10.DAT, AxVz0.DAT
C                (velocities in the profile files are perturbation
C                 velocities Vz-1, Vr as in Conway-main.f)
C
C  A.P. Schaffarczyk / Claude, HAW Kiel, September 2026
C=======================================================================
      PROGRAM UNIFORM
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=201)
      INTEGER NN,NITER,NPROF,L,I,KK,II(15),mode
      REAL*8 ZN(NMAX),RN(NMAX),GN(NMAX),ZK(0:NMAX),RK(0:NMAX),
     *       R2K(0:NMAX),VZM(NMAX)
      REAL*8 X2,GINF,RINF,pi,b,w,relax,relaxg,psie,z1
      REAL*8 z,rho,R0,G,Rq,Psi,Psi1,Vz,Vr,Re,crit1(1000),crit2(1000)
      REAL*8 rhoc,th,V2,pr,dl,Vn,Frp,Fzp,Frm,Fzm,Sv2
      INTEGER NTH
      REAL*8 Rnew(NMAX),gnew(NMAX),zz(15),Rinfex,psimt,a,cT,cP
      REAL*8 RADIUS,DRADIUS,GAMSH
      COMMON/SHEET/X2,GINF,RINF,ZN,RN,GN,ZK,RK,R2K,NN
      COMMON/FUN/pi
      EXTERNAL RADIUS,DRADIUS,GAMSH
C
      pi = 4.d0*DATAN(1.d0)
C
C  Numerical parameters
C
      X2     = 40.d0
      NN     = 40
      READ(5,*) z1
      NITER  = 30
      relax  = 0.7d0
      relaxg = 0.7d0
      NPROF  = 10
C
C  axial locations for the radial profiles (as in Conway-main.f)
C
      zz(1)  = -2.50d0
      zz(2)  = -2.0d0
      zz(3)  = -1.5d0
      zz(4)  = -1.0d0
      zz(5)  = -0.5d0
      zz(6)  =  0.0d0
      zz(7)  =  0.5d0
      zz(8)  =  1.0d0
      zz(9)  =  1.5d0
      zz(10) =  2.0d0
      zz(11) =  2.5d0
      II(1)  = 1
      II(2)  = 2
      II(3)  = 3
      II(4)  = 4
      II(5)  = 7
      II(6)  = 8
      II(7)  = 9
      II(8)  = 10
      II(9)  = 11
      II(10) = 12
      II(11) = 13
C
      OPEN(UNIT=14,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Stream.DAT')
      OPEN(UNIT=15,FORM='FORMATTED',STATUS='UNKNOWN',
     *     FILE='Solution.DAT')
C
      WRITE(6,*)'Enter b (= h-h0 inside slipstream, cT = 2b): '
      READ(5,*) b
      IF(b.LE.-0.5d0) THEN
         WRITE(6,*)'b must be > -0.5 (cT < 1 for the wind turbine)'
         STOP
      END IF
      w = DSQRT(1.d0+2.d0*b) - 1.d0
      WRITE(6,*)'Enter X2, NN, NITER (0 0 0 = defaults 40 40 30): '
      READ(5,*) z,I,L
      IF(z.GT.0.d0) X2 = z
      IF(I.GT.0.AND.I.LE.NMAX) NN = I
      IF(L.GT.0.AND.L.LE.1000) NITER = L
C
C  Nodes on the sheet (geometric), initial guess: cylinder R = 1,
C  gamma = far-wake value w
C
      DO 1 I = 1,NN
         ZN(I) = z1*(X2/z1)**(DFLOAT(I-1)/DFLOAT(NN-1))
         RN(I) = 1.d0
         GN(I) = w
    1 CONTINUE
      CALL SETSPL
      psie = 0.5d0*(1.d0+0.5d0*w)
C
C  Main iteration loop
C
      DO 100 L = 1,NITER
      WRITE(6,*)'Iteration number: ',L
      WRITE(14,*)'Iteration number: ',L
C
C  stream function at the disk edge = start of the sheet
C
      CALL FIELD(1.d0,0.d0,2,psie,Vz,Vr)
      WRITE(14,150) 0.d0,1.d0
C
C  streamline condition Psi(R(z),z) = Psi_e  -->  new R(z)
C
      DO 10 I = 1,NN
         z  = ZN(I)
         R0 = RN(I)
         CALL FIELD(R0,z,2,Psi,Vz,Vr)
         G  = (Psi - 0.5d0*R0*R0)*2.d0/R0
         Rq = -0.5d0*G + DSQRT(0.25d0*G*G + 2.d0*psie)
         Rnew(I) = R0 + relax*(Rq-R0)
         WRITE(14,150) z,Rq
   10 CONTINUE
      crit1(L) = 0.d0
      DO 11 I = 1,NN
         crit1(L) = crit1(L) + (Rnew(I)-RN(I))**2
         RN(I) = Rnew(I)
   11 CONTINUE
      crit1(L) = DSQRT(crit1(L)/NN)
      CALL SETSPL
C
C  Bernoulli jump condition: gamma = b / Vzm
C
      DO 20 I = 1,NN
         z  = ZN(I)
         R0 = RN(I)
         CALL FIELD(R0,z,3,Psi,Vz,Vr)
         VZM(I)  = Vz
         gnew(I) = b/Vz
   20 CONTINUE
      crit2(L) = 0.d0
      DO 21 I = 1,NN
         crit2(L) = crit2(L) + (gnew(I)-GN(I))**2
         GN(I) = GN(I) + relaxg*(gnew(I)-GN(I))
   21 CONTINUE
      crit2(L) = DSQRT(crit2(L)/NN)
      GINF = GN(NN)
      WRITE(6,160) psie,RN(NN),crit1(L),crit2(L)
  160 FORMAT(1X,'psie =',F10.6,'  R(X2) =',F10.6,'  err1 =',E10.3,
     *       '  err2 =',E10.3)
      IF(crit1(L).LT.1.d-6.AND.crit2(L).LT.1.d-6) GO TO 101
  100 CONTINUE
      L = NITER
  101 CONTINUE
      NITER = L
C
C  Output convergence history and global parameters
C
      WRITE(15,120)
      WRITE(15,125)
      DO 110 L = 1,NITER
      WRITE(15,130)L,crit1(L),crit2(L)
  110 CONTINUE
      Rinfex = DSQRT(2.d0*psie/(1.d0+w))
      psimt  = 0.5d0*(1.d0+0.5d0*w)
      cT = 2.d0*b
      cP = 4.d0*b*psie
      IF(b.GE.0.d0) THEN
         WRITE(15,*)'***** Propeller Case *****'
      ELSE
         WRITE(15,*)'***** Wind Turbine Case *****'
      END IF
      WRITE(15,*)'Uniform load: h - h0 = b inside slipstream, J = 0'
      WRITE(15,200) b
      WRITE(15,205) relax
      WRITE(15,210) psie
      WRITE(15,215) psimt
      WRITE(15,220) 2.d0*psie
      WRITE(15,225) w
      WRITE(15,230) Rinfex
      WRITE(15,235) RN(NN),X2
      WRITE(15,240) GINF
      WRITE(15,245) cT
      WRITE(15,250) cP
      WRITE(15,255) cT*(1.d0+0.5d0*w)
      IF(b.LT.0.d0) THEN
         a = -0.5d0*w
         WRITE(15,260) -cT,-cP,a,-4.d0*a*(1.d0-a)**2
      END IF
  120 FORMAT (1X,'CONVERGENCE HISTORY')
  125 FORMAT (1X,'ITERATION      ERROR1(R)    ERROR2(gamma)')
  130 FORMAT (I8,3X,E12.4,1X,E12.4)
  150 FORMAT (2F10.6)
  155 FORMAT (3F10.6)
  200 FORMAT(1X,'b                    = ',F12.6)
  205 FORMAT(1X,'relax                = ',F12.6)
  210 FORMAT(1X,'psi-edge             = ',F12.6)
  215 FORMAT(1X,'psi-edge mom. theory = ',F12.6)
  220 FORMAT(1X,'mean Vz at disk 2psie= ',F12.6)
  225 FORMAT(1X,'w (far wake, exact)  = ',F12.6)
  230 FORMAT(1X,'R_inf = sqrt(2psie/(1+w)) = ',F12.6)
  235 FORMAT(1X,'R(X2)                = ',F12.6,'   X2 =',F8.2)
  240 FORMAT(1X,'gamma(X2)  (-> w)    = ',F12.6)
  245 FORMAT(1X,'cT = 2b              = ',F12.6)
  250 FORMAT(1X,'cP = 4 b psie        = ',F12.6)
  255 FORMAT(1X,'cP mom. theory       = ',F12.6)
  260 FORMAT(1X,'wind turbine sign: cT = ',F10.6,'  cP = ',F10.6,
     *       /,1X,'   a = ',F10.6,'  cP(mom.th.) = -4a(1-a)^2 = ',F10.6)
C
      DO 111 KK = 1,NPROF
      WRITE(15,'(a24,f10.4)')'Radial Profile at z/Ra =',zz(KK)
  111 CONTINUE
      CLOSE(UNIT=14)
      CLOSE(UNIT=15)
C
C  Sheet data: z, R, dR/dz, gamma, Vzm
C
      OPEN(UNIT=16,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Gamma.DAT')
      WRITE(16,*)'#     z          R        dR/dz      gamma      Vzm'
      DO 30 I = 1,NN
         WRITE(16,165) ZN(I),RN(I),DRADIUS(ZN(I)),GN(I),VZM(I)
   30 CONTINUE
  165 FORMAT(5F11.6)
      CLOSE(UNIT=16)
C
C
C  Momentum flux through small circles around the disk edge (1,0)
C  (meridional plane, per unit rim length): "edge force" test
C     F = int [ p n + V (V.n) ] dl ,  p = h - V^2/2,  h = 1/2 + b inside
C
      OPEN(UNIT=17,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Edge.DAT')
      OPEN(UNIT=18,FORM='FORMATTED',STATUS='UNKNOWN',FILE='EdgeV.DAT')
      WRITE(18,*)'#  rho   theta   r   z   Vr   Vz   |V|   Psi-Psie'
      WRITE(17,*)'#  rho   Fr   Fz   Fr(pres)   Fz(pres)',
     *           '   Fr(mom)    Fz(mom)   int|V|^2dl'
      DO 400 KK = 1,6
         IF(KK.EQ.1) rhoc = 0.05d0
         IF(KK.EQ.2) rhoc = 0.02d0
         IF(KK.EQ.3) rhoc = 0.01d0
         IF(KK.EQ.4) rhoc = 0.005d0
         IF(KK.EQ.5) rhoc = 0.002d0
         IF(KK.EQ.6) rhoc = 0.001d0
         NTH = 720
         Frp = 0.d0
         Fzp = 0.d0
         Frm = 0.d0
         Fzm = 0.d0
         Sv2 = 0.d0
         DO 410 I = 1,NTH
            th  = 2.d0*pi*(I-0.5d0)/NTH
            rho = 1.d0 + rhoc*DCOS(th)
            z   = rhoc*DSIN(th)
            CALL FIELD(rho,z,0,Psi,Vz,Vr)
            V2 = Vz*Vz + Vr*Vr
            pr = 0.5d0*(1.d0 - V2)
            IF(Psi.LT.psie) pr = pr + b
            dl = 2.d0*pi*rhoc/NTH
            Vn = Vr*DCOS(th) + Vz*DSIN(th)
            Frp = Frp + pr*DCOS(th)*dl
            Fzp = Fzp + pr*DSIN(th)*dl
            Frm = Frm + Vr*Vn*dl
            Fzm = Fzm + Vz*Vn*dl
            Sv2 = Sv2 + V2*dl
            IF(MOD(I,60).EQ.30) WRITE(18,416) rhoc,th*180.d0/pi,rho,z,
     *          Vr,Vz,DSQRT(V2),Psi-psie
  416       FORMAT(F8.4,F8.1,2F10.5,4F11.5)
  410    CONTINUE
         WRITE(17,415) rhoc,Frp+Frm,Fzp+Fzm,Frp,Fzp,Frm,Fzm,Sv2
         WRITE(6,415) rhoc,Frp+Frm,Fzp+Fzm,Frp,Fzp,Frm,Fzm,Sv2
  415    FORMAT(F8.4,7E12.4)
  400 CONTINUE
      CLOSE(UNIT=17)
      NPROF = 0
C
C  Stream function on radial profiles (as Psi.DAT of Conway-main.f)
C
      OPEN(UNIT=21,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Psi.DAT')
      DO 95 KK = 1,NPROF
         z = zz(KK)
         IF(z.LE.0.d0) THEN
            Re = 1.d0
         ELSE
            Re = RADIUS(z)
         END IF
         WRITE(21,*)
         WRITE(21,155) 0.d0,z,0.d0
         DO 90 I = 1,151
            rho = 0.01d0*I*Re
            mode = 0
            IF(z.GE.0.d0.AND.DABS(rho-RADIUS(z)).LT.1.d-9) mode = 2
            CALL FIELD(rho,z,mode,Psi,Vz,Vr)
            WRITE(21,155) rho,z,Psi
   90    CONTINUE
   95 CONTINUE
      CLOSE(UNIT=21)
C
C  Axial velocity profiles (perturbation Vz - 1)
C
      OPEN(UNIT=1,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vz1.DAT')
      OPEN(UNIT=2,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vz2.DAT')
      OPEN(UNIT=3,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vz3.DAT')
      OPEN(UNIT=4,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vz4.DAT')
      OPEN(UNIT=7,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vz5.DAT')
      OPEN(UNIT=8,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vz6.DAT')
      OPEN(UNIT=9,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vz7.DAT')
      OPEN(UNIT=10,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vz8.DAT')
      OPEN(UNIT=11,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vz9.DAT')
      OPEN(UNIT=12,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vz10.DAT')
      OPEN(UNIT=13,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vz11.DAT')
      DO 40 KK = 1,NPROF
         z = zz(KK)
         IF(z.LE.0.d0) THEN
            Re = 1.d0
         ELSE
            Re = RADIUS(z)
         END IF
         WRITE(6,*)'Vz profile at z = ',z
         DO 35 I = 0,170
            rho = 0.01d0*I*Re
            mode = 0
            IF(z.GE.0.d0.AND.DABS(rho-RADIUS(z)).LT.1.d-9) mode = 1
            IF(mode.EQ.1.AND.z.EQ.0.d0) THEN
C  disk edge: Vz log-singular, mean value from finite difference of Psi
C  (as in Conway-main.f)
               CALL FIELD(rho*1.01d0,z,0,Psi,Vz,Vr)
               CALL FIELD(rho*0.99d0,z,0,Psi1,Vz,Vr)
               Vz = (Psi-Psi1)/(0.02d0*rho*rho)
            ELSE
               CALL FIELD(rho,z,mode,Psi,Vz,Vr)
            END IF
            WRITE(II(KK),150) rho,Vz-1.d0
   35    CONTINUE
         CLOSE(UNIT=II(KK))
   40 CONTINUE
C
C  Radial velocity profiles
C
      OPEN(UNIT=1,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vr1.DAT')
      OPEN(UNIT=2,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vr2.DAT')
      OPEN(UNIT=3,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vr3.DAT')
      OPEN(UNIT=4,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vr4.DAT')
      OPEN(UNIT=7,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vr5.DAT')
      OPEN(UNIT=8,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vr6.DAT')
      OPEN(UNIT=9,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vr7.DAT')
      OPEN(UNIT=10,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vr8.DAT')
      OPEN(UNIT=11,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vr9.DAT')
      OPEN(UNIT=12,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vr10.DAT')
      OPEN(UNIT=13,FORM='FORMATTED',STATUS='UNKNOWN',FILE='Vr11.DAT')
      DO 70 KK = 1,NPROF
         z = zz(KK)
         WRITE(6,*)'Vr profile at z = ',z
         DO 60 I = 1,170
            rho = 0.01d0*I
            mode = 0
            IF(z.GT.0.d0.AND.DABS(rho-RADIUS(z)).LT.1.d-9) mode = 1
            IF(z.EQ.0.d0.AND.DABS(rho-1.d0).LT.1.d-9) THEN
C  disk edge: Vr log-singular, value not defined
               Vr = 0.d0
            ELSE
               CALL FIELD(rho,z,mode,Psi,Vz,Vr)
            END IF
            WRITE(II(KK),150) rho,Vr
   60    CONTINUE
         CLOSE(UNIT=II(KK))
   70 CONTINUE
C
C  Axial velocity on the axis (exact, no extrapolation needed)
C
      OPEN(UNIT=1,FORM='FORMATTED',STATUS='UNKNOWN',FILE='AxVz0.DAT')
      DO 270 I = 1,801
         z = 0.01d0*(I-401)
         CALL FIELD(0.d0,z,0,Psi,Vz,Vr)
         WRITE(1,150) z,Vz-1.d0
  270 CONTINUE
      CLOSE(UNIT=1)
      STOP
      END
C=======================================================================
C  Field quantities at (rho,z): stream function Psi, total axial
C  velocity Vz, radial velocity Vr.
C   mode = 0 : general point off the sheet (Psi, Vz, Vr)
C   mode = 1 : point ON the sheet, rho = R(z), 0 < z <= X2.  The
C              log|z-z'| and 1/(z-z') singularities are subtracted
C              analytically; Vz is then the mean velocity on the sheet
C              (principal value) and Vr = dR/dz * Vz (streamline).
C   mode = 2 : as 1, but Psi only (also used at the disk edge z = 0,
C              where Vz is log-singular for dR/dz(0+) <> 0)
C   mode = 3 : as 1, but Vz only
C=======================================================================
      SUBROUTINE FIELD(rho,z,mode,Psi,Vz,Vr)
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=201)
      INTEGER NN,mode,kind
      REAL*8 ZN(NMAX),RN(NMAX),GN(NMAX),ZK(0:NMAX),RK(0:NMAX),
     *       R2K(0:NMAX)
      REAL*8 X2,GINF,RINF,pi,rho,z,Psi,Vz,Vr
      REAL*8 rhof,zf,CL,CP,zeta,g0,R0,Rz0,res,PT,VZT,VRT,zc,xlnx
      REAL*8 GAMSH,DRADIUS,FUNC
      COMMON/SHEET/X2,GINF,RINF,ZN,RN,GN,ZK,RK,R2K,NN
      COMMON/FUN/pi
      COMMON/FLD/rhof,zf,CL,CP,kind
      EXTERNAL FUNC,GAMSH,DRADIUS
C
      rhof = rho
      zf   = z
      zc   = DMIN1(DMAX1(z,0.d0),X2)
      zeta = X2 - z
      CALL CYLINT(RINF,rho,zeta,PT,VZT,VRT)
      PT  = GINF*RINF*PT
      VZT = GINF*RINF*VZT
      VRT = -GINF*RINF*VRT
      g0  = 0.d0
      R0  = rho
      Rz0 = 0.d0
      IF(mode.GE.1) THEN
         g0  = GAMSH(z)
         IF(z.LT.X2) Rz0 = DRADIUS(z)
      END IF
      Psi = 0.d0
      Vz  = 0.d0
      Vr  = 0.d0
C
C  stream function
C
      IF(mode.NE.3) THEN
      kind = 1
      CL = 0.d0
      CP = 0.d0
      IF(mode.GE.1) CL = -g0/pi
      CALL QSING(FUNC,zc,res)
      IF(mode.GE.1) res = res + CL*(xlnx(z)-z+xlnx(X2-z)-(X2-z))
      Psi = 0.5d0*rho*rho + 0.5d0*rho*(res+PT)
      END IF
      IF(mode.EQ.2) RETURN
C
C  axial velocity
C
      kind = 2
      CL = 0.d0
      CP = 0.d0
      IF(mode.GE.1) THEN
         CL = -g0/(2.d0*pi*R0)
         CP = -g0*Rz0/(pi*(1.d0+Rz0*Rz0))
      END IF
      CALL QSING(FUNC,zc,res)
      IF(mode.GE.1) THEN
         res = res + CL*(xlnx(z)-z+xlnx(X2-z)-(X2-z))
         IF(CP.NE.0.d0) res = res + CP*DLOG(z/(X2-z))
      END IF
      Vz = 1.d0 + 0.5d0*(res+VZT)
C
C  radial velocity
C
      IF(mode.GE.1) THEN
         Vr = Rz0*Vz
      ELSE IF(rho.LT.1.d-10) THEN
         Vr = 0.d0
      ELSE
         kind = 3
         CL = 0.d0
         CP = 0.d0
         CALL QSING(FUNC,zc,res)
         Vr = 0.5d0*(res+VRT)
      END IF
      RETURN
      END
C=======================================================================
      FUNCTION xlnx(x)
      IMPLICIT NONE
      REAL*8 xlnx,x
      IF(x.GT.0.d0) THEN
         xlnx = x*DLOG(x)
      ELSE
         xlnx = 0.d0
      END IF
      RETURN
      END
C=======================================================================
C  Integrand for the axial integration (z' = x) at the field point
C  (rhof,zf):  kind = 1 Psi, 2 Vz, 3 Vr ; CL, CP subtraction terms
C=======================================================================
      FUNCTION FUNC(x)
      IMPLICIT NONE
      INTEGER kind
      REAL*8 FUNC,x,rhof,zf,CL,CP,Rx,gx,zeta,P,VZ,VR,f
      REAL*8 RADIUS,GAMSH
      COMMON/FLD/rhof,zf,CL,CP,kind
      EXTERNAL RADIUS,GAMSH
   
