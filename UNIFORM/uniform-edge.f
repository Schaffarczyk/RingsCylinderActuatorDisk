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
      Rx   = RADIUS(x)
      gx   = GAMSH(x)
      zeta = zf - x
      IF(zeta.EQ.0.d0) THEN
         FUNC = 0.d0
         RETURN
      END IF
      CALL RING(Rx,rhof,zeta,P,VZ,VR)
      IF(kind.EQ.1) THEN
         f = gx*Rx*P
      ELSE IF(kind.EQ.2) THEN
         f = gx*Rx*VZ
      ELSE
         f = DSIGN(1.d0,zeta)*gx*Rx*VR
      END IF
      IF(CL.NE.0.d0) f = f - CL*DLOG(DABS(zeta))
      IF(CP.NE.0.d0) f = f - CP/zeta
      FUNC = f
      RETURN
      END
C=======================================================================
C  Ring-vortex kernels (Bessel-Laplace integrals, Conway's notation
C  I(lambda,mu,nu)(R,rho,z) = int_0^inf s^lambda J_mu(sR) J_nu(s rho)
C  exp(-s|z|) ds ) in closed form:
C     P  = I(0,1,1)   stream function kernel
C     VZ = I(1,1,0)   axial velocity kernel
C     VR = I(1,1,1)   radial velocity kernel
C  (ring of radius R, unit circulation, field point (rho,z):
C   Psi = rho R/2 * P,  Vz = R/2 * VZ,  Vr = sgn(z) R/2 * VR)
C=======================================================================
      SUBROUTINE RING(R,rho,z,P,VZ,VR)
      IMPLICIT NONE
      REAL*8 R,rho,z,P,VZ,VR,pi,D2,d2s,kp2,k2,k,Dd,AK,AE,rf,rd
      COMMON/FUN/pi
      IF(rho.LT.1.d-10) THEN
         P  = 0.d0
         VZ = R/(R*R+z*z)**1.5d0
         VR = 0.d0
         RETURN
      END IF
      D2  = (R+rho)**2 + z*z
      d2s = (R-rho)**2 + z*z
      kp2 = d2s/D2
      k2  = 4.d0*R*rho/D2
      k   = DSQRT(k2)
      Dd  = DSQRT(D2)
      AK  = rf(0.d0,kp2,1.d0)
      AE  = AK - k2*rd(0.d0,kp2,1.d0)/3.d0
      P   = ((2.d0/k-k)*AK - 2.d0/k*AE)/(pi*DSQRT(R*rho))
      VZ  = (AK + (R*R-rho*rho-z*z)/d2s*AE)/(pi*R*Dd)
      VR  = DABS(z)/(pi*R*rho*Dd)*(-AK + (R*R+rho*rho+z*z)/d2s*AE)
      RETURN
      END
C=======================================================================
C  Semi-infinite vortex cylinder (rings of radius R from axial
C  distance zeta >= 0 to infinity, unit strength per unit length):
C     PM  = I(-1,1,1)(R,rho,zeta) = int_zeta^inf I(0,1,1) dt
C     P010= I( 0,1,0)             = int_zeta^inf I(1,1,0) dt
C     P011= I( 0,1,1)             = int_zeta^inf I(1,1,1) dt
C  Closed forms with K, E and Heuman's Lambda function; identical to
C  IN(-1,1,1), IN(0,1,0), IN(0,1,1) of GENINT (math.f) but without the
C  recursions (no overflow for large zeta / small rho).
C=======================================================================
      SUBROUTINE CYLINT(R,rho,zeta,PM,P010,P011)
      IMPLICIT NONE
      REAL*8 R,rho,zeta,PM,P010,P011,pi
      REAL*8 D2,d2s,kp2,k2,k,kp,Dd,AK,AE,F0,Beta,EF,FF,alam,w
      REAL*8 I000,I011,I010,I001,I100,I022,I122,I021,rf,rd,ellf,elle
      COMMON/FUN/pi
      IF(rho.LT.1.d-10) THEN
         PM   = 0.d0
         P010 = (1.d0 - zeta/DSQRT(R*R+zeta*zeta))/R
         P011 = 0.d0
         RETURN
      END IF
      IF(zeta.LT.1.d-12.AND.DABS(rho-R).LT.1.d-10) THEN
C  field point on the rim of the cylinder: PM = 1/2, P010 = 1/(2R)
         PM   = 0.5d0
         P010 = 0.5d0/R
         P011 = 0.d0
         RETURN
      END IF
      D2  = (R+rho)**2 + zeta*zeta
      d2s = (R-rho)**2 + zeta*zeta
      kp2 = d2s/D2
      k2  = 4.d0*R*rho/D2
      k   = DSQRT(k2)
      kp  = DSQRT(kp2)
      Dd  = DSQRT(D2)
      AK  = rf(0.d0,kp2,1.d0)
      AE  = AK - k2*rd(0.d0,kp2,1.d0)/3.d0
      w   = (R*R+rho*rho+zeta*zeta)/(2.d0*rho*R)
C  Heuman's Lambda function
      F0  = DSQRT(d2s)
      Beta = DASIN(zeta/F0)
      EF = elle(Beta,kp)
      FF = ellf(Beta,kp)
      alam = 2.d0*(AK*EF-(AK-AE)*FF)/pi
C  I(0,0,0), I(0,1,1)
      I000 = k*AK/(pi*DSQRT(rho*R))
      I011 = ((2.d0-k2)*AK-2.d0*AE)/(k*pi*DSQRT(rho*R))
C  I(0,1,0), I(0,0,1)
      IF(rho.LT.R) THEN
         I010 = (-zeta*k*AK/(2.d0*pi*DSQRT(rho*R)) - 0.5d0*alam + 1.d0)
     *          /R
      ELSE
         I010 = (-zeta*k*AK/(2.d0*pi*DSQRT(rho*R)) + 0.5d0*alam)/R
      END IF
      IF(R.LT.rho) THEN
         I001 = (-zeta*k*AK/(2.d0*pi*DSQRT(rho*R)) - 0.5d0*alam + 1.d0)
     *          /rho
      ELSE
         I001 = (-zeta*k*AK/(2.d0*pi*DSQRT(rho*R)) + 0.5d0*alam)/rho
      END IF
C  I(1,0,0), I(0,2,2), I(1,2,2), I(0,2,1)  (recursions of GENINT)
      I100 = zeta*k**3*AE/(4.d0*pi*kp2*(rho*R)**1.5d0)
      I022 = (4.d0*w*I011 - I000)/3.d0
      I122 = -0.125d0*zeta*k**4/(R*rho*kp2)*3.d0*(w*I022-I011)
      I021 = (rho/R)*I010 + 0.5d0*rho*(I122-I100)
C
      PM   = 0.5d0*R*(I021+I001)
      P010 = I010
      P011 = I011
      RETURN
      END
C=======================================================================
C  Integration of FUNC over 0 <= x <= X2 with (near-)singular point at
C  x = zc: subintervals growing geometrically away from zc, Romberg
C  (open midpoint rule, qromo) on each of them.  On the two subintervals
C  adjacent to zc the substitution x = zc -/+ d*t**2 removes the
C  remaining (x-zc)*log|x-zc| kink of the subtracted integrand.
C=======================================================================
      SUBROUTINE QSING(func,zc,res)
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=201)
      INTEGER NN
      REAL*8 ZN(NMAX),RN(NMAX),GN(NMAX),ZK(0:NMAX),RK(0:NMAX),
     *       R2K(0:NMAX)
      REAL*8 X2,GINF,RINF,func,zc,res,d,D0,xl,xr,s
      REAL*8 zc0,dd,sgn,FUNCT
      COMMON/SHEET/X2,GINF,RINF,ZN,RN,GN,ZK,RK,R2K,NN
      COMMON/QTR/zc0,dd,sgn
      EXTERNAL func,midpnt,FUNCT
      PARAMETER (D0 = 0.05d0)
      res = 0.d0
      zc0 = zc
C  left of zc
      d  = D0
      xr = zc
      IF(xr.GT.0.d0) THEN
         xl = DMAX1(xr-d,0.d0)
         dd  = xr-xl
         sgn = -1.d0
         CALL qromo(FUNCT,0.d0,1.d0,s,midpnt)
         res = res + s
         xr = xl
         d  = 4.d0*d
      END IF
    1 IF(xr.GT.0.d0) THEN
         xl = DMAX1(xr-d,0.d0)
         CALL qromo(func,xl,xr,s,midpnt)
         res = res + s
         xr = xl
         d  = 4.d0*d
         GO TO 1
      END IF
C  right of zc
      d  = D0
      xl = zc
      IF(xl.LT.X2) THEN
         xr = DMIN1(xl+d,X2)
         dd  = xr-xl
         sgn = 1.d0
         CALL qromo(FUNCT,0.d0,1.d0,s,midpnt)
         res = res + s
         xl = xr
         d  = 4.d0*d
      END IF
    2 IF(xl.LT.X2) THEN
         xr = DMIN1(xl+d,X2)
         CALL qromo(func,xl,xr,s,midpnt)
         res = res + s
         xl = xr
         d  = 4.d0*d
         GO TO 2
      END IF
      RETURN
      END
C=======================================================================
      FUNCTION FUNCT(t)
      IMPLICIT NONE
      REAL*8 FUNCT,t,zc0,dd,sgn,FUNC
      COMMON/QTR/zc0,dd,sgn
      EXTERNAL FUNC
      FUNCT = FUNC(zc0+sgn*dd*t*t)*2.d0*dd*t
      RETURN
      END
C=======================================================================
C  Slipstream radius R(z): natural cubic spline through (0,1) and the
C  nodes (ZN(i),RN(i)); constant continuation outside [0,X2].
C=======================================================================
      SUBROUTINE SETSPL
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=201)
      INTEGER NN,I
      REAL*8 ZN(NMAX),RN(NMAX),GN(NMAX),ZK(0:NMAX),RK(0:NMAX),
     *       R2K(0:NMAX)
      REAL*8 X2,GINF,RINF
      COMMON/SHEET/X2,GINF,RINF,ZN,RN,GN,ZK,RK,R2K,NN
      ZK(0) = 0.d0
      RK(0) = 1.d0
      DO 1 I = 1,NN
         ZK(I) = ZN(I)
         RK(I) = RN(I)
    1 CONTINUE
      CALL spline(ZK(0),RK(0),NN+1,1.d31,1.d31,R2K(0))
      RINF = RN(NN)
      RETURN
      END
C=======================================================================
      FUNCTION RADIUS(x)
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=201)
      INTEGER NN
      REAL*8 ZN(NMAX),RN(NMAX),GN(NMAX),ZK(0:NMAX),RK(0:NMAX),
     *       R2K(0:NMAX)
      REAL*8 X2,GINF,RINF,RADIUS,x,y,dy
      COMMON/SHEET/X2,GINF,RINF,ZN,RN,GN,ZK,RK,R2K,NN
      IF(x.LE.0.d0) THEN
         RADIUS = 1.d0
      ELSE IF(x.GE.X2) THEN
         RADIUS = RN(NN)
      ELSE
         CALL splint(ZK(0),RK(0),R2K(0),NN+1,x,y,dy)
         RADIUS = y
      END IF
      RETURN
      END
C=======================================================================
      FUNCTION DRADIUS(x)
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=201)
      INTEGER NN
      REAL*8 ZN(NMAX),RN(NMAX),GN(NMAX),ZK(0:NMAX),RK(0:NMAX),
     *       R2K(0:NMAX)
      REAL*8 X2,GINF,RINF,DRADIUS,x,y,dy
      COMMON/SHEET/X2,GINF,RINF,ZN,RN,GN,ZK,RK,R2K,NN
      IF(x.LE.0.d0.OR.x.GE.X2) THEN
         DRADIUS = 0.d0
      ELSE
         CALL splint(ZK(0),RK(0),R2K(0),NN+1,x,y,dy)
         DRADIUS = dy
      END IF
      RETURN
      END
C=======================================================================
C  Sheet strength gamma(z): piecewise linear between the nodes,
C  constant outside
C=======================================================================
      FUNCTION GAMSH(x)
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=201)
      INTEGER NN,klo,khi,k
      REAL*8 ZN(NMAX),RN(NMAX),GN(NMAX),ZK(0:NMAX),RK(0:NMAX),
     *       R2K(0:NMAX)
      REAL*8 X2,GINF,RINF,GAMSH,x
      COMMON/SHEET/X2,GINF,RINF,ZN,RN,GN,ZK,RK,R2K,NN
      IF(x.LE.ZN(1)) THEN
         GAMSH = GN(1)
      ELSE IF(x.GE.ZN(NN)) THEN
         GAMSH = GN(NN)
      ELSE
         klo = 1
         khi = NN
    1    IF(khi-klo.GT.1) THEN
            k = (khi+klo)/2
            IF(ZN(k).GT.x) THEN
               khi = k
            ELSE
               klo = k
            END IF
            GO TO 1
         END IF
         GAMSH = GN(klo) + (GN(khi)-GN(klo))*(x-ZN(klo))
     *                     /(ZN(khi)-ZN(klo))
      END IF
      RETURN
      END
C=======================================================================
C  Cubic spline (Numerical Recipes), splint extended by the derivative
C=======================================================================
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      IMPLICIT NONE
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C=======================================================================
      SUBROUTINE splint(xa,ya,y2a,n,x,y,dy)
      IMPLICIT NONE
      INTEGER n
      REAL*8 x,y,dy,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) write(*,*) 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.d0
      dy=(ya(khi)-ya(klo))/h
     *   +(-(3.d0*a*a-1.d0)*y2a(klo)+(3.d0*b*b-1.d0)*y2a(khi))*h/6.d0
      return
      END
C=======================================================================
C  Elliptic integrals (Carlson forms, Numerical Recipes) as in
C  Conway-main.f
C=======================================================================
      FUNCTION ellf(phi,ak)
      IMPLICIT NONE
      REAL*8 ellf,ak,phi
CU    USES rf
      REAL*8 s,rf
      s=dsin(phi)
      ellf=s*rf(dcos(phi)**2,(1.d0-s*ak)*(1.d0+s*ak),1.d0)
      return
      END
C=======================================================================
      FUNCTION elle(phi,ak)
      IMPLICIT NONE
      REAL*8 elle,ak,phi
CU    USES rd,rf
      REAL*8 cc,q,s,rd,rf
      s=dsin(phi)
      cc=dcos(phi)**2
      q=(1.d0-s*ak)*(1.d0+s*ak)
      elle=s*(rf(cc,q,1.d0)-((s*ak)**2)*rd(cc,q,1.d0)/3.d0)
      return
      END
C=======================================================================
      FUNCTION rd(x,y,z)
      IMPLICIT NONE
      REAL*8 rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
      PARAMETER (ERRTOL=.0015d0,TINY=1.d-25,BIG=4.5d21,C1=3.d0/14.d0,
     *C2=1.d0/6.d0,C3=9.d0/22.d0,C4=3.d0/26.d0,C5=.25d0*C3,C6=1.5d0*C4)
      REAL*8 alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
     *sqrtz,sum,xt,yt,zt
      if(min(x,y).lt.0.d0.or.min(x+y,z).lt.TINY.or.max(x,y,
     *z).gt.BIG)write(*,*) 'invalid arguments in rd'
      xt=x
      yt=y
      zt=z
      sum=0.d0
      fac=1.d0
1     continue
        sqrtx=dsqrt(xt)
        sqrty=dsqrt(yt)
        sqrtz=dsqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=.25d0*fac
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        ave=.2d0*(xt+yt+3.d0*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6.d0*eb
      ee=ed+ec+ec
      rd=3.d0*sum+fac*(1.d0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*
     *(-C3*ec+delz*C4*ea)))/(ave*dsqrt(ave))
      return
      END
C=======================================================================
      FUNCTION rf(x,y,z)
      IMPLICIT NONE
      REAL*8 rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.0025d0,TINY=1.5d-38,BIG=3.d37,THIRD=1.d0/3.d0,
     *C1=1.d0/24.d0,C2=.1d0,C3=3.d0/44.d0,C4=1.d0/14.d0)
      REAL*8 alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0.d0.or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *z).gt.BIG)write(*,*) 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=dsqrt(xt)
        sqrty=dsqrt(yt)
        sqrtz=dsqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(dabs(delx),dabs(dely),dabs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/dsqrt(ave)
      return
      END
C=======================================================================
C  Romberg integration on an open interval (Numerical Recipes) as in
C  Conway-main.f
C=======================================================================
      SUBROUTINE qromo(func,a,b,ss,choose)
      IMPLICIT NONE
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func,choose
      PARAMETER (EPS=1.d-8, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.d0
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (dabs(dss).le.EPS*max(1.d0,dabs(ss))) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.d0
11    continue
      if (dabs(dss).gt.1.d-5*max(1.d0,dabs(ss)))
     *   WRITE(6,*)'too many steps in qromo',a,b,ss,dss
      END
C=======================================================================
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      IMPLICIT NONE
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=dabs(x-xa(1))
      do 11 i=1,n
        dift=dabs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0)write(*,*) 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C=======================================================================
      SUBROUTINE midpnt(func,a,b,s,n)
      IMPLICIT NONE
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5d0*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.d0*tnm)
        ddel=del+del
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.d0
      endif
      return
      END
