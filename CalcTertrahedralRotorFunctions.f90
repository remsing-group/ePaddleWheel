Program TetraRotorFunctions
Implicit None
  Integer, Parameter                   :: nsite = 594,nframes = 80000
  Integer, Parameter                   :: natoms = 108
  Real*8, Parameter                    :: conv=0.529177 ! Bohr to Angstrom
  Real*8, Parameter                    :: cut=2.35, cutl=1.2 ! in Angstrom
  Real*8, Parameter                    :: Lx=0.2883061E+02*conv, Ly=0.2883061E+02*conv, Lz=0.2883061E+02*conv
  Real*8, allocatable                  :: xA(:,:), yA(:,:), zA(:,:)
  Real*8, allocatable                  :: xWC(:,:), yWC(:,:), zWC(:,:)
  Real*8, allocatable                  :: newx(:,:), newy(:,:), newz(:,:)
  Character*4                          :: dum1, dum2, mol
Real*8                                 :: drAWC, rAWC, dxAWC, dyAWC, dzAWC, xshift, yshift, zshift
Real*8                                 :: M2, time, xjk, yjk, zjk, rjk
Integer                                :: ii, jj, kk, i, k, j, lc, ts, count, lcnt

!!!!!!!!!!!Allocate***************************************************************************************
allocate (xA(nframes,natoms)); allocate (yA(nframes,natoms)); allocate (zA(nframes,natoms))
allocate (xWC(nframes,nsite-natoms)); allocate (yWC(nframes,nsite-natoms)); allocate (zWC(nframes,nsite-natoms))
allocate (newx(nframes,nsite)); allocate (newy(nframes,nsite)); allocate (newz(nframes,nsite))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Open(10,file = 'TotalWannier.xyz')
Open(50,file = 'M-all.dat')

!----------------read coordinates--------------------------------------

xA(:,:)=0.0; yA(:,:)=0.0; zA(:,:)=0.0
xWC(:,:)=0.0; yWC(:,:)=0.0; zWC(:,:)=0.0

do k = 1, nframes
   Read(10,*)
   Read(10,*)
   do ii = 1, natoms
      Read(10,*)  dum1, xA(k,ii),yA(k,ii), zA(k,ii)
   enddo
   do jj = 1, nsite-natoms
      Read(10,*)  dum1, xWC(k,jj),yWC(k,jj), zWC(k,jj)
   enddo
enddo

!----------------search for I WCs--------------------------------------

Write(*,*) "I-MLWFC bond search"

do k = 1, nframes
  Write(*,*) "timestep:", k
  count = 1
  do i = 3, natoms, 4
  do j = 0, 1
  ii = i + j
  lcnt = 0
  do jj = 1, nsite-natoms

    dxAWC = xWC (k, jj) - xA (k, ii)
    dxAWC = dxAWC - (nint(dxAWC / Lx)) * Lx;    xshift = (nint(dxAWC / Lx)) * Lx
    dyAWC = yWC (k, jj) - yA (k, ii)
    dyAWC = dyAWC - (nint(dyAWC / Ly)) * Ly;    yshift = (nint(dyAWC / Ly)) * Ly
    dzAWC = zWC (k, jj) - zA (k, ii)
    dzAWC = dzAWC - (nint(dzAWC / Lz)) * Lz;    zshift = (nint(dzAWC / Lz)) * Lz

    drAWC = (dxAWC * dxAWC) + (dyAWC * dyAWC) + (dzAWC * dzAWC)
    rAWC = sqrt (drAWC)
    
    If (rAWC.le.cutl) Then

      If (lcnt.eq.0) Then

      newx(k, count) = xA (k, ii)
      newy(k, count) = yA (k, ii)
      newz(k, count) = zA (k, ii)
      count = count+1

      End If

      newx(k, count) = xWC (k, jj) + xshift
      newy(k, count) = yWC (k, jj) + yshift
      newz(k, count) = zWC (k, jj) + zshift
      count=count+1
      lcnt=lcnt+1

    End IF
  enddo
  If (lcnt.ne.4) Then
  End If
  enddo
  enddo
enddo

!----------------Calculate Tetrahedral Rotor Functions--------------------------------

do i = 1, nframes
  do j = 1, count-1, 5
    M2=0.0

    do k =1,4

      xjk = newx(i, j+k) - newx(i, j)
      xjk = xjk - (nint(xjk / Lx)) * Lx
      yjk = newy(i, j+k) - newy(i, j)
      yjk = yjk - (nint(yjk / Ly)) * Ly
      zjk = newz(i, j+k) - newz(i, j)
      zjk = zjk - (nint(zjk / Lz)) * Lz

      rjk = sqrt(xjk**2 + yjk**2 + zjk**2)

      xjk = xjk/rjk
      yjk = yjk/rjk
      zjk = zjk/rjk

      M2 = M2 + (5.0*xjk**3 - 3.0*xjk)

    enddo

    M2 = M2 * (3.0/40.0)*sqrt(5.0)
    
    write (50,*) M2

  enddo
enddo

END
