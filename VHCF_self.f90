Program dist
Implicit None
  Integer, Parameter                   :: natoms =108,nframes =1000
  Real*8, Parameter                    :: Lx=15.2565, Ly=15.2565, Lz=15.2565
  Real*8, allocatable                  :: xA(:,:), yA(:,:), zA(:,:)
  Real*8, allocatable                  :: xB(:,:), yB(:,:), zB(:,:)
  Character*4                          :: dum1,dum2
Integer, Parameter                     :: nmolAg = 54, nmolI = 54
Integer                                :: maxbinr, moltotal, nbinAA
Real*8                                 :: pi, delr, drAA, rAA, dxAA, dyAA, dzAA
Integer                                :: i, k, ii, kk, jj, binr, bint, delt
Real*8                                 :: numdensity, ntotalr, ru, rl, r, C
Real*8, Allocatable, Dimension (:,:)     :: histrAA, GrAA
Integer, Dimension (6)                 :: t_lag, countA

pi = 3.1415926536D0
delr = 0.05
maxbinr = int (Lx / (1 * delr)) ! 2 if truncating at L/2
moltotal = nmolAg

!!!!!!!!!!!Allocate***************************************************************************************
allocate (xA(nframes,nmolAg)); allocate (yA(nframes,nmolAg)); allocate (zA(nframes,nmolAg))
allocate (xB(nframes,nmolI)); allocate (yB(nframes,nmolI)); allocate (zB(nframes,nmolI))

allocate (histrAA (maxbinr,6)); allocate (GrAA (maxbinr,6))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Open(50,file = 'VHCF-self_Ag.dat')

!----------------read coordinates--------------------------------------

Open(10,file = '../AgI-MD-pos-750K.xyz')

xA(:,:)=0.0; yA(:,:)=0.0; zA(:,:)=0.0
xB(:,:)=0.0; yB(:,:)=0.0; zB(:,:)=0.0

do k = 1, nframes
   Read(10,*)
   Read(10,*)
   do ii = 1, natoms/2, 2
      Read(10,*)  dum1, xA(k,ii),yA(k,ii), zA(k,ii)
      Read(10,*)  dum1, xA(k,ii+1),yA(k,ii+1), zA(k,ii+1)
      Read(10,*)  dum2, xB(k,ii),yB(k,ii), zB(k,ii)
      Read(10,*)  dum2, xB(k,ii+1),yB(k,ii+1), zB(k,ii+1)
    !  write(20,'(a4,3F20.10)')  dum2, xB(k,ii+1),yB(k,ii+1), zB(k,ii+1)
   enddo
enddo
write(*,*) k
close(10)

Write(*,*) "coordinates Read"

histrAA(:,:) = 0
countA(:) = 0

t_lag = (/100, 250, 500, 1000, 2000, 4000 /)

!dzAA = dzAA - (nint(dzAA / Lz)) * Lz

do i = 1, 6
Do kk = 1, nframes-t_lag(i)
  jj = kk + t_lag(i)
  Do ii = 1, nmolAg
  countA(i) = countA(i) + 1
	dxAA = xA (kk, ii) - xA (jj, ii)
	dyAA = yA (kk, ii) - yA (jj, ii)
	dzAA = zA (kk, ii) - zA (jj, ii)
	drAA = (dxAA * dxAA) + (dyAA * dyAA) + (dzAA * dzAA)
	rAA = sqrt (drAA)
	binr = (int (rAA / delr)) + 1
	If (binr .LE. maxbinr) Then
	histrAA (binr,i) = histrAA (binr,i) + 1
 	End If

  End Do
End Do
End Do

write(*,*) countA

Write (50, *) 0. , 0. , 0. , 0., 0. , 0., 0. 

Do binr = 1, maxbinr
rl = real (binr - 1) * delr
ru = rl + delr
C = 4.0 * pi * (ru**3 - rl**3) / 3.0   
r = real (binr - 0.5) * delr

GrAA (binr,1) = 4*pi*r*r * real (histrAA (binr,1)) / (real (countA(1)*C))
GrAA (binr,2) = 4*pi*r*r * real (histrAA (binr,2)) / (real (countA(2)*C))
GrAA (binr,3) = 4*pi*r*r * real (histrAA (binr,3)) / (real (countA(3)*C))
GrAA (binr,4) = 4*pi*r*r * real (histrAA (binr,4)) / (real (countA(4)*C))
GrAA (binr,5) = 4*pi*r*r * real (histrAA (binr,5)) / (real (countA(5)*C))
GrAA (binr,6) = 4*pi*r*r * real (histrAA (binr,6)) / (real (countA(6)*C))

Write (50, *) r, GrAA (binr,1), GrAA (binr,2), GrAA (binr,3), GrAA (binr,4), GrAA (binr,5), GrAA (binr,6)

End Do
End
