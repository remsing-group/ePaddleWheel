Program TetraRotorFunctions
Implicit None
  Integer, Parameter                   :: nframes = 80000
  Integer, Parameter                   :: natoms = 54
  Real*8, allocatable                  :: m2(:,:)
  Real*8, allocatable                  :: corm2(:,:), nori(:,:)
  Real*8, allocatable                  :: cornm2(:,:)
  Real*8, allocatable                  :: m2m20(:)
  Real*8, allocatable                  :: avgm2(:)
  Real*8, allocatable                  :: avgTCF2(:), time(:)
  Character*4                          :: dum1, dum2, mol
Real*8                                 :: tttt, s2
Integer                                :: ii, jj, kk, i, k, j, TCFcnt

!!!!!!!!!!!Allocate***************************************************************************************
allocate (m2(natoms,nframes))
allocate(corm2(natoms,nframes)); allocate(nori(natoms,nframes))
allocate(cornm2(natoms,nframes))
allocate (m2m20(natoms))
allocate (avgm2(natoms))
allocate (time(nframes))
allocate (avgTCF2(nframes))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Open(10,file = 'M-all.dat')

!----------------read Rotor Functions--------------------------------------

avgTCF2(:)=0.0
m2(:,:)=0.0

write(*,*) natoms

do ii = 1, nframes
   do k = 1, natoms
      Read(10,*)  m2(k,ii)
   enddo
   time(ii) = (ii-1)*0.0005
enddo

!----------------compute TCF for each I-atom--------------------------------------

avgm2(:)=0.0
m2m20(:)=0.0

do i = 1, natoms
  do j = 1, nframes

		corm2(i,j)=0.0; nori(i,j)=0.0
		m2m20(i) = m2m20(i) + m2(i,j)*m2(i,j)
		avgm2(i) = avgm2(i) + m2(i,j)
  
	enddo

	m2m20(i)=m2m20(i)/nframes
	avgm2(i)= avgm2(i)/nframes

	write(*,*) avgm2(i)

!-------------Compute the correlation functions-----------------------
  
  do j = 1, nframes
    do k = 1, nframes-j+1  
      corm2(i,k) = corm2(i,k) + (m2(i,j))*(m2(i,j+k-1))
			nori(i,k) = nori(i,k) + 1.0
		enddo
	enddo

  do j = 1, nframes 

    corm2(i,j) = corm2(i,j) / nori(i,j)
    cornm2(i,j) = corm2(i,j) / corm2(i,1)

  enddo

!--------------Print output-------------------------
  
  do j= 1, nframes 
		tttt = time(j) - time(1)
    avgTCF2(j) = avgTCF2(j) + corm2(i,j)
	enddo
	
 TCFcnt = TCFcnt + 1.0

enddo

!--------Average the TCFs over all I- to get the averaged TCF-------------

open(60,file="AvgMTCFs.dat")

avgTCF2(:)=avgTCF2(:)/TCFcnt

do j=1, nframes
  tttt = time(j) - time(1)

	If(j==1) Then
		s2=avgTCF2(1)
  end If

  avgTCF2(j)=avgTCF2(j)/s2

	write(60,*) tttt, avgTCF2(j)
enddo

end
