Program TetraRotorFunctions
Implicit None
  Integer, Parameter                   :: numAg = 54,nframes = 80000
  Integer                              :: listAg(nframes,numAg,10)
  Integer                              :: res1, res2
  Real*8, allocatable                  :: corhAB(:,:), nori(:,:), cornhAB(:,:)
  Real*8, allocatable                  :: avgTCF(:), time(:)

  Real*8                               :: tttt, s1
  Integer                              :: i, j, k, l, kmax, TCFcnt, ngap, nrun

  character(len=10)                    :: file_id
  character(len=50)                    :: file_name
  Character*4                          :: dum1, dum2, mol

!!!!!!!!!!!Allocate***************************************************************************************
allocate(corhAB(nframes,numAg)); allocate(cornhAB(nframes,numAg)); allocate(nori(nframes,numAg))
allocate (time(nframes))
allocate (avgTCF(nframes))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ngap = 1
nrun = nframes

!----------------read coordinates--------------------------------------
Open(10,file = 'complete_neighbour_list.dat')

listAg(:,:,:)=0.0

do k = 1, nframes
    time(k) = (k-1)*0.0005
  do i = 1, numAg
    Read(10,*) listAg(k,i,:)
  enddo
    !write(20,*) listAg(k,1,:)
enddo

close(10)

Write(*,*) "Data Reading Complete"

res1 = compare_arrays(listAg(1000,2,:), listAg(2,1,:))

!-------------Compute the correlation functions-----------------------

avgTCF(:)=0.0
corhAB(:,:)=0.0; nori(:,:)=0.0

do i = 1, numAg
  write(*,*) "Ag ion number = ", i
  do j = 1, nframes, ngap
    kmax = min(nframes, j+nrun)
    do k = j, kmax
      res1 = compare_arrays(listAg(j,i,:), listAg(j,i,:))
      res2 = compare_arrays(listAg(j,i,:), listAg(k,i,:))
		  corhAB(k-j+1,i) = corhAB(k-j+1,i) + res1*res2
		  nori(k-j+1,i) = nori(k-j+1,i) + 1
	  enddo
  enddo

  do k = 1, nframes
		corhAB(k,i) = corhAB(k,i) / nori(k,i)
    cornhAB(k,i) = corhAB(k,i) / corhAB(1,i)
  enddo

  do k= 1, nframes   !nconf-1
		tttt = time(k) - time(1)
		avgTCF(k) = avgTCF(k) + corhAB(k,i)    ! Sum for all Ag ions
	enddo
	TCFcnt = TCFcnt + 1.0
enddo

!--------Average the TCFs over all Ag+ to get the averaged TCF-------------

open(60,file="AvgRTCFs-i.dat")
avgTCF(:)=avgTCF(:)/TCFcnt
do k=1, nframes
  tttt = time(k) - time(1)
	If(k==1) Then
		s1=avgTCF(1)
  end If
	avgTCF(k)=avgTCF(k)/s1
	write(60,*) tttt, avgTCF(k)
enddo

CONTAINS

FUNCTION compare_arrays(a, b) RESULT(res)
  Implicit None
    INTEGER, INTENT(IN) :: a(:), b(:)
    INTEGER :: res
    IF (ALL(a == b)) THEN
        res = 1
    ELSE
        res = 0
    END IF
END FUNCTION compare_arrays

END
