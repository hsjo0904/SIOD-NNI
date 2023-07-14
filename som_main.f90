PROGRAM som_main
!******************************************************************
!
! Self-Organizing Feature Map (SOFM) Analysis Main
!
!   INFORMS : mhkwon@kiost.ac.kr     8Jun2019
!   REVISED : 1st release            8Aug2019
!             ensemble init.        19Aug2019
!             2nd release           21Aug2019
! 
!******************************************************************
  IMPLICIT NONE
  REAL, ALLOCATABLE :: x(:,:), w(:,:,:), wm(:,:,:)
  INTEGER, ALLOCATABLE ::  clust1(:), clust2(:)
  CHARACTER(LEN=200) :: ifile, wfile, cfile
  INTEGER :: nd, nt, m1, m2, nter, nens
  INTEGER :: iens, it
  REAL :: alpha0, taua, sigma0, taus, dw_tol
  REAL :: qerror, qerrorm
  NAMELIST /files/ ifile, wfile, cfile, nd, nt, m1, m2
  NAMELIST /somparameters/ nter, nens, alpha0, taua, &
                           sigma0, taus, dw_tol

  WRITE(*,'(A)') ''
  WRITE(*,'(A)') '> Self-Organizing Feature Map Analysis'
  WRITE(*,'(A)') '>   Last Updated : 19Aug2019'
  WRITE(*,'(A)') ''

! namelists & array allocations
  OPEN(10,FILE='som.nml',STATUS='OLD')
  READ(10,nml=files)
  READ(10,nml=somparameters)
  CLOSE(10)
  ALLOCATE(x(nd,nt),w(nd,m1,m2),wm(nd,m1,m2),clust1(nt),clust2(nt))

! input data
  OPEN(11,FILE=TRIM(ifile),RECL=nd*nt*4, &
  ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD')
  READ(11,REC=1) x
  CLOSE(11)

  WRITE(*,'(A)') '> input file = '//TRIM(ifile)
  WRITE(*,'(A,I5,A,I5)') '>   nd = ',nd,', nt = ',nt
  WRITE(*,'(A)') '> weight-vector output file = '//TRIM(wfile)
  WRITE(*,'(A,I5,A,I3,A,I3)') '>   nd = ',nd,', m1 = ',m1,', m2 = ',m2
  WRITE(*,'(A)') '> cluster output file = '//TRIM(cfile)
  WRITE(*,'(A,I5,A,I3,A,I3)') '>   nd = ',nd,', m1 = ',m1,', m2 = ',m2
  WRITE(*,'(A)') ''

  WRITE(*,'(A,I8)')    '> max. # of iterations  = ',nter
  WRITE(*,'(A,I8)')    '> number of ensembles   = ',nens
  WRITE(*,'(A,F10.5)') '> init. learning rate   = ',alpha0
  WRITE(*,'(A,F10.5)') '> tau (learning rate)   = ',taua
  WRITE(*,'(A,F10.5)') '> init. influence radius= ',sigma0
  WRITE(*,'(A,F10.5)') '> tau (influence radius)= ',taus

! som w/ ens
  qerrorm=1E+8
  WRITE(*,'(A)',ADVANCE='NO') '> finding the weight vectors'
  DO iens=1,nens

    IF ( 10*iens/nens==10.*iens/nens ) THEN
      WRITE(*,'(A)',ADVANCE='NO') '.'
    ENDIF

    CALL som(x,nd,nt,w,m1,m2,nter,alpha0,taua,sigma0,taus,dw_tol)

    CALL cluster(x,w,nd,nt,m1,m2,clust1,clust2,qerror)

    IF ( qerror<qerrorm ) THEN
      qerrorm=qerror
      wm=w
    ENDIF

  ENDDO !iens
  WRITE(*,'(A)') ''
  WRITE(*,'(A,F15.10)') '> Minimum QE = ',qerrorm

! clustering for the minimum quatization error
  CALL cluster(x,wm,nd,nt,m1,m2,clust1,clust2,qerror)

! output weight vectors
  OPEN(21,FILE=TRIM(wfile),RECL=nd*m1*m2*4, &
  ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN')
  WRITE(21,REC=1) wm
  CLOSE(21)
  WRITE(*,'(A)') '> '//TRIM(wfile) &
    //' (output weight vectors) has been generated.'

! output clusters
  OPEN(22,FILE=TRIM(cfile),STATUS='UNKNOWN')
  WRITE(22,'(A)') 'samples  clust1  clust2'
  DO it=1,nt
    WRITE(22,'(I7,1X,I7,1X,I7)') it,clust1(it),clust2(it)
  ENDDO
  CLOSE(22)
  WRITE(*,'(A)') '> '//TRIM(cfile) &
    //' (output clusters) has been generated.'
  WRITE(*,'(A)') ''

END PROGRAM som_main
!******************************************************************

SUBROUTINE som(x,nd,nt,w,m1,m2,nter,alpha0,taua,sigma0,taus,dw_tol)
!******************************************************************
!
! Self-Organizing Feature Map (SOFM or SOM) Analysis 
!
!   INFORMS : mhkwon@kiost.ac.kr   8Jun2019
!   REVISED : 1st release          8Aug2019
!             dw                  19Aug2019
!   ACCEPTS : x(nd,nt),m1,m2,nter,alpha0,taua,sigma0,taus,dw_tol
!   RETURNS : w(nd,m1,m2)
!
!   USES normalize, initwt, bmu, edis, dw, cluster
!
!   SOM PARAMETERS
!    nter   : maximum number of iterations
!    alpha0 : initial learning rate
!    taua   : time constant for learning rate
!    sigma0 : initial influence radius
!    taus   : time constant for influence radius
!    dw_tol : tolerence value for the weight vector differences
!
!******************************************************************
  IMPLICIT NONE
  INTEGER :: nd, nt, m1, m2
  REAL, DIMENSION(nd,nt) :: x
  REAL, DIMENSION(nd,m1,m2) :: w, w0
  REAL, PARAMETER :: alphamin=1E-3 ! lower limit of alpha
  INTEGER :: nter
  REAL :: alpha, alpha0, taua, sigma, sigma0, taus
  REAL :: t, rnd, d, h, dw, dw_tol
  INTEGER :: iter, id, it, i, j, i0, j0

! data normalization
  DO id=1,nd
    CALL normalize(x(id,:),nt)
  ENDDO

! iterations
  CALL initwt(w,nd,m1,m2) ! random initial weights
  CALL RANDOM_SEED
  DO iter=1,nter
    t=iter-1. ! time for iterations
    IF ( iter<=nt ) THEN
      it=iter
    ELSE
      CALL RANDOM_NUMBER(rnd)
      it=INT(rnd*nt)+1 ! select a training sample at random
    ENDIF
    w0=w
    CALL bmu(w,x(:,it),nd,m1,m2,i0,j0) ! best matching unit
    alpha=MAX(alpha0*EXP(-t/taua),alphamin) ! learning rate
    sigma=MAX(sigma0*EXP(-t/taus),1.)       ! influence radius
    DO j=1,m2
    DO i=1,m1
      d=SQRT((i0-i*1.)**2+(j0-j*1.)**2)
      ! 2D Euclidean distance for rectangular grids
      IF ( d<sigma ) THEN
        h=EXP(-0.5*d**2/sigma**2) ! Gaussian neigborhood function
        w(:,i,j)=w(:,i,j)+alpha*h*(x(:,it)-w(:,i,j))
        ! update weight vectors within influence radius
      ENDIF
    ENDDO
    ENDDO
    IF ( dw(w0,w,nd,m1,m2)<dw_tol ) EXIT
  ENDDO !iter

END SUBROUTINE som
!************************************************************

SUBROUTINE normalize(x,n)
  IMPLICIT NONE
  INTEGER :: n
  REAL, DIMENSION(n) :: x
  INTEGER :: i
  REAL :: xm, xv
  xm=0.
  DO i=1,n
    xm=xm+x(i)/n
  ENDDO
  xv=0.
  DO i=1,n
    xv=xv+(x(i)-xm)**2/n
  ENDDO
  DO i=1,n
    x(i)=(x(i)-xm)/SQRT(xv)
  ENDDO
END SUBROUTINE normalize

SUBROUTINE initwt(w,nd,m1,m2)
! initialize weights
  IMPLICIT NONE
  INTEGER :: nd, m1, m2
  REAL, DIMENSION(nd,m1,m2) :: w
  INTEGER, DIMENSION(8) :: values
  INTEGER, DIMENSION(1000) :: seed
  INTEGER :: i, j, id
  REAL :: rnd
  CALL DATE_AND_TIME(values=values)
  seed(:)=values(8)
  CALL RANDOM_SEED(PUT=seed)
  DO j=1,m2
  DO i=1,m1
  DO id=1,nd
    CALL RANDOM_NUMBER(rnd)
    w(id,i,j)=2.*rnd-1.
  ENDDO
  ENDDO
  ENDDO
END SUBROUTINE initwt

SUBROUTINE bmu(w,x1,nd,m1,m2,i0,j0)
! find the Best Matching Unit (BMU)
  IMPLICIT NONE
  INTEGER :: nd, m1, m2
  REAL, DIMENSION(nd,m1,m2) :: w
  REAL, DIMENSION(nd) :: x1
  INTEGER :: i0, j0
  REAL :: dismin, dis, edis
  INTEGER :: i, j
  dismin=1E+20
  DO j=1,m2
  DO i=1,m1
    dis=edis(w(:,i,j),x1,nd)
    IF ( dis<dismin ) THEN
      dismin=dis
      i0=i
      j0=j
    ENDIF
  ENDDO
  ENDDO
END SUBROUTINE bmu

FUNCTION edis(x1,x2,nd)
! Euclidean distance squared devided by dim.
  IMPLICIT NONE
  REAL :: edis
  INTEGER :: nd
  REAL, DIMENSION(nd) :: x1, x2
  INTEGER :: id
  edis=0.
  DO id=1,nd
    edis=edis+(x1(id)-x2(id))**2/nd
  ENDDO
END FUNCTION edis

FUNCTION dw(w0,w,nd,m1,m2)
  IMPLICIT NONE
  REAL :: dw
  INTEGER :: nd, m1, m2
  REAL, DIMENSION(nd,m1,m2) :: w0, w
  REAL :: dis, edis
  INTEGER :: i, j
  dw=0.
  DO j=1,m2
  DO i=1,m1
    dis=edis(w0(:,i,j),w(:,i,j),nd)
    IF ( dis>dw ) dw=dis
  ENDDO
  ENDDO
END FUNCTION dw

SUBROUTINE cluster(x,w,nd,nt,m1,m2,clust1,clust2,qerror)
  IMPLICIT NONE
  INTEGER :: nd, nt, m1, m2
  REAL, DIMENSION(nd,nt) :: x
  REAL, DIMENSION(nd,m1,m2) :: w
  INTEGER, DIMENSION(nt) :: clust1, clust2
  REAL :: edis, dis, dismin, qerror
  INTEGER :: it, i, j, i0, j0
  qerror=0.
  DO it=1,nt
    dismin=1E+20
    DO j=1,m2
    DO i=1,m1
      dis=edis(w(:,i,j),x(:,it),nd)
      IF ( dis<dismin ) THEN
        i0=i
        j0=j
        dismin=dis
      ENDIF
    ENDDO
    ENDDO
    clust1(it)=i0
    clust2(it)=j0
    qerror=qerror+dismin/nt
  ENDDO !it
END SUBROUTINE cluster
