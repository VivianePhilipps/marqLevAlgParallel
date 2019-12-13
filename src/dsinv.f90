

!------------------------------------------------------------
!                            DSINV
!------------------------------------------------------------


subroutine dsinv(A,N,EPS,IER,DET)

  !
  !     INVERSION D'UNE MATRICE SYMETRIQUE DEFINIE POSITIVE :
  !
  !     MATRICE = TRANSPOSEE(T)*T
  !     INERSE(MATRICE) = INVERSE(T)*INVERSE(TRANSPOSEE(T))
  !
  !     A : TABLEAU CONTENANT LA PARTIE SUPERIEURE DE LA MATRICE A INVERSER
  !         STOCKEE COLONNE PAR COLONNE
  !     DIM. MATRICE A INVERSER = N
  !     DIM. TABLEAU A = N*(N+1)/2
  !
  !     EPS : SEUIL DE TOLERANCE AU-DESSOUS DUQUEL UN PIVOT EST CONSIDERE
  !           COMME NUL
  !
  !     IER : CODE D'ERREUR
  !         IER=0 PAS D'ERREUR
  !         IER=-1 ERREUR SUR LA DIM.N OU MATRICE PAS DEFINIE POSITIVE
  !         IER=1 PERTE DE SIGNIFICANCE, LE CALCUL CONTINUE
  !
  implicit none

  integer,intent(in)::n
  integer,intent(inout)::ier
  double precision,intent(inout)::eps      
  double precision,intent(inout),optional::det     
  double precision,dimension(n*(n+1)/2),intent(inout)::A     
  double precision::din,work
  integer::ind,ipiv,i,j,k,l,min,kend,lhor,lver,lanf

  !
  !     FACTORIZE GIVEN MATRIX BY MEANS OF SUBROUTINE DMFSD
  !     A=TRANSPOSE(T) * T
  !

  call dmfsd(A,n,eps,ier)

  det=0.d0

  if (ier.lt.0) goto 9
  if (ier.ge.0) det=0.d0
  !
  !     INVERT UPPER TRIANGULAR MATRIX T
  !     PREPARE INVERSION-LOOP
  !
  !
  ! calcul du log du determinant    

  do i=1,n
     det=det+dlog(A(i*(i+1)/2))
  end do
  det=2*det
  ipiv=n*(n+1)/2
  ind=ipiv
  !
  !     INITIALIZE INVERSION-LOOP
  !
  do i=1,n
     din=1.d0/A(ipiv)
     A(ipiv)=din
     min=n
     kend=i-1
     lanf=n-kend
     if (kend.le.0) goto 5
     if (kend.gt.0) j=ind
     !
     !     INITIALIZE ROW-LOOP
     !
     do k=1,kend
        work=0.d0
        min=min-1
        lhor=ipiv
        lver=j
        !
        !     START INNER LOOP
        !
        do l=lanf,min 
           lver=lver+1
           lhor=lhor+l
           work=work+A(lver)*A(lhor)
        end do
        !
        !     END OF INNER LOOP
        !
        A(j)=-work*din
        j=j-min
     end do

     !
     !     END OF ROW-LOOP
     !
5    ipiv=ipiv-min 
     ind=ind-1
  end do

  !
  !     END OF INVERSION-LOOP
  !
  !     CALCULATE INVERSE(A) BY MEANS OF INVERSE(T)
  !     INVERSE(A) = INVERSE(T) * TRANSPOSE(INVERSE(T))
  !     INITIALIZE MULTIPLICATION-LOOP
  !
  do i=1,n
     ipiv=ipiv+i
     j=ipiv
     !
     !     INITIALIZE ROW-LOOP
     !
     do k=i,n
        work=0.d0
        lhor=j
        !
        !     START INNER LOOP
        !
        do l=k,n
           lver=lhor+k-i
           work=work+A(lhor)*A(lver)
           lhor=lhor+l
        end do
        !
        !     END OF INNER LOOP
        !       
        A(j)=work
        j=j+k
     end do
  end do

  !
  !     END OF ROW-AND MULTIPLICATION-LOOP
  !
9 return
end subroutine dsinv


