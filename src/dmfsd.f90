

subroutine dmfsd(a,n,eps,ier)
  !
  !   FACTORISATION DE CHOLESKY D'UNE MATRICE SDP
  !   MATRICE = TRANSPOSEE(T)*T
  !   ENTREE : TABLEAU A CONTENANT LA PARTIE SUPERIEURE STOCKEE COLONNE
  !            PAR COLONNE DE LA METRICE A FACTORISER
  !   SORTIE : A CONTIENT LA PARTIE SUPPERIEURE DE LA MATRICE triangulaire T
  !
  !   SUBROUTINE APPELE PAR DSINV
  !
  !   N : DIM. MATRICE
  !   EPS : SEUIL DE TOLERANCE
  !   IER = 0 PAS D'ERREUR
  !   IER = -1 ERREUR
  !   IER = K COMPRIS ENTRE 1 ET N, WARNING, LE CALCUL CONTINUE
  !
  implicit none

  integer,intent(in)::n
  integer,intent(inout)::ier
  double precision,intent(in)::eps 
  double precision,dimension(n*(n+1)/2),intent(inout)::A
  double precision :: dpiv,dsum,tol
  integer::i,k,l,kpiv,ind,lend,lanf,lind

  !
  !   TEST ON WRONG INPUT PARAMETER N
  !
  dpiv=0.d0
  if (n-1.lt.0) goto 12
  if (n-1.ge.0) ier=0
  !
  !   INITIALIZE DIAGONAL-LOOP
  !
  kpiv=0
  do k=1,n
     kpiv=kpiv+k
     ind=kpiv
     lend=k-1
     !
     !   CALCULATE TOLERANCE
     !
     tol=dabs(eps*sngl(A(kpiv)))
     !
     !   START FACTORIZATION-LOOP OVER K-TH ROW
     !
     do i=k,n
        dsum=0.d0
        if (lend.lt.0) goto 2
        if (lend.eq.0) goto 4
        if (lend.gt.0) goto 2
        !
        !   START INNER LOOP
        !
2       do l=1,lend
           lanf=kpiv-l
           lind=ind-l
           dsum=dsum+A(lanf)*A(lind)
        end do

        !     
        !   END OF INNEF LOOP
        !
        !   TRANSFORM ELEMENT A(IND)
        ! 	
4       dsum=A(ind)-dsum
        if (i-k.ne.0) goto 10
        if (i-k.eq.0) goto 5
        !   TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE
        !	


5       if (sngl(dsum)-tol.le.0) goto 6
        if (sngl(dsum)-tol.gt.0) goto 9
6       if (dsum.le.0) goto 12 
        if (dsum.gt.0) goto 7
7       if (ier.le.0) goto 8
        if (ier.gt.0) goto 9
8       ier=k-1
        !
        !   COMPUTE PIVOT ELEMENT
        !
9       dpiv=dsqrt(dsum)
        A(kpiv)=dpiv
        dpiv=1.D0/dpiv
        goto 11
        !
        !   CALCULATE TERMS IN ROW
        !
10      A(ind)=dsum*dpiv
11      ind=ind+i
     end do
  end do

  !
  !   END OF DIAGONAL-LOOP
  !
  return
12 ier=-1
  return

end subroutine dmfsd



