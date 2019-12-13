


!------------------------------------------------------------
!                         DCHOLE
!------------------------------------------------------------

subroutine dchole(a,k,nq,idpos)

  implicit none

  integer,intent(in)::k,nq
  integer,intent(inout)::idpos
  double precision,dimension(k*(k+3)/2),intent(inout)::a

  integer::i,ii,i1,i2,i3,m,j,k2,jmk
  integer::ijm,irm,jji,jjj,l,jj,iil,jjl,il
  integer,dimension(k)::is
  double precision::term,xn,diag,p
  equivalence (term,xn)


  !      ss programme de resolution d'un systeme lineaire symetrique
  !
  !       k ordre du systeme /
  !       nq nombre de seconds membres
  !
  !       en sortie les seconds membres sont remplaces par les solutions
  !       correspondantes
  !
  i2=0
  ii=0
  idpos=0
  k2=k+nq
  !     calcul des elements de la matrice
  do i=1,k   
     ii=i*(i+1)/2
     !       elements diagonaux
     diag=a(ii)
     i1=ii-i
     if(i-1.ne.0) goto 1
     if(i-1.eq.0) goto 4
1    i2=i-1
     do l=1,i2
        m=i1+l
        p=a(m)
        p=p*p
        if(is(l).lt.0) goto 2
        if(is(l).ge.0) goto 3
2       p=-p
3       diag=diag-p
     end do

4    if(diag.lt.0) goto 5
     if(diag.eq.0) goto 50
     if(diag.gt.0) goto 6
5    is(i)=-1
     idpos=idpos+1
     diag=-dsqrt(-diag)
     a(ii)=-diag
     goto 7
6    is(i)=1
     diag=dsqrt(diag)
     a(ii)=diag
     !       elements non diagonaux
7    i3=i+1
     do j=i3,k2
        jj=j*(j-1)/2+i
        jmk=j-k-1
        if(jmk.le.0) goto 9
        if(jmk.gt.0) goto 8
8       jj=jj-jmk*(jmk+1)/2
9       term=a(jj)
        if(i-1.ne.0) goto 10
        if(i-1.eq.0) goto 13 
10      do l=1,i2
           iil=ii-l
           jjl=jj-l
           p=a(iil)*a(jjl)
           il=i-l
           if(is(il).lt.0) goto 11
           if(is(il).ge.0) goto 12
11         p=-p
12         term=term-p
        end do
13      a(jj)=term/diag
     end do
  end do

  !       calcul des solutions
  jj=ii-k+1
  do l=1,nq
     jj=jj+k
     i=k-1
14   jji=jj+i
     xn=a(jji)
     if(i-k+1.lt.0) goto 20
     if(i-k+1.ge.0) goto 22
20   j=k-1
21   jjj=jj+j
     ijm=i+1+j*(j+1)/2
     xn=xn-a(jjj)*a(ijm)
     if(j-i-1.le.0) goto 22
     if(j-i-1.gt.0) goto 30
30   j=j-1
     goto 21
22   irm=(i+1)*(i+2)/2
     a(jji)=xn/a(irm)
     if(i.le.0) cycle
     if(i.gt.0) goto 40
40   i=i-1
     go to 14
  end do
50 continue
  return
end subroutine dchole


