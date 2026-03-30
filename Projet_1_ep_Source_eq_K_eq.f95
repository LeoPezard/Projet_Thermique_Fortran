Program batterie_multicouche
Implicit None

!========================
! Déclarations
!========================
Integer :: N, i, j, Ncouches, k, sumN
Integer, Parameter :: N_rep = 50 ! Nombre de couches dans la batterie 
Integer, Allocatable :: Nvol(:)

Real*8 :: L_tot
Real*8, Allocatable :: K_couche(:), S_couche(:), e_couche(:), xcou(:)
Real*8, Allocatable :: K_vol(:), S_vol(:)

Real*8 :: hw, T_ext_w, qw, he, T_ext_e, qe, alpha
Real*8 :: R_tot_cell, K_eq_cell, L_cell, S_eq

Real*8, Allocatable :: T(:), ap(:), ae(:), aw(:), b(:)
Real*8, Allocatable :: X(:), dw(:), de(:), dx(:)
Real*8, Allocatable :: P(:), Q(:)

Real*8 :: denom, ke, kw

!========================
! Lecture données
!========================
Open(15,file="data_lente.txt")

Read(15,*) 
Read(15,*) N
Read(15,*) Ncouches

Allocate(e_couche(Ncouches), K_couche(Ncouches), S_couche(Ncouches))
Allocate(Nvol(Ncouches), xcou(Ncouches))


Do j = 1, Ncouches
    Read(15,*) e_couche(j), K_couche(j), S_couche(j)
End Do

Read(15,*) hw
Read(15,*) T_ext_w
Read(15,*) qw
Read(15,*) he
Read(15,*) T_ext_e
Read(15,*) qe
Read(15,*) alpha

Close(15)

! Conversion en Kelvin
T_ext_w = T_ext_w + 273.15d0
T_ext_e = T_ext_e + 273.15d0

!========================
! Épaisseur totale
!========================
L_tot = sum(e_couche)

!========================
! Répartition des volumes (CORRIGÉ)
!========================
sumN = 0
Do j = 1, Ncouches-1
    Nvol(j) = Max(1, Nint(dble(N)*e_couche(j)/L_tot))
    sumN = sumN + Nvol(j)
End Do
Nvol(Ncouches) = N - sumN

!========================
! Allocation
!========================
Allocate(X(N), T(N), ap(N), aw(N), ae(N), b(N))
Allocate(dw(N), de(N), dx(N), K_vol(N), S_vol(N))
Allocate(P(N), Q(N))

!========================
! Maillage
!========================
i = 0
Do j = 1, Ncouches
    Block
        Real*8 :: dx_j
        dx_j = e_couche(j)/dble(Nvol(j))
        Do k = 1, Nvol(j)
            i = i + 1
            If (i == 1) Then
                X(i) = dx_j/2.d0
            Else
                X(i) = X(i-1) + dx_j
            End If
        End Do
    End Block
End Do

!========================
! Interfaces couches
!========================
xcou(1) = e_couche(1)
Do j = 2, Ncouches
    xcou(j) = xcou(j-1) + e_couche(j)
End Do

!========================
! Propriétés locales (ROBUSTE)
!========================

R_tot_cell = 0.d0
Do j = 1, Ncouches
    R_tot_cell = R_tot_cell + e_couche(j) / K_couche(j)
End Do

L_cell = L_tot / dble(N_rep)   ! épaisseur d'une cellule (1 répétition sur N_rep)
K_eq_cell = L_cell / (R_tot_cell / dble(N_rep))
! Simplifié : K_eq_cell = L_tot / R_tot_cell  (indépendant de N_rep)
K_eq_cell = L_tot / R_tot_cell

S_eq = 0.d0
Do j = 1, Ncouches
    S_eq = S_eq + S_couche(j) * e_couche(j)
End Do
S_eq = S_eq / L_tot

Do i = 1, N
    k = 1
    Do while (X(i) > xcou(k) .and. k < Ncouches)
        k = k + 1
    End do
    K_vol(i) = K_eq_cell
    S_vol(i) = S_eq
End Do



!========================
! Distances
!========================
dw(1) = X(1)
de(1) = X(2) - X(1)
dx(1) = (dw(1)+de(1))/2.d0

Do i = 2, N-1
    dw(i) = X(i)-X(i-1)
    de(i) = X(i+1)-X(i)
    dx(i) = (dw(i)+de(i))/2.d0
End Do

dw(N) = X(N)-X(N-1)
de(N) = L_tot - X(N)
dx(N) = (dw(N)+de(N))/2.d0

!========================
! Noeuds internes
!========================
Do i = 2, N-1

    kw = 2.d0*K_vol(i)*K_vol(i-1)/(K_vol(i)+K_vol(i-1))
    ke = 2.d0*K_vol(i)*K_vol(i+1)/(K_vol(i)+K_vol(i+1))

    aw(i) = kw/dw(i)
    ae(i) = ke/de(i)
    ap(i) = aw(i) + ae(i)

    If (abs(alpha) > 1.d-10) Then
        b(i) = S_vol(i)/alpha * &
        (exp(-alpha*(X(i)-dw(i)/2.d0)) - exp(-alpha*(X(i)+de(i)/2.d0)))
    Else
        b(i) = S_vol(i)*dx(i)
    End If

End Do

!========================
! CL GAUCHE (CORRIGÉ)
!========================
ke = 2.d0*K_vol(1)*K_vol(2)/(K_vol(1)+K_vol(2))
ae(1) = ke/de(1)
aw(1) = 0.d0

ap(1) = ae(1) + hw
b(1)  = hw*T_ext_w + qw

If (abs(alpha) > 1.d-10) Then
    b(1) = b(1) + S_vol(1)/alpha * &
    (1.d0 - exp(-alpha*(X(1)+de(1)/2.d0)))
Else
    b(1) = b(1) + S_vol(1)*dx(1)
End If

!========================
! CL DROITE (CORRIGÉ)
!========================
kw = 2.d0*K_vol(N)*K_vol(N-1)/(K_vol(N)+K_vol(N-1))
aw(N) = kw/dw(N)
ae(N) = 0.d0

ap(N) = aw(N) + he
b(N)  = he*T_ext_e - qe

If (abs(alpha) > 1.d-10) Then
    b(N) = b(N) + S_vol(N)/alpha * &
    (exp(-alpha*(X(N)-dw(N)/2.d0)) - exp(-alpha*L_tot))
Else
    b(N) = b(N) + S_vol(N)*dx(N)
End If

!========================
! TDMA
!========================
denom = ap(1)
P(1) = ae(1)/denom
Q(1) = b(1)/denom

Do i = 2, N
    denom = ap(i) - aw(i)*P(i-1)
    P(i) = ae(i)/denom
    Q(i) = (b(i) + aw(i)*Q(i-1))/denom
End Do

T(N) = Q(N)
Do i = N-1,1,-1
    T(i) = P(i)*T(i+1) + Q(i)
End Do

!========================
! Output
!========================
Open(20,file="sol_multicouche.txt")
Do i = 1, N
    Write(20,*) X(i), T(i), T(i)-273.15
End Do
Close(20)

Print*, "Resultats dans sol_multicouche.txt"

End Program batterie_multicouche