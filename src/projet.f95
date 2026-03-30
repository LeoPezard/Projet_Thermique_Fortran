Program batterie_multicouche
Implicit None
 
!========================
! Déclarations
!========================
Integer :: N, i, j, Ncouches
Integer, Parameter :: N_rep = 100 ! A changer pour tester le nombre de cellules voulues
Real*8 :: L_tot, L_batterie, dx_bat
Real*8 :: R_tot_cell, K_eq_cell, S_eq
Real*8 :: e_contact, K_contact, S_contact, K_eq_asm, R_tot_asm
Real*8, Allocatable :: K_couche(:), S_couche(:), e_couche(:)
Real*8, Allocatable :: K_vol(:), S_vol(:)
 
Real*8 :: hw, T_ext_w, qw, he, T_ext_e, qe, alpha
 
Real*8, Allocatable :: T(:), ap(:), ae(:), aw(:), b(:)
Real*8, Allocatable :: X(:), dw(:), de(:), dx(:)
Real*8, Allocatable :: P(:), Q(:)
 
Real*8 :: denom, ke, kw
 
!========================
! Lecture données
!========================
Open(15, file="data_lente.txt")
!Open(15, file="data_rapide.txt")
 
Read(15,*)
Read(15,*) N
Read(15,*) Ncouches
 
Allocate(e_couche(Ncouches), K_couche(Ncouches), S_couche(Ncouches))
 
Do j = 1, Ncouches
    Read(15,*) e_couche(j), K_couche(j), S_couche(j)
End Do

Read(15,*) e_contact, K_contact, S_contact


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
! Épaisseurs
!========================
L_tot      = sum(e_couche)
L_batterie = L_tot * dble(N_rep)
 
!========================
! Conductivité équivalente (résistances en série)
!========================
R_tot_cell = 0.d0
Do j = 1, Ncouches
    R_tot_cell = R_tot_cell + e_couche(j) / K_couche(j)
End Do
K_eq_cell = L_tot / R_tot_cell
 

R_tot_asm = dble(N_rep) * L_tot / K_eq_cell + dble(N_rep - 1) * ( e_contact / K_contact)
K_eq_asm  = L_batterie / R_tot_asm


!========================
! Source équivalente (moyenne pondérée par épaisseur)
!========================
S_eq = 0.d0
Do j = 1, Ncouches
    S_eq = S_eq + S_couche(j) * e_couche(j)
End Do
S_eq = S_eq / L_tot
 
Print*, "K_eq_cell (W/mK) = ", K_eq_cell
Print*, "S_eq      (W/m3) = ", S_eq
Print*, "L_batterie  (m)  = ", L_batterie
 
!========================
! Allocation
!========================
Allocate(X(N), T(N), ap(N), aw(N), ae(N), b(N))
Allocate(dw(N), de(N), dx(N), K_vol(N), S_vol(N))
Allocate(P(N), Q(N))
 
!========================
! Maillage uniforme sur L_batterie
!========================
dx_bat = L_batterie / dble(N)
Do i = 1, N
    X(i) = (dble(i) - 0.5d0) * dx_bat
End Do
 
!========================
! Propriétés locales (homogène équivalent)
!========================
Do i = 1, N
    K_vol(i) = K_eq_asm
    S_vol(i) = S_eq
End Do
 
!========================
! Distances (maillage uniforme)
!========================
dw(1) = X(1)
de(1) = dx_bat
dx(1) = (dw(1) + de(1)) / 2.d0
 
Do i = 2, N-1
    dw(i) = dx_bat
    de(i) = dx_bat
    dx(i) = dx_bat
End Do
 
dw(N) = dx_bat
de(N) = L_batterie - X(N)
dx(N) = (dw(N) + de(N)) / 2.d0
 
!========================
! Noeuds internes
!========================
Do i = 2, N-1
 
    kw = 2.d0*K_vol(i)*K_vol(i-1) / (K_vol(i)+K_vol(i-1))
    ke = 2.d0*K_vol(i)*K_vol(i+1) / (K_vol(i)+K_vol(i+1))
 
    aw(i) = kw / dw(i)
    ae(i) = ke / de(i)
    ap(i) = aw(i) + ae(i)
 
    If (abs(alpha) > 1.d-10) Then
        b(i) = S_vol(i)/alpha * &
               (exp(-alpha*(X(i)-dw(i)/2.d0)) - exp(-alpha*(X(i)+de(i)/2.d0)))
    Else
        b(i) = S_vol(i) * dx(i)
    End If
 
End Do
 
!========================
! CL GAUCHE
!========================
ke = 2.d0*K_vol(1)*K_vol(2) / (K_vol(1)+K_vol(2))
ae(1) = ke / de(1)
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
! CL DROITE
!========================
kw = 2.d0*K_vol(N)*K_vol(N-1) / (K_vol(N)+K_vol(N-1))
aw(N) = kw / dw(N)
ae(N) = 0.d0
 
ap(N) = aw(N) + he
b(N)  = he*T_ext_e - qe
 
If (abs(alpha) > 1.d-10) Then
    b(N) = b(N) + S_vol(N)/alpha * &
           (exp(-alpha*(X(N)-dw(N)/2.d0)) - exp(-alpha*L_batterie))
Else
    b(N) = b(N) + S_vol(N)*dx(N)
End If
 
!========================
! TDMA
!========================
denom = ap(1)
P(1)  = ae(1) / denom
Q(1)  = b(1)  / denom
 
Do i = 2, N
    denom = ap(i) - aw(i)*P(i-1)
    P(i)  = ae(i) / denom
    Q(i)  = (b(i) + aw(i)*Q(i-1)) / denom
End Do
 
T(N) = Q(N)
Do i = N-1, 1, -1
    T(i) = P(i)*T(i+1) + Q(i)
End Do
 
!========================
! Output
!========================
Open(20, file="sol_multicouche.txt")
Write(20,*) "# X(m)   T(K)   T(C)"
Do i = 1, N
    Write(20,*) X(i), T(i), T(i)-273.15d0
End Do
Write(20,*) N_rep   ! en tête du fichier
Close(20)
 
Print*, "Resultats dans sol_multicouche.txt"
 
End Program batterie_multicouche