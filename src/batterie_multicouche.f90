 Program batterie_multicouche
!========================================================
! Mur 1D Multicouche Stationnaire - Volumes Finis
! Modele thermique d'une cellule de batterie Li-Ion
! Methode: Volumes Finis avec conductivites harmoniques
!========================================================
      Implicit None

!     Declarations
      Integer :: N, i, j, Ncouches
      Integer, Allocatable :: Nvol(:)        ! nb de volumes par couche
      Real*8  :: L_tot
      Real*8, Allocatable :: K_couche(:), S_couche(:), e_couche(:)
      Real*8, Allocatable :: K_vol(:)        ! conductivite de chaque volume
      Real*8, Allocatable :: S_vol(:)        ! source de chaque volume
      Real*8  :: hw, T_ext_w, he, T_ext_e
      Real*8, Allocatable :: T(:), ap(:), ae(:), aw(:), b(:)
      Real*8, Allocatable :: X(:), dw(:), de(:), dx(:)
      Real*8, Allocatable :: P(:), Q(:)
      Real*8  :: denom, ke, kw
      ! Variables pour fichier de sortie
      Character(len=50) :: mode_charge
      Real*8  :: T_max, T_min, T_moy
      Integer :: i_max

!     ---- Lecture des donnees ----
      Open(unit=15, file="data_multicouche.dat", status='old')
      Read(15,*) mode_charge
      Read(15,*) N          ! Nombre total de volumes finis
      Read(15,*) Ncouches   ! Nombre de couches

      Allocate(e_couche(Ncouches), K_couche(Ncouches), S_couche(Ncouches), Nvol(Ncouches))

      Do j = 1, Ncouches
          Read(15,*) e_couche(j), K_couche(j), S_couche(j)
      End Do

      Read(15,*) hw
      Read(15,*) T_ext_w
      Read(15,*) he
      Read(15,*) T_ext_e
      Close(15)

!     ---- Repartition des volumes par couche (proportionnelle a l epaisseur) ----
      L_tot = 0.0d0
      Do j = 1, Ncouches
          L_tot = L_tot + e_couche(j)
      End Do

      Do j = 1, Ncouches
          Nvol(j) = Max(1, Nint(dble(N) * e_couche(j) / L_tot))
      End Do
      ! Ajustement pour que la somme = N
      N = 0
      Do j = 1, Ncouches
          N = N + Nvol(j)
      End Do

      Print*, "Mode de charge     : ", Trim(mode_charge)
      Print*, "Nombre total de VF : ", N
      Print*, "Epaisseur totale   : ", L_tot*1.0d6, " µm"

!     ---- Allocation ----
      Allocate(X(N), T(N), ap(N), aw(N), ae(N), b(N))
      Allocate(dw(N), de(N), dx(N), K_vol(N), S_vol(N))
      Allocate(P(N), Q(N))

!     ---- Affectation des proprietes par volume ----
      i = 0
      Do j = 1, Ncouches
          Do While (i < Sum(Nvol(1:j)))
              i = i + 1
              K_vol(i) = K_couche(j)
              S_vol(i)  = S_couche(j)
          End Do
      End Do

!     ---- Maillage uniforme par couche ----
      i = 0
      X(1) = 0.0d0
      Do j = 1, Ncouches
          Dim_j: Block
              Real*8 :: dx_j
              Integer :: k, i_start
              dx_j   = e_couche(j) / dble(Nvol(j))
              i_start = i
              Do k = 1, Nvol(j)
                  i = i + 1
                  If (i == 1) Then
                      X(1) = dx_j / 2.0d0
                  Else
                      X(i) = X(i-1) + dx_j
                  End If
              End Do
          End Block Dim_j
      End Do

!     ---- Distances geometriques ----
      dw(1) = X(1)
      de(1) = X(2) - X(1)
      dx(1) = (dw(1) + de(1)) / 2.0d0
      Do i = 2, N-1
          dw(i) = X(i) - X(i-1)
          de(i) = X(i+1) - X(i)
          dx(i) = (dw(i) + de(i)) / 2.0d0
      End Do
      dw(N) = X(N) - X(N-1)
      de(N) = L_tot - X(N)
      dx(N) = (dw(N) + de(N)) / 2.0d0

!     ---- Coefficients noeuds internes ----
      Do i = 2, N-1
          ! Conductance harmonique entre volumes adjacents
          kw = 2.0d0 * K_vol(i) * K_vol(i-1) / (K_vol(i) + K_vol(i-1))
          ke = 2.0d0 * K_vol(i) * K_vol(i+1) / (K_vol(i) + K_vol(i+1))
          aw(i) = kw / dw(i)
          ae(i) = ke / de(i)
          ap(i) = aw(i) + ae(i)
          b(i)  = S_vol(i) * dx(i)
      End Do

!     ---- Condition limite Gauche (convection) ----
      ke = 2.0d0 * K_vol(1) * K_vol(2) / (K_vol(1) + K_vol(2))
      ae(1) = ke / de(1)
      aw(1) = 0.0d0
      If (hw > 1.0d-10) Then
          ap(1) = ae(1) + hw
          b(1)  = S_vol(1) * dx(1) + hw * T_ext_w
      Else
          ! Adiabatique
          ap(1) = ae(1)
          b(1)  = S_vol(1) * dx(1)
      End If

!     ---- Condition limite Droite (convection) ----
      kw = 2.0d0 * K_vol(N) * K_vol(N-1) / (K_vol(N) + K_vol(N-1))
      aw(N) = kw / dw(N)
      ae(N) = 0.0d0
      If (he > 1.0d-10) Then
          ap(N) = aw(N) + he
          b(N)  = S_vol(N) * dx(N) + he * T_ext_e
      Else
          ! Adiabatique
          ap(N) = aw(N)
          b(N)  = S_vol(N) * dx(N)
      End If

!     ---- Resolution TDMA ----
      denom = ap(1)
      P(1)  = ae(1) / denom
      Q(1)  = b(1)  / denom
      Do i = 2, N
          denom = ap(i) - aw(i) * P(i-1)
          If (abs(denom) < 1.0d-12) Then
              Print*, "Erreur pivot nul au noeud", i; Stop
          End If
          P(i) = ae(i) / denom
          Q(i) = (b(i) + aw(i) * Q(i-1)) / denom
      End Do
      T(N) = Q(N)
      Do i = N-1, 1, -1
          T(i) = P(i) * T(i+1) + Q(i)
      End Do

!     ---- Ecriture resultats ----
      Open(unit=20, file='sol_multicouche.dat', status='replace')
      Write(20,'(A)') '# x(m)   T(K)   T(C)'
      Do i = 1, N
          Write(20,'(3E15.6)') X(i), T(i), T(i)-273.15d0
      End Do
      Close(20)

!     ---- Statistiques ----
      T_max = T(1); T_min = T(1); T_moy = 0.0d0; i_max = 1
      Do i = 1, N
          If (T(i) > T_max) Then; T_max = T(i); i_max = i; End If
          If (T(i) < T_min) T_min = T(i)
          T_moy = T_moy + T(i)
      End Do
      T_moy = T_moy / dble(N)

      Print*, "--- Resultats ---"
      Print'(A,F8.3,A)', " T_max = ", T_max-273.15d0, " C"
      Print'(A,F8.3,A)', " T_min = ", T_min-273.15d0, " C"
      Print'(A,F8.3,A)', " T_moy = ", T_moy-273.15d0, " C"
      Print'(A,F8.3,A)', " DeltaT= ", T_max-T_min,    " K"
      Print*, "Resultats dans sol_multicouche.dat"

      Deallocate(X,T,ap,aw,ae,b,dw,de,dx,K_vol,S_vol,P,Q)
      Deallocate(e_couche,K_couche,S_couche,Nvol)
      End Program batterie_multicouche
