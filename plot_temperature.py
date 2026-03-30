import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# ── Chargement du fichier ──────────────────────────────────────────────────
fichier = sys.argv[1] if len(sys.argv) > 1 else "sol_multicouche.txt"

if not os.path.exists(fichier):
    print(f"Erreur : fichier '{fichier}' introuvable.")
    sys.exit(1)

# Lecture robuste : ignore les lignes commentaires (#) et les en-têtes
data = np.loadtxt(fichier, comments='#')

x_m  = data[:, 0]          # Position (m)
T_K  = data[:, 1]          # Température (K)
T_C  = data[:, 2]          # Température (°C)

x_um = x_m * 1e6           # Conversion en µm pour l'affichage

print(f"Fichier       : {fichier}")
print(f"Nb de points  : {len(x_m)}")
print(f"x : [{x_um[0]:.3f} ; {x_um[-1]:.3f}] µm")
print(f"T : [{T_C.min():.4f} ; {T_C.max():.4f}] °C   ΔT = {T_C.max()-T_C.min():.4f} K")

# ── Tracé ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(13, 5))
fig.suptitle(f"Profil de température ",
             fontsize=12, fontweight='bold')

# --- Bandes de couleur par couche ---
couches = [
    (0,   20,  '#B0BEC5', 'Al cathode'),
    (20,  120, '#A5D6A7', 'NMC actif'),
    (120, 145, '#FFF9C4', 'Séparateur'),
    (145, 235, '#FFCC80', 'Graphite actif'),
    (235, 250, '#EF9A9A', 'Cu anode'),
]
for x_start, x_end, couleur, label in couches:
    ax.axvspan(x_start, x_end, alpha=0.3, color=couleur, label=label)

ax.plot(x_um, T_C, color='#E53935', lw=2.0, marker='s',
        markersize=3)
ax.set_xlabel("Position x (µm)", fontsize=11)
ax.set_ylabel("Température (C°)", fontsize=11)
ax.set_title("Profil T(x) en C°", fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_xlim(x_um[0], x_um[-1])
ax.legend(loc='lower right', fontsize=9, framealpha=0.8)
plt.ticklabel_format(style='plain')
plt.tight_layout()

# Sauvegarde automatique
nom_fig = os.path.splitext(fichier)[0] + "_plot.png"
plt.savefig(nom_fig, dpi=150, bbox_inches='tight')
print(f"Figure sauvegardée : {nom_fig}")
plt.show()
