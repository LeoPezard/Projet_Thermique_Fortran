import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def lire_fichier(fichier):
    if not os.path.exists(fichier):
        print(f"Erreur : fichier '{fichier}' introuvable.")
        sys.exit(1)
    with open(fichier, 'r') as f:
        lines = [l for l in f.readlines() if not l.strip().startswith('#') and l.strip()]
    N_cell = int(lines[-1].strip())
    data = np.array([list(map(float, l.split())) for l in lines[:-1]])
    x_m  = data[:, 0]
    T_C  = data[:, 2]
    return x_m * 1e6, T_C, N_cell

# ── Chargement des deux fichiers ──────────────────────────────────────────
fichier_lente  = sys.argv[1] if len(sys.argv) > 1 else "donnes_lente.txt"
fichier_rapide = sys.argv[2] if len(sys.argv) > 2 else "sol_multicouche.txt"

x_lente,  T_lente,  N_cell_l = lire_fichier(fichier_lente)
x_rapide, T_rapide, N_cell_r = lire_fichier(fichier_rapide)

# ── Tracé ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(13, 5))
fig.suptitle("Profil de température — Charge lente vs Charge rapide",
             fontsize=12, fontweight='bold')

# --- Bandes de couleur par couche (basées sur N_cell du fichier lente) ---
couches = [
    (0,   20,  '#B0BEC5', 'Al cathode'),
    (20,  120, '#A5D6A7', 'NMC actif'),
    (120, 145, '#FFF9C4', 'Séparateur'),
    (145, 235, '#FFCC80', 'Graphite actif'),
    (235, 250, '#EF9A9A', 'Cu anode'),
]
L_cell = 250  # µm

for rep in range(N_cell_l):
    offset = rep * L_cell
    for x_start, x_end, couleur, label in couches:
        lbl = label if rep == 0 else None
        ax.axvspan(offset + x_start, offset + x_end,
                   alpha=0.3, color=couleur, label=lbl)

# --- Seuils de température ---
ax.axhline(45, color='orange', linestyle='--', linewidth=1.2, label='Limite conception (45°C)')
ax.axhline(60, color='red',    linestyle='--', linewidth=1.2, label='Limite absolue (60°C)')

# --- Profils ---
ax.plot(x_lente,  T_lente,  color='steelblue', lw=2.0, marker='s',
        markersize=2, label='Charge lente')
ax.plot(x_rapide, T_rapide, color='#E53935',   lw=2.0, marker='s',
        markersize=2, label='Charge rapide')

ax.set_xlabel("Position x (µm)", fontsize=11)
ax.set_ylabel("Température (°C)", fontsize=11)
ax.set_title("Profil T(x) en °C", fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_xlim(x_lente[0], x_lente[-1])
ax.legend(loc='upper right', fontsize=9, framealpha=0.8)
plt.ticklabel_format(style='plain')
plt.tight_layout()

# ── Sauvegarde ────────────────────────────────────────────────────────────
nom_fig = "sol_comparaison_plot.png"
plt.savefig(nom_fig, dpi=150, bbox_inches='tight')
print(f"Figure sauvegardée : {nom_fig}")
plt.show()