import matplotlib.pyplot as plt
import numpy as np

# Données d'exemple (remplacez par vos vraies données)
versions = ['version0', 'version1', 'version2', 'version3', 'version4', 'version5']
temps_execution = [2.69058, 0.2932, 0.144592 , 0.0845588, 0.0811896, 0.072584]  

# Création du graphique
plt.figure(figsize=(10, 6))
plt.bar(versions, temps_execution, color='steelblue', edgecolor='black')
plt.xlabel('Version du Solver', fontsize=12)
plt.ylabel('Temps d\'exécution (secondes)', fontsize=12)
plt.title('Comparaison des temps d\'exécution des différentes versions', fontsize=14)
plt.grid(axis='y', alpha=0.3)

# Afficher les valeurs sur les barres
for i, v in enumerate(temps_execution):
    plt.text(i, v + 0.2, f'{v:.2f}s', ha='center', fontsize=10)

plt.tight_layout()
plt.savefig('timing_comparison.png', dpi=300)

# Graphique du speedup
plt.figure(figsize=(10, 6))
speedup = [temps_execution[0] / t for t in temps_execution]
plt.plot(versions, speedup, marker='o', linewidth=2, markersize=8, color='green')
plt.xlabel('Version du Solver', fontsize=12)
plt.ylabel('Speedup (par rapport à version0)', fontsize=12)
plt.title('Accélération relative des optimisations', fontsize=14)
plt.grid(True, alpha=0.3)
plt.axhline(y=1, color='r', linestyle='--', label='Baseline (version0)')
plt.legend()

for i, v in enumerate(speedup):
    plt.text(i, v + 0.1, f'{v:.2f}x', ha='center', fontsize=10)

plt.tight_layout()
plt.savefig('speedup_comparison.png', dpi=300)

print("Graphiques sauvegardés: timing_comparison.png et speedup_comparison.png")