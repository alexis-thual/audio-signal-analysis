# Fenêtrage et Calcul de Transformée de Fourier rapide

1. On a défini $N = 0,7 \times F_s$, $F_{min} = 50$Hz et $F_{max} = 900$Hz, donc 0,7 seconde de signal suffit amplement pour détecter ces fréquences.

2. L'offset sert à ne pas prendre en compte les toutes premières fréquences du signal afin de ne pas ce qui permet d'éviter éventuellement une partie du bruit que contient le signal.

3. La précision est exactement $F_s$.

# Estimation de la fréquence fondamentale par la méthode du produit spectral

2. 2. $H \frac{R-1}{N_{fft}}F_s \leq F_s / 2 \Leftrightarrow R \leq \frac{N_{fft}}{2H} + 1$

# Soustraction du son correspondant à la fréquence fondamentale détectée

3. 1. Pour le piano, on s'attend à ce que la fréquence fondamentale soit bien démarquée (percusion), il est raisonnable de penser qu'une imprécision de l'ordre de 1% suffit. Pour le hautbois, 5% pourrait probablement suffire.

Frames par seconde $F_s$

TTF : $N_{fft}$ fréquences : $\forall 0 \leq k \leq N_{fft}, \quad f=\frac{k}{N_{fft}} \mod 2\pi$

Produit : $R$ valeurs avec $H(R-1) \leq N_{fft}$
