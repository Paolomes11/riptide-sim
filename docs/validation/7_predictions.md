# Previsioni della fase 7 — scritte prima dei run

Versione 1 (2026-11-25). Va committata prima di qualsiasi riga di codice del Caso B, quindi prima
della checklist 7a e dei run di simulazione della fase 7. I risultati (`7_results.md`)
citeranno l'hash di quel commit. Dopo i run le soglie non si toccano.

Versione 2 (2026-11-25): aggiunge la §4.1, scritta dopo la 7a e prima del codice e dei run della
7b. Precisa come si leggono le soglie della §4 e non cambia nessun numero.

Versione 3 (2026-11-26): aggiunge la §4.2, la diagnosi 7c del pull rms dell'obliqua a N = 300,
scritta dopo la 7b e prima dei run della 7c. Le soglie della §4 non cambiano.

La fase 7 porta il Caso B del toy in C++: Ω_n incognita e condivisa, E_n nuisance per evento con
prior piatto su [10, 20] MeV, stima MAP su una griglia di Fibonacci con calotta raffinata
(D41, D42). Come nelle fasi 5 → 6, prima si valida il port contro il toy con pseudo-esperienze
(7a, D43), poi lo si applica a eventi Geant4 con tre direzioni fisse (7b, D45).

Fonti dei numeri:
- la calibrazione del Caso B nel toy (`toy_model_findings.md` §4 e la tabella completa del
  resoconto del toy, riportata qui sotto);
- lo stress test del prior sul nuisance (`toy_model_findings.md` §5);
- per i run, i tempi e lo spazio dei run isotropi da 10⁷ della fase 3 (`3_results.md`);
- per la parte Geant4, i risultati del Caso A (`6_results.md`) e le ipotesi di D45.

Etichette: **(a)** fatto stabilito (sorgenti, dati, costruzione del codice); **(b)** stima con
ipotesi esplicite; **(c)** verificato con Monte Carlo; **(d)** ipotesi da mettere alla prova.

## Riferimento: il toy (c)

σ_Ep = 0.10 MeV, σ_θ = 0.08 rad, verità isotropa, E_n ~ U(10, 20), griglie di D41–D42.

| N | M | errore medio (°) | rms (°) | risoluzione dichiarata (°) | pull rms | HPD 68/90/95 |
|---|---|---|---|---|---|---|
| 10 | 400 | 4.30 | 4.94 | 4.76 | 1.02 | 0.68 / 0.90 / 0.95 |
| 30 | 300 | 2.31 | 2.60 | 2.57 | 1.00 | 0.67 / 0.90 / 0.94 |
| 100 | 200 | 1.21 | 1.36 | 1.37 | 0.99 | 0.66 / 0.92 / 0.95 |
| 300 | 100 | 0.76 | 0.84 | 0.79 | 1.07 | 0.63 / 0.84 / 0.92 |
| 1000 | 40 | 0.35 | 0.37 | 0.43 | 0.85 | 0.80 / 0.98 / 1.00 |
| 1000 | 240 | — | 0.42 | — | 0.98 | 0.70 / 0.90 / 0.95 |

rms·√N = 0.27 / 0.25 / 0.24 / 0.25 / 0.23 rad; pendenza log-log −0.53 su 10–1000, −0.51 su
30–1000. A N = 10 la mediana di pull·√2 è 1.17 (Rayleigh: 1.18) e la frazione sopra 3 è 0.018.
A N = 300 e a N = 1000 con M = 40 gli scarti dal nominale sono compatibili con la statistica di
M (con M = 240 tornano nominali).

## 1. Oracoli e test di limite (Catch2)

| Test | Atteso | |
|---|---|---|
| Ω_n ∥ z, stesso prior su θ_p | il posterior marginale di E_n del Caso B per evento coincide col Caso A a meno della differenza fra kernel vMF e gaussiana 1D, O(σ_θ²) relativo; tolleranza dichiarata nel test prima di eseguirlo | (a) |
| Kernel K(θ_obs \| θ) | normalizzato a 1 su [0, π] entro l'errore di quadratura dichiarato; per κ grande tende alla gaussiana in θ_obs − θ | (a) |
| I0e | valori di riferimento di `scipy.special.i0e` riprodotti entro 1·10⁻¹² relativo su tutto l'intervallo usato | (a) |
| Un evento, prior piatto su E_n | L_k(Ω) dipende da Ω solo tramite θ_obs e, lontano dai bordi del dominio, ha la forma del toy (L ∝ π(E*)/cos θ, circa costante sull'emisfero anteriore) | (a) |
| `sphereGrid` | n direzioni unitarie distinte; area uguale: i conteggi in calotte casuali seguono l'area entro gli errori di Poisson | (a) |
| `directionCapGrid` | tutti i punti entro il raggio; cos α uniforme (KS) | (a) |
| Generatore B | CDF dell'angolo traccia–Ω vero (prima dello smearing) uguale a quella del modello (KS a p > 0.01); Ep ~ U(0, E_n) | (c) |
| Posterior | finito su tutta la sfera per ogni esperienza; altrimenti eccezione | (a) |
| Ottimizzazioni | numeri invariati, o scarti entro il limite dimostrato, rilanciando gli oracoli dopo ogni `perf` | (a) |

## 2. Port contro toy (7a)

Le stesse coppie (N, M) del toy, più N = 1000 con M = 240 su id disgiunti. Seme `20260907`. I
numeri non sono confrontabili bit per bit (RNG e linguaggio diversi), ma solo entro gli errori
statistici.

| Grandezza | Soglia, a ogni N | |
|---|---|---|
| Pull rms | 1 ± 0.07 per M ≥ 200; per M = 100 e M = 40, entro 2σ statistici dell'rms (σ ≈ 1/√(2M)) | (c) nel toy |
| Copertura HPD 68/90/95 | entro 2σ binomiali dal nominale; con M = 40 conta l'esito con M = 240 | (c) nel toy |
| rms·√N | fra 0.22 e 0.27 rad | (c) nel toy |
| Pendenza log-log dell'rms su 30–1000 | −0.5 ± 0.05 | (c) nel toy |
| Mediana di pull·√2 a N = 10 | 1.17 ± 0.05 | (c) nel toy |
| Frazione di pull·√2 > 3 a N = 10 | ≤ 0.04 (Rayleigh: 0.011; toy: 0.018) | (c) nel toy |
| `truth_in_cap` | 1 in tutte le esperienze | (b) |
| Errore medio, rms, risoluzione dichiarata | compatibili con la tabella del toy entro 2σ statistici | (c) nel toy |

## 3. Robustezza al prior (7a)

A N = 100, 300, 1000, sulle stesse esperienze della §2.

| Variante | Atteso | |
|---|---|---|
| Prior vMF su Ω, κ = 1, centrato a 90° dalla verità | spostamento del MAP < 0.3 σ (σ = risoluzione dichiarata) a N ≥ 100: il prior sul parametro condiviso si dimentica. Stima: σ·gradiente del log-prior ≈ 0.02–0.09 σ | (b) |
| Prior su E_n piatto su [10, 30] invece di [10, 20] | spostamento a media nulla (media entro 2σ statistici da 0) con dispersione ≈ 0.75 σ a ogni N: il prior sul nuisance non si dimentica, è rumore a media nulla. Pull e copertura restano calibrati | (c) nel toy |

## 4. Geant4 (7b)

Tre run a campo lontano piatti in 10–19.5 MeV, 10⁷ eventi ciascuno (D45). Catena
`single_np,has_proton,contained`, conversione HP, smearing D31. Modello ideale (r ≡ 1) salvo
dove indicato. Per ogni direzione si confronta con il toy gemello (stessa verità, stessi N e M).

| Grandezza | +X | +Z | (3, 2, 1)/√14 | |
|---|---|---|---|---|
| Bias vettoriale di Ω̂ a N = 1000 | entro 2σ da 0 | entro 2σ da 0 | < 0.5 σ_c (σ_c = risoluzione per esperienza) | (a) per simmetria su +X e +Z; (d) sull'obliqua |
| Risoluzione G4 / toy gemello | ≤ 1.15 | ≤ 1.15 | ≤ 1.15 | (d) |
| Pull rms | 1 ± 0.1 | 1 ± 0.1 | 1 ± 0.1 | (d) |
| Copertura HPD 68/90/95 | entro 3σ binomiali | entro 3σ binomiali | entro 3σ binomiali | (d) |

Le soglie Geant4 si applicano al modello ideale e a quello con `--acceptance` (r di D40,
approssimata per +Z e per l'obliqua). Se falliscono: diagnosi con una decisione a parte, nessuna
correzione delle soglie.

### 4.1 Precisazioni prima dei run (versione 2)

- **Natura del confronto.** Il toy gemello è un controllo di coerenza, non un bersaglio. Una
  differenza Geant4 − toy misura il costo del realismo del rivelatore e va spiegata. Non si
  regolano né le soglie né il modello per avvicinare Geant4 al toy. Un fallimento apre una
  diagnosi con una decisione a parte.
- **Righe.** Le stesse della 7a, con M per riga della fase 5 (D37): (N, M) = (10, 400),
  (30, 300), (100, 200), (300, 100), (1000, 40), più la ripetizione (1000, 240) su indici
  disgiunti. Il toy gemello usa gli stessi N e M.
- **N = 1000.** Bias, rapporto di risoluzione, pull e copertura si valutano sulla ripetizione
  M = 240. La riga M = 40 si riporta senza esito, per il motivo visto nella 7a (`7_results.md`
  §4.2).
- **Bias vettoriale.** Per ogni esperienza, le coordinate tangenti (e_θ, e_φ) di Ω̂ sulla base
  locale alla verità. Il bias è la loro media, con errore = deviazione standard / √M per
  componente.
  - "Entro 2σ da 0" (+X, +Z) vale componente per componente.
  - "< 0.5 σ_c" (obliqua) vale per la norma del vettore medio, con σ_c = media della risoluzione
    dichiarata nella riga.
- **Risoluzione G4 / toy.** Rapporto degli errori angolari rms, a ogni N. Il rapporto delle
  risoluzioni dichiarate medie si riporta accanto, senza soglia.
- **Pull rms e copertura HPD.** Si applicano a ogni N, con M della riga.
- **Valori del toy gemello.** Quelli non arrotondati del suo record `case_b.json`.
- **Run +X.** Si riusa `6_ff_flat_10_19p5MeV_x_1e7` della fase 6 (config
  `6_ff_flat_10_19p5MeV_x.json`, c1 PASS): coincide con `7_ff_flat_10_19p5MeV_x` di D45, che
  quindi non si rilancia. Nuovi run solo per +Z e per l'obliqua.
- **Accettanza.** `--acceptance 6acc_s20260908`: la r(E_n, θ_p) di D40, misurata lungo +X sul
  run indipendente `6_ff_flat_10_19p5MeV_x_1e7_s20260908`. Per +Z e l'obliqua è
  un'approssimazione dichiarata (D45).

### 4.2 Diagnosi 7c: l'obliqua con un altro split (versione 3)

Nella 7b l'unico fallimento è il pull rms dell'obliqua a N = 300, M = 100: 0.859 senza
accettanza e 0.845 con l'accettanza, sulle stesse esperienze (`7_results.md` §6.4). L'ipotesi
(d) è una fluttuazione dello split. Qui la si mette alla prova.

- **Chiamate.** Le stesse quattro della 7b per l'obliqua, cambiando solo `--split-seed 20260908`.
  Due chiamate Geant4 sul run `7_ff_flat_10_19p5MeV_obl_1e7`, senza e con
  `--acceptance 6acc_s20260908`, e i due toy gemelli con lo stesso numero di eventi (222456).
  Restano invariati `--seed 20260907` (smearing e toy), la catena dei tagli, le griglie e le
  risoluzioni. I record sono `stats/7c_obl{,_acc}{,_toy}`.
- **Cosa cambia.** Cambia solo la permutazione del pool, quindi l'assegnazione degli eventi alle
  esperienze; gli eventi e il loro smearing sono gli stessi. Le 30 000 righe di N = 300
  condividono con quelle della 7b in media 30 000² / 222 456 ≈ 4.0·10³ eventi, cioè circa il
  13.5% **(b)**. Le esperienze sono quindi quasi indipendenti, ma non lo è il pool: un difetto
  legato a questo campione di eventi non si vede con questo test.
- **Soglie.** Quelle della §4 lette come nella §4.1, senza modifiche: pull rms 1 ± 0.1,
  risoluzione G4 / toy ≤ 1.15, copertura HPD entro 3σ binomiali, bias < 0.5 σ_c a N = 1000 con
  la ripetizione (qui M = 182). Si valutano tutte le 25 voci per ciascuna delle due chiamate
  Geant4.
- **Attese (b).**
  - Con M = 100 l'errore del pull rms è ≈ 1/√(2M) = 0.071, quindi una riga corretta esce da
    1 ± 0.1 con probabilità ≈ 16%.
  - Sulle 25 voci di una chiamata si attende ≈ 0.4 fallimenti di pull da falso allarme (righe
    con M da 100 a 400). Le due chiamate condividono le esperienze e contano come una.
- **Criterio di esito della diagnosi.**
  - **Fluttuazione (d → c):** il pull rms a N = 300 cade in 1 ± 0.1 in tutte e due le
    chiamate e rms·√N torna nella banda 0.22–0.27 rad della 7a.
  - **Non conclusiva:** il pull rms a N = 300 cade in 1 ± 0.1 in una sola chiamata, oppure il
    fallimento si sposta su un'altra riga con un valore compatibile con il falso allarme atteso.
    La si riporta così, senza ripetere lo split cercando l'esito voluto.
  - **Difetto possibile:** di nuovo sotto 0.9 a N = 300. Allora si apre una diagnosi del
    modello con una decisione a parte; la prima candidata è la parte della selezione che
    dipende dalla direzione, rimasta fuori dal modello (D45).
- Il fallimento della 7b resta registrato com'è: la 7c lo spiega, non lo sostituisce.

## 5. Tempi e risorse

| Voce | Atteso | |
|---|---|---|
| Un'esperienza a N = 1000, thread singolo, `release` | ≤ 2 s (D44) | (d) |
| Checklist 7a completa (~1.0·10⁵ eventi più le varianti di prior) | ≤ 10 min con 4 thread, RSS < 1 GB | (b) |
| Risultati con 1 e 4 thread | identici bit per bit | (a) per costruzione |
| Run sim da 10⁷ eventi | ~150 s e ~1.2 GB su SSD ciascuno, RSS < 0.2 GB (run isotropi della fase 3: 145 s, 180 MB) | (b) |
