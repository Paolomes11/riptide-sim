# Risultati ereditati dal toy model

Sintesi di quanto stabilito dal toy model Python (`github.com/Paolomes11/RIPTIDE_toy_model`, v0.26,
2026-11-14), da usare come specifica e come termine di confronto per `riptide-sim`. Il dettaglio
(tabelle complete, script, previsioni scritte prima di ogni run) resta in `docs/resoconto.md` e
`docs/roadmap.md` di quel repository.

Etichette: **(a)** fatto consolidato · **(b)** stima toy con assunzioni esplicite · **(c)** verificato
con Monte Carlo (nel toy) · **(d)** ipotesi da testare.

**Avvertenza di fondo.** Nel toy il generatore degli eventi e il modello di ricostruzione sono lo
stesso modello (cinematica a due corpi, isotropia in CM, risoluzioni gaussiane): le etichette (c)
dicono che il motore statistico è *autoconsistente*, non che descrive un rivelatore reale. È il
motivo di questo progetto: con Geant4 il generatore diventa indipendente, e ogni risultato (c) va
rimisurato.

## 1. Modello

Un neutrone di energia E_n diffonde elasticamente su un protone; si misura il protone di rinculo
(energia E_p, direzione della traccia). θ_p = angolo fra neutrone incidente e protone in lab.

| Caso | Parametro condiviso | Nuisance per evento | Assunzioni |
|---|---|---|---|
| A | E_n | θ_p | 1 sorgente unica, monoenergetica, direzione nota |
| B | Ω_n (direzione del neutrone) | E_n (prior largo) | 1 + 2 campo lontano |
| C | Ω_n e (μ_E, σ_E) | E_n^(k) ~ N(μ_E, σ_E) troncata al dominio | 1 + 2 + 3 energie simili |

Assunzioni: **1** sorgente unica (necessaria, da testare) · **2** direzione ~costante (campo
lontano) · **3** energie simili. L'indipendenza condizionata degli eventi si testa (pull,
coverage), non si assume.

Parametri del toy (b): dominio E_n ∈ [10, 20] MeV; σ_Ep = 0.10 MeV, σ_θ = 0.08 rad (gaussiane);
prior sempre propri: E_n piatto su [10, 20], σ_E log-uniforme su [0.01, 50] MeV, θ_p con densità
sin 2θ_p (Caso A). Seme comune `20260907`. Massa del protone 938.272 MeV (CODATA); masse n = p
nella cinematica (differenza ~5·10⁻⁷ sull'energia massima di rinculo, trascurabile).

Combinazione degli eventi: log p(ψ | dati) = log π(ψ) + Σ_k log L_k(ψ), prior contato **una volta**.

## 2. Cinematica

| Risultato | Valore | |
|---|---|---|
| Energia del protone (relativistica, masse uguali m) | E_p = 2m E_n cos²θ_p / (2m + E_n sin²θ_p), θ_p ≤ 90° | (a) |
| Limite non relativistico | E_p = E_n cos²θ_p; a 20 MeV e 45° lo scarto è 0.05 MeV = 0.5 σ_Ep | (a) |
| Esempio | E_n = 20 MeV, θ_p = 45° ⇒ E_p = 9.947 MeV (non rel.: 10.000) | (a) |
| Isotropia in CM ⇒ spettro del protone | E_p ~ U(0, E_n) esattamente, anche relativistico | (a) |
| Isotropia in CM ⇒ angolo di rinculo | p(θ_p) = sin 2θ_p nel limite non relativistico | (a) |
| Densità delle tracce attorno a Ω_n | cos θ_p / π per sr (non rel.); correzione ±1% a 20 MeV | (a) |
| Anisotropia vera n–p a 10–20 MeV (ENDF/B-VIII.0, H-1) | presente, effetto trascurabile sulla ricostruzione (§6) | (c) |

**Prima verifica non circolare per la fase 1:** con `G4NeutronHP` lo spettro E_p del primo
scattering n–p deve risultare piatto su [0, E_n] e la relazione E_p(θ_p) deve seguire la formula
relativistica sopra; le deviazioni misurano l'anisotropia vera nei dati G4NDL.

## 3. Casi di riferimento (oracoli)

Da riportare come test Catch2 nel motore C++ (fase 5); i valori sono indipendenti dal dominio
10–20 MeV.

| Oracolo | Ingresso | Atteso | |
|---|---|---|---|
| R1, singolo evento | Ê_p = 1.40 MeV, θ̂_p = 0.50 rad, σ_Ep = 0.10, σ_θ = 0.08, prior **piatto** su E_n (griglia 0.5–6 MeV) | E_n = 1.82 ± 0.21 MeV (toy: σ pesata 0.219; tolleranza 0.01 su σ) | (c) |
| R2, due eventi | log-verosimiglianze gaussiane 2.70 ± 0.45 e 2.35 ± 0.20 MeV, prior piatto | 2.40 ± 0.18 MeV (toy: 2.408 ± 0.183) | (c) |
| R3, ricostruzione difettosa | scala +3%, offset +0.05 MeV, σ dichiarate 30% più strette della vera σ(E_n) = 0.08 E_n (rumore applicato **dopo** l'errore di scala) | pull medio ≈ 0.9, larghezza ≈ 1.4, coverage 68% → ~45% (toy: 0.90, 1.44, 0.433) | (c) |

Test di limite (da mantenere): Ω_n ∥ z ⇒ Caso A · σ_E → 0 ⇒ A · σ_E → ∞ ⇒ B · un evento con
prior piatto ⇒ L_k(Ω_n) circa costante sull'emisfero anteriore.

## 4. Calibrazione nel toy (10–20 MeV)

**Caso A** (c). Checklist 20 000 eventi con verità dal prior: pull −0.000 / 1.001, coverage
0.680 / 0.898 / 0.950. Risoluzione per evento su E_n ~1.3 MeV (dominata da σ_θ), larga rispetto al
dominio ⇒ forte attrazione verso il centro del prior (bias per bin +1.63 … −1.83 MeV). Contrazione
std·√N = 1.38 / 0.98 / 0.93 MeV a N = 10 / 100 / 1000: ≈ 1/√N da N ≈ 100, nessun plateau.
**Il prior su θ_p conta:** col prior piatto (incoerente con la cinematica) a N = 1000 si ha bias
−0.101 ± 0.009 MeV, pull −3.9 (plateau); col prior fisico sin 2θ_p bias −0.006 ± 0.008. Uno shift
sistematico del 2% su E_p produce un plateau a ~0.02·E_n.

**Caso B** (c). Calibrato a ogni N (10 … 1000): pull rms 0.98–1.02, coverage HPD nominale;
rms·√N ≈ 0.24 rad, pendenza di contrazione −0.51/−0.53. A N = 10 code un po' più pesanti di una
gaussiana 2D. Stima = MAP su calotta raffinata di una griglia di Fibonacci.

**Caso C** (c). Ricostruzione a due stadi (Ω_n come B, poi griglia 2D su μ_E, σ_E con prior
predittivo), mai griglia 4D. Calibrato a N ≥ 150; a N = 50 la coverage di log σ_E è bassa
(0.61 / 0.81 / 0.87), probabile effetto di prior (d). La correlazione Ω_n–(μ_E, σ_E) trascurata dai
due stadi costa ≤ 1.5% sulle incertezze in media (≤ 8% nel caso peggiore): fattorizzazione
giustificata. La calibrazione richiede che le energie vere siano generate dalla gaussiana
**troncata a [10, 20] MeV**: senza troncatura la frazione fuori dominio arriva a 0.107 e la
calibrazione di μ_E e log σ_E si rompe, sempre peggio al crescere di N.

## 5. Stress test delle assunzioni (toy)

| Stress | Esito | |
|---|---|---|
| Prior sui parametri condivisi (E_n in A, Ω_n in B, μ_E in C) | dimenticato con N | (c) |
| Prior sul nuisance per evento (E_n nel Caso B) | **non** dimenticato: Ω̂_n si sposta ~0.75 σ a ogni N con prior largo; è rumore a media nulla (relazione tipo Hausman): costa precisione, non calibrazione | (c) |
| Prior su σ_E (uniforme vs log-uniforme) | spostamento di log σ_E +0.76 / +0.13 / +0.07 σ a N = 50 / 300 / 1000 | (c) |
| Forma dello spettro (Caso C) | μ_E e Ω_n robusti; la coverage di log σ_E dipende dalla curtosi (sovra-copertura per spettri piatti, sotto-copertura per code pesanti) | (c) |
| Spettro che sborda dal dominio | danno che **cresce con N** (coverage 68% di μ_E 0.38–0.55 a N = 1000); meccanismo non indagato | (c), causa (d) |
| Sorgente unica (seconda sorgente o fondo isotropo, frazione f, separazione Δ) | Ω̂_n spostata ≈ f·Δ; con f = 0.05 a 30° coverage 68% di Ω_n 0.07 a N = 300; diagnostica KS cieca per Δ ≲ 30°; nel Caso C μ_E spinto al bordo superiore, collasso del posterior in 12/2800 esperimenti | (c) |
| Anisotropia vera in CM | trascurabile: bias aggiunto su μ_E 0.15 σ a N = 1000 (0.0128 MeV), da ricontrollare oltre N ≈ 1000 | (c) |

## 6. Limiti dichiarati dal toy (da affrontare qui)

- **Sezione d'urto n–p decrescente con E_n**: lo spettro rivelato non è quello incidente
  (efficienza). Geant4 la include da sé: va messa nella verosimiglianza (fase 9).
- **Soglia su E_p basse**: assente nel toy; nel rivelatore reale c'è (fase 9).
- **Eventi n–C, scattering multipli, protoni che escono dal volume**: assenti nel toy; Geant4 li
  produce, la selezione dei candidati a singolo scattering diventa un passo esplicito (fase 4).
- **Sorgenti multiple / fondo**: un modello a mistura che stimi f è il passo naturale (fuori dal
  perimetro del toy).
- **Profondità d'interazione z** (nuisance del Caso A): nessun osservabile definito nel toy.
- **Approssimazioni non relativistiche residue** nel prior su θ_p e nella CDF angolare delle tracce
  (correzioni ~1% a 20 MeV, b).

Prossimi passi del toy rimasti aperti e che qui diventano naturali: σ_np(E) da ENDF/B-VIII.0
(MF3/MT2), anisotropia maggiorata, coverage di log σ_E a N = 50 con prior ristretto, pull medio
comune di μ_E (+0.053 ± 0.024, d).

## 7. Prestazioni del toy (riferimento da battere)

| Operazione | Toy (numpy) | |
|---|---|---|
| log-verosimiglianza Caso A, 1 evento | 0.67 ms | (c) |
| log-verosimiglianza Caso A, 1000 eventi | 0.71 s (target ridefinito da 50 ms a 2 s: limite del numpy denso) | (c) |
| Esperimento Caso C a N = 300 | ~2 s | (c) |
| Checklist completa Caso C | ~8 min con 3 processi | (c) |

Lezioni d'implementazione: niente cicli sugli eventi nel path caldo; griglie contigue riusate;
profilare prima di ottimizzare; dopo ogni ottimizzazione rilanciare gli oracoli (la velocità non
deve cambiare i numeri: nel toy un passaggio a float32 ha cambiato i risultati di 1.8·10⁻⁶
relativo, accettabile solo perché misurato).

## 8. Checklist di validazione (ordine obbligatorio)

Nessun risultato è "finale" senza: bias → risoluzione → pull N(0, 1) → coverage (68/90/95%) →
contrazione ≈ 1/√N (un plateau indica un sistematico condiviso) → robustezza al prior. Con M
esperimenti indipendenti per N, riportando l'errore binomiale sulla coverage. Le previsioni si
scrivono prima del run.
