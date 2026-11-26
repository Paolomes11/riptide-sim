# Previsioni della fase 8 — scritte prima dei run

Versione 1 (2026-11-26). Va committata prima di qualsiasi riga di codice del Caso C, quindi
prima della checklist 8a e delle chiamate Geant4 della 8b. I risultati (`8_results.md`)
citeranno l'hash di quel commit. Dopo i run le soglie non si toccano.

Versione 2 (2026-11-26): aggiunge la §2.1, la diagnosi 8c di μ_E a N = 150, scritta dopo la 8a e
prima del run della 8c. Le soglie delle §2–3 non cambiano.

La fase 8 porta il Caso C del toy in C++: Ω_n condivisa e incognita, E_n^(k) ~ N(μ_E, σ_E)
troncata al dominio, stima a due stadi (D46). Come nella fase 7, prima il port contro il toy con
pseudo-esperienze (8a, D47), poi i campioni Geant4 ricampionati dai run piatti della 7b
(8b, D48).

Fonti dei numeri:
- la checklist del Caso C nel toy, ultima versione dei risultati (2026-10-13), riletta per questo
  documento (tabella qui sotto);
- `toy_model_findings.md` §4–5 per la fattorizzazione, la troncatura e il prior su σ_E;
- `7_results.md` per i pool selezionati e il comportamento del Caso B su Geant4.

Etichette: **(a)** fatto stabilito (sorgenti, dati, costruzione del codice); **(b)** stima con
ipotesi esplicite; **(c)** verificato con Monte Carlo; **(d)** ipotesi da mettere alla prova.

## Riferimento: il toy (c)

σ_Ep = 0.10 MeV, σ_θ = 0.08 rad, dominio [10, 20] MeV, μ_E ~ U(13, 17) MeV,
σ_E ~ log-U(0.8, 2.4) MeV, Ω_n isotropa, E_n dalla gaussiana troncata. Prior baseline
uniforme in (μ_E, log σ_E). Stime: μ_E = media, log σ_E = mediana, Ω_n = MAP. Gli errori sono
"± errore della media".

| Voce | N = 50 (M = 200) | N = 150 (M = 200) | N = 300 (M = 100) | N = 1000 (M = 40) |
|---|---|---|---|---|
| bias μ_E (MeV) | +0.079 ± 0.039 | +0.029 ± 0.025 | +0.015 ± 0.021 | −0.014 ± 0.014 |
| bias log σ_E | −0.056 ± 0.057 | +0.004 ± 0.019 | +0.012 ± 0.011 | −0.002 ± 0.010 |
| rms μ_E (MeV) (std a posteriori media) | 0.553 (0.781) | 0.350 (0.328) | 0.207 (0.179) | 0.092 (0.086) |
| rms log σ_E (std media) | 0.805 (0.657) | 0.264 (0.220) | 0.113 (0.125) | 0.066 (0.066) |
| rms Ω_n (°) (risoluzione dichiarata) | 1.85 (1.75) | 0.88 (0.95) | 0.69 (0.66) | 0.34 (0.36) |
| pull μ_E, media / larghezza | +0.09 / 0.88 | +0.05 / 1.00 | +0.18 / 0.92 | −0.14 / 1.04 |
| pull log σ_E, media / larghezza | +0.07 / 0.97 | +0.11 / 0.95 | +0.11 / 0.89 | +0.00 / 1.00 |
| pull rms Ω_n | 1.07 | 0.93 | 1.04 | 0.98 |
| coverage μ_E 68/90/95 | 0.72 / 0.91 / 0.97 | 0.68 / 0.90 / 0.95 | 0.68 / 0.90 / 0.95 | 0.55 / 0.90 / 0.93 |
| coverage log σ_E 68/90/95 | 0.61 / 0.81 / 0.87 | 0.67 / 0.87 / 0.94 | 0.69 / 0.95 / 0.98 | 0.65 / 0.90 / 0.95 |
| coverage HPD Ω_n 68/90/95 | 0.64 / 0.86 / 0.92 | 0.73 / 0.95 / 0.98 | 0.63 / 0.93 / 0.97 | 0.65 / 0.93 / 0.98 |

- rms·√N: μ_E 3.9 / 4.3 / 3.6 / 2.9 MeV, log σ_E 5.7 / 3.2 / 2.0 / 2.1, Ω_n 0.23 / 0.19 / 0.21 /
  0.19 rad. Pendenze log-log su 150–1000: μ_E −0.71, log σ_E −0.73, Ω_n −0.50. μ_E e log σ_E
  contraggono più in fretta di 1/√N fra 150 e 1000: a N piccolo pesa il prior, e a N = 1000 il
  toy ha solo M = 40 (errore relativo dell'rms ≈ 11%).
- Robustezza a N = 150, prior uniforme in σ_E contro la baseline, in unità della deviazione
  standard a posteriori: spostamento di μ_E medio +0.01, massimo 0.32; di log σ_E medio
  +0.21, massimo 1.08. Coverage di log σ_E con il prior uniforme: 0.62 / 0.86 / 0.92.
- In `toy_model_findings.md` §5 lo spostamento di log σ_E per il prior è +0.76 / +0.13 /
  +0.07 σ a N = 50 / 300 / 1000: il prior su σ_E si dimentica lentamente.
- Un riassunto testuale precedente (2026-09-26) del toy riporta altri numeri: rms di μ_E più
  piccoli di un fattore ~7, cioè verità diverse da quelle della checklist attuale. Vale la
  tabella qui sopra, riletta dai risultati del 2026-10-13, che coincide con la coverage di
  log σ_E a N = 50 citata nei findings.

## 1. Oracoli e test di limite (Catch2)

| Test | Atteso | |
|---|---|---|
| Prior troncato | integrale trapezoidale 1 entro 1·10⁻¹², anche con μ_E fuori dal dominio e σ_E = 0.01 MeV | (a) |
| σ_E → 0 | marginale per evento = verosimiglianza del Caso B a E_n = μ_E, entro l'errore di quadratura dichiarato nel test prima di eseguirlo | (a) |
| σ_E → ∞ | base sommata = Caso B con prior piatto sul dominio, entro 1·10⁻⁹ relativo | (a) |
| Base per evento | uguale a `energyLogPosterior` del Caso B con log prior ≡ 0 (identità di costruzione) | (a) |
| R1 | log post −∞ su tutta la griglia → eccezione | (a) |
| Finestra dello stadio 2 | contiene il massimo globale e tutti i punti con log post > max − 10 | (a) |
| Generatore C | E_n dalla gaussiana troncata (KS a p > 0.01); senza troncatura la frazione fuori dominio segue Φ | (c) |
| Toy | un'esperienza fissata riproduce media e deviazione standard di μ_E e log σ_E del toy entro 1·10⁻³ relativo, sulla stessa base per evento | (b) |
| Ottimizzazioni | numeri invariati, o scarti entro il limite dimostrato, rilanciando gli oracoli dopo ogni `perf` | (a) |

## 2. Port contro toy (8a)

Righe (D47): (50, 200), (150, 200), (300, 100), (1000, 40), più la ripetizione (1000, 240) su
semi disgiunti. A N = 1000 l'esito si legge sulla ripetizione; la riga M = 40 si riporta per il
confronto. Dominio [10, 20] MeV, come il toy.

| Voce, a ogni N | Soglia | |
|---|---|---|
| Bias di μ_E e di log σ_E | \|bias\| < 0.25 × deviazione standard a posteriori media | (c) nel toy |
| rms di μ_E, log σ_E, Ω_n: C++ / toy | entro 1 ± 2√(1/M) (± 0.14 a M = 200, ± 0.20 a M = 100; a N = 1000, ± 0.24 per M = 240 contro il toy a M = 40) | (c) nel toy |
| Larghezza del pull di μ_E e log σ_E, pull rms di Ω_n | 1 ± 0.07 per M ≥ 200; per M = 100 entro 2σ statistici (σ ≈ 1/√(2M)) | (c) nel toy |
| Coverage μ_E, log σ_E, Ω_n a 0.68 / 0.90 / 0.95 | entro 2σ binomiali dal nominale | (c) nel toy |
| Pendenza log-log dell'rms di Ω_n su 150–1000 | −0.5 ± 0.07 | (c) nel toy |
| Pendenze di μ_E e log σ_E su 150–1000 | entro ± 0.15 dal toy (−0.71, −0.73) | (c) nel toy |
| `truth_in_cap` | 1 in tutte le esperienze | (b) |

**Eccezioni dichiarate.**
- Coverage di log σ_E a N = 50: bassa anche nel toy (0.61 / 0.81 / 0.87), probabile effetto del
  prior (d). Si riporta senza esito; atteso entro 2σ binomiali dal **toy**.
- Pull medio di μ_E: nel toy c'è un pull medio comune +0.053 ± 0.024 (d). Non è una voce con
  soglia; si riporta la media su tutte le righe.

**Falsi allarmi attesi (b).** Le voci con esito sono 17 per N (2 bias, 3 rms, 3 pull, 9 coverage),
meno le 3 dell'eccezione: 65 in tutto, più 3 pendenze. A 2σ la probabilità di un falso allarme è
circa 0.05 per voce, ma le coverage dello stesso parametro sono correlate. Attesi 2–4 fallimenti
senza difetto. Con più di 7 fallimenti, o con fallimenti concentrati su un parametro o un N,
si apre una diagnosi (8c).

### 2.1 Diagnosi 8c: μ_E a N = 150 con altri semi (versione 2)

Nella 8a (`stats/8a_w`) tre fallimenti su quattro cadono sulla riga N = 150, M = 200, tutti su
μ_E: larghezza del pull 1.0702, coverage 0.590 a 0.68 e 0.855 a 0.90 (`8_results.md` §2.1).
La std a posteriori media di μ_E è 0.306 MeV, contro 0.328 nel toy, con lo stesso rms (0.344
contro 0.350). L'ipotesi (d) è una fluttuazione delle 200 esperienze. Qui la si mette alla
prova.

- **Run.** Solo la riga N = 150, M = 200 della checklist, con il prior baseline, cambiando solo
  `--seed 20260908`. Restano invariati la sezione della riga, le verità (μ_E ~ U(13, 17) MeV,
  σ_E ~ log-U(0.8, 2.4) MeV, Ω_n isotropa), le risoluzioni, le griglie e la build. Il record è
  `stats/8c_n150`.
- **Cosa cambia.** Il seme entra in `perEventSeed` attraverso splitMix64, quindi verità, eventi e
  rumore sono flussi nuovi: le 200 esperienze sono indipendenti da quelle della 8a **(b)**.
- **Soglie.** Quelle della §2 per N = 150 e M = 200, senza modifiche: 17 voci (2 bias, 3 rms
  contro il toy, 3 pull, 9 coverage).
- **Attese (b).**
  - Con M = 200 l'errore della larghezza del pull è ≈ 1/√(2M) = 0.05, quindi ± 0.07 vale circa
    1.4σ. Una riga corretta esce dalla soglia con probabilità ≈ 16% per ciascuna delle tre voci
    di pull, ≈ 5% per ogni coverage e ogni rms.
  - Sulle 17 voci si attende circa un falso allarme.
- **Criterio di esito della diagnosi.**
  - **Fluttuazione (d → c):** passano le tre voci di μ_E fallite nella 8a (larghezza del pull,
    coverage a 0.68 e a 0.90), e la riga ha al più 2 fallimenti in tutto.
  - **Difetto possibile:** la larghezza del pull di μ_E supera di nuovo 1.07 **e** la coverage
    di μ_E a 0.68 cade di nuovo sotto la banda dei 2σ binomiali. Allora si apre una diagnosi
    del modello con una decisione a parte. La prima candidata è la larghezza del posterior di
    μ_E, cioè la std media più stretta di quella del toy: finestra e griglia fine dello stadio 2,
    poi la predittiva.
  - **Non conclusiva:** ogni altro esito, per esempio una sola delle due condizioni del difetto,
    oppure più di 2 fallimenti su altre voci. La si riporta così, senza cambiare seme cercando
    l'esito voluto.
- **Di riferimento, senza esito.**
  - La std a posteriori media di μ_E, contro 0.306 (8a) e 0.328 (toy).
  - Le voci di μ_E sulle 400 esperienze della 8a e della 8c insieme.
- I fallimenti della 8a restano registrati come sono: la 8c li spiega, non li sostituisce.

## 3. Robustezza al prior e troncatura (8a)

| Variante | Atteso | |
|---|---|---|
| Prior uniforme in σ_E, N = 150 | spostamento massimo di μ_E < 0.5 σ; spostamento medio di log σ_E +0.21 ± 0.10 σ (non si dimentica a N = 150); coverage di log σ_E entro 2σ binomiali dal nominale | (c) nel toy |
| *Bordo*: μ_E ~ U(17.5, 19), verità troncata, N = 300 e 1000 | pull e coverage di μ_E, log σ_E, Ω_n con le soglie della §2: il modello contiene la troncatura | (b) |
| *Sborda*: verità non troncata, N = 300 e 1000 | calibrazione rotta, peggio a N = 1000: coverage 68% di μ_E fra 0.38 e 0.55 a N = 1000. Test di sensibilità: si riporta, non conta come fallimento | (c) nel toy, causa (d) |

## 4. Geant4 (8b)

Pool selezionati dei tre run della 7b (+X 799076, +Z 193863, obliqua 222456 eventi), catena e
smearing di D45. Ricampionamento sullo spettro rivelato (D48), dominio [10, 19.5] MeV in tutto.
Modello ideale e con `--acceptance 6acc_s20260908`. Per ogni direzione, toy gemello con le
stesse verità, N e M.

| Voce | +X | +Z | obliqua | |
|---|---|---|---|---|
| Bias di μ_E e log σ_E a N = 1000 | < 0.5 σ a posteriori | < 0.5 σ | < 0.5 σ | (d) |
| Bias vettoriale di Ω̂ a N = 1000 | entro 2σ da 0 | entro 2σ da 0 | < 0.5 σ_c | (a) per simmetria su +X e +Z; (d) sull'obliqua |
| rms G4 / toy gemello, per μ_E, log σ_E, Ω_n | ≤ 1.15 | ≤ 1.15 | ≤ 1.15 | (d) |
| Larghezza del pull e pull rms | 1 ± 0.1 | 1 ± 0.1 | 1 ± 0.1 | (d) |
| Coverage 68/90/95 | entro 3σ binomiali | entro 3σ binomiali | entro 3σ binomiali | (d) |

- Le definizioni sono quelle della 7b (`7_predictions.md` §4.1): bias vettoriale sulla base
  tangente alla verità, valori del toy gemello non arrotondati dal suo record.
- L'eccezione di log σ_E a N = 50 vale anche qui: esito contro il toy gemello, non contro il
  nominale.
- **Modello ideale contro accettanza (d).** Nel Caso A senza accettanza il bias a 14 MeV era
  −0.023 MeV a N = 1000 (fase 6). Qui ci si aspetta un effetto dello stesso ordine su μ_E,
  piccolo rispetto alla deviazione standard a posteriori (~0.09 MeV): entrambi i modelli
  dovrebbero passare. Se fallisce solo il modello ideale, la causa attesa è l'accettanza.
- **Pool.** Con l'accettazione f/(c ĝ) di D48 servono più eventi che nella 7b. Prima dei run,
  in una versione 2 di questo documento, si scrive il numero di esperienze disponibili per riga
  e per direzione, calcolato dal solo pool (nessuna ricostruzione). Se una riga non si riempie
  vale la regola della 7b: soglie con le M disponibili, registrate nel record.
- Un fallimento apre una diagnosi (8c) con una decisione a parte. Le soglie non si correggono.

## 5. Prestazioni (b)

Thread singolo, build `release`: un'esperienza a N = 1000 in ≤ 5 s (D47). Nel toy
un'esperienza a N = 300 costa ~2 s e la checklist ~8 min con 3 processi. I tempi vanno in
`8_results.md`, non in un test con soglia.
