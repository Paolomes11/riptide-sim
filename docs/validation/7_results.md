# Risultati della fase 7 — Caso B

Data: 2026-11-25. Previsioni e soglie: `docs/validation/7_predictions.md`, versione 1 in
`09634a8`, committata **prima** del codice del Caso B e dei run. Dopo i run le soglie non sono
state toccate. Le precisazioni della versione 2 (§4.1, `4797f09`) sono state scritte prima del
codice e dei run della 7b, senza cambiare numeri. Le §1–5 coprono la 7a (port contro toy), la §6
la 7b (Geant4).

**Esito 7a.**
- La checklist passa **70/76** soglie lette alla lettera.
- Nessuno dei sei fallimenti indica un difetto del port; per ognuno la causa è misurata (§4):
  - due vengono dall'arrotondamento della tabella del toy;
  - due dal punto N = 1000 con M = 40, che fallisce allo stesso modo anche nel toy;
  - due da soglie che il testo lascia ambigue, lette nel modo più severo.
- Il port del Caso B è validato contro il toy; si passa alla 7b.

**Esito 7b** (§6).
- La checklist passa **184/186** soglie lette alla lettera, nelle tre direzioni di D45 con e senza
  accettanza. Bias, risoluzione G4 / toy, copertura e `truth_in_cap` passano ovunque.
- I due fallimenti sono un solo fatto: il pull rms dell'obliqua a N = 300 (M = 100) vale 0.859
  senza accettanza e 0.845 con, cioè −2.0σ e −2.2σ, sulle stesse esperienze.
- Ipotesi (d): fluttuazione dello split. Non è verificata; la verifica è una diagnosi da
  decidere a parte (§6.4).

**Esito 7c** (§6.6, previsioni §4.2, versione 3 in `d7c6222` e `7d4e4c2`).
- Le stesse chiamate dell'obliqua con `--split-seed 20260908` passano **62/62** soglie.
- A N = 300 il pull rms vale 0.974 senza accettanza e 0.971 con, e rms·√N vale 0.234 / 0.233,
  nella banda 0.22–0.27. Per il criterio della §4.2 il fallimento della 7b è una fluttuazione
  dello split **(c)**.

Etichette: **(a)** fatto stabilito; **(b)** stima con ipotesi esplicite; **(c)** verificato con
Monte Carlo; **(d)** ipotesi da mettere alla prova.

## 1. Run

Build `release` a `c45a140` (`git describe` v0.8.0-15-gc45a140, `project_version` 0.8.0).

```
scripts/run_capped.sh -m 5G -- /usr/bin/time -v \
    build/release/analysis/apps/riptide_case_b_calibrate --calibration-id 7a --threads 4
```

- Impostazioni del record:
  - seme `20260907`, σ_Ep = 0.10 MeV, σ_θ = 0.08 rad;
  - dominio di E_n [10, 20] MeV su 500 punti;
  - griglia di 3000 direzioni sulla sfera, calotta di 8000 punti (Δlog = 10);
  - 400 × 801 punti per la tabella θ / θ_obs.
- Prior largo su E_n: [10, 30] MeV su 999 punti.
- Output sull'SSD: `stats/7a/case_b_calibration.json`.
- La valutazione è fatta con script fuori dal repo, che applicano le formule di
  `7_predictions.md` così come sono scritte: `validation/7a_eval/eval_7a.{py,txt}`, con
  `toy_exact_sigma.{py,txt}` e `notes.txt`.
- Un pilot (`stats/7a_pilot/`, conteggi ridotti) ha preceduto il run; serve solo da controllo
  di sanità.

## 2. Port contro toy (§2 delle previsioni)

| N | M | errore medio (°) | rms (°) | risoluzione dichiarata (°) | pull rms | HPD 68/90/95 | rms·√N (rad) |
|---|---|---|---|---|---|---|---|
| 10 | 400 | 4.18 | 4.83 | 4.80 | 0.997 | 0.700 / 0.897 / 0.943 | 0.266 |
| 30 | 300 | 2.21 | 2.52 | 2.56 | 0.982 | 0.700 / 0.887 / 0.957 | 0.241 |
| 100 | 200 | 1.26 | 1.39 | 1.36 | 1.022 | 0.625 / 0.895 / 0.950 | 0.243 |
| 300 | 100 | 0.74 | 0.84 | 0.785 | 1.073 | 0.670 / 0.890 / 0.910 | 0.254 |
| 1000 | 40 | 0.31 | 0.35 | 0.434 | 0.811 | 0.725 / 0.975 / 1.000 | 0.195 |
| 1000 | 240 | 0.38 | 0.43 | 0.433 | 0.984 | 0.700 / 0.900 / 0.946 | 0.235 |

- **Pull rms:** a ogni N entro soglia (1 ± 0.07 per M ≥ 200; 2/√(2M) per M = 100, 40) **(c)**.
- **Copertura HPD:** entro 2σ binomiali a ogni N; a N = 1000 conta il punto M = 240 **(c)**.
- **N = 10:**
  - mediana di pull·√2 1.145 (1.17 ± 0.05);
  - frazione oltre 3: 0.0175 (≤ 0.04; toy 0.018, Rayleigh 0.011) **(c)**.
- **`truth_in_cap`:** 1 in tutte le esperienze di tutte le righe **(c)**.
- **Errore medio, rms, risoluzione dichiarata contro il toy:**
  - entro 2σ a ogni N, salvo la risoluzione a N = 300 e 1000 contro la tabella arrotondata (§4.1);
  - scelta di valutazione: σ della differenza = √2 volte il nostro errore standard (toy con lo
    stesso M); errore di un rms = rms / (2√M), errori gaussiani 2D.
- **rms·√N:** fra 0.22 e 0.27 rad a ogni N, salvo il punto N = 1000 con M = 40 (§4.2).
- **Pendenza log-log su 30–1000:**
  - −0.552 con il punto M = 40, fuori da −0.5 ± 0.05;
  - −0.503 con il punto M = 240 (§4.2).

## 3. Robustezza al prior (§3 delle previsioni)

Stesse esperienze della §2, rigenerate dagli stessi id. Spostamento in unità della risoluzione
dichiarata di riferimento; lo spostamento vettoriale è sulla base (e_θ, e_φ) alla verità.

| N | M | vMF: media / max | E_n [10, 30]: spostamento scalare | E_n: vettore medio (e_θ, e_φ) | E_n: pull rms | E_n: HPD 68/90/95 |
|---|---|---|---|---|---|---|
| 100 | 200 | 0.022 / 0.233 | 0.68 ± 0.02 (max 1.69) | −0.02 ± 0.04, −0.08 ± 0.04 | 0.928 | 0.730 / 0.940 / 0.970 |
| 300 | 100 | 0.013 / 0.324 | 0.69 ± 0.04 (max 1.71) | −0.00 ± 0.06, −0.00 ± 0.05 | 0.979 | 0.690 / 0.880 / 0.960 |
| 1000 | 40 | 0.000 / 0.000 | 0.81 ± 0.07 (max 1.97) | +0.03 ± 0.08, +0.01 ± 0.13 | 0.834 | 0.825 / 0.975 / 0.975 |

- **Prior vMF (κ = 1, asse a 90° dalla verità):**
  - spostamento medio 0.00–0.02 σ, contro una stima di 0.02–0.09 σ: il prior sul parametro
    condiviso si dimentica **(c)**;
  - il massimo supera 0.3 σ in una sola esperienza su 100 a N = 300 (§4.4).
- **Prior largo su E_n:**
  - il vettore medio è entro 2σ da 0 a ogni N, quindi lo spostamento è a media nulla;
  - lo spostamento scalare è 0.68–0.81 σ, come il ~0.75 σ del toy: il prior sul nuisance non si
    dimentica ed è rumore a media nulla **(c)**;
  - copertura entro 2σ binomiali; il pull rms passa a N = 300 e 1000, mentre a N = 100 è 0.928
    contro 0.93 (§4.3);
  - la copertura leggermente sopra il nominale a N = 100 è coerente con un posterior più largo
    sotto un prior più largo del vero **(b)**.

## 4. I sei fallimenti

### 4.1 Risoluzione dichiarata a N = 300 e 1000 (2 voci)

0.785° contro 0.79 e 0.434° contro 0.43, con una tolleranza di ±0.004° e ±0.002°.
- La tabella del toy in `7_predictions.md` è arrotondata a 0.01°, quindi è più grossolana
  dell'errore statistico della media delle σ.
- I valori non arrotondati, ricalcolati dai risultati per esperienza del toy
  (`outputs/caso_B_checklist_results.pkl` del toy), sono 0.7862 ± 0.0017° e 0.4328 ± 0.0010°.
- Con questi valori lo scarto è z = −0.42 e +0.71; a tutti gli N |z| ≤ 1.88
  (`toy_exact_sigma.txt`).
- È un artefatto della tabella, non del port **(c)**.

### 4.2 rms·√N e pendenza con il punto N = 1000, M = 40 (2 voci)

rms·√N = 0.195 rad (soglia 0.22–0.27) e pendenza −0.552 (soglia −0.5 ± 0.05).
- Il resoconto del toy ricava rms·√N = 0.23 a N = 1000 e la pendenza −0.51 su 30–1000 dal punto
  **M = 240**.
- Il punto M = 40 del toy stesso darebbe 0.204 rad e −0.546: fallirebbe allo stesso modo
  (`notes.txt`).
- Con il punto M = 240 la 7a dà 0.235 rad e −0.503, entrambi entro soglia.
- La soglia non dice quale punto usare a N = 1000. Per la copertura lo dice ("con M = 40 conta
  l'esito con M = 240"), per queste due no; qui la voce è riportata come fallita alla lettera.
- È la fluttuazione nota a M = 40, già vista nel toy (rms 0.37° contro 0.42° con M = 240) **(c)**.

### 4.3 Pull rms con il prior largo su E_n a N = 100 (1 voce)

0.928 contro 1 ± 0.07.
- La §3 dice solo "pull e copertura restano calibrati", senza un numero. La soglia applicata è
  quella della §2 per M = 200, ed è la lettura più severa.
- Lo scarto dal limite è 0.002, cioè 0.2σ statistici dell'rms (σ ≈ 1/√(2M) = 0.05).
- La copertura della stessa riga passa.

### 4.4 Spostamento massimo con il prior vMF a N = 300 (1 voce)

0.324 σ in un'esperienza su 100, con una media di 0.013 ± 0.006 σ.
- La soglia "spostamento del MAP < 0.3 σ" non dice se si riferisce a ogni esperienza o alla
  media. Qui è applicata al massimo, la lettura più severa; la stima delle previsioni
  (0.02–0.09 σ) descrive lo spostamento tipico, che passa largamente.
- Media quasi nulla, pochi spostamenti grandi e massimo esattamente 0 a N = 1000: è l'andamento
  di un MAP su una griglia discreta, che resta sullo stesso punto della calotta o salta a un
  punto vicino. Che il salto valga un passo della calotta non è verificato **(d)**.

## 5. Tempi e risorse (§5 delle previsioni)

| Voce | Atteso | Misurato | |
|---|---|---|---|
| Un'esperienza a N = 1000, thread singolo, `release` | ≤ 2 s | 1.88–1.97 s (`riptide_stats_bench --case b`) | PASS **(a)** |
| Checklist 7a completa, 4 thread | ≤ 10 min, RSS < 1 GB | 5 min 11 s, RSS 51 MB | PASS **(a)** |
| Risultati con 1 e 4 thread | identici bit per bit | per costruzione; test Catch2 sulle righe della checklist | PASS **(a)** |

- La §5 conta "~1.0·10⁵ eventi più le varianti di prior": sono le righe della checklist
  (1.03·10⁵ eventi).
- Il run completo contiene in più la ripetizione a N = 1000 con M = 240, cioè 2.4·10⁵ eventi.
  Le varianti di prior rigenerano 0.9·10⁵ eventi e li analizzano due volte.
- In tutto 4.3·10⁵ eventi generati e 5.2·10⁵ analizzati. La ripetizione da sola prende 146 s
  dei 311.

## 6. Geant4 (7b, §4 e §4.1 delle previsioni)

Build `release` a `0eb2fac` (`git describe` v0.8.0-25-g0eb2fac), Geant4 11.4.1, `QGSP_BIC_HP`,
BC-408 60 × 20 × 20 mm. Il campo lontano è piatto in 10–19.5 MeV, 10⁷ eventi, seme 20260907.
Valutazione in `validation/7b_eval/` sulla radice dei dati (`eval_c1.py`, `eval_7b.py`,
`eval_7b.txt`).

### 6.1 Run e verifiche `c1`

| Run | Direzione | p χ² (u / z / frazioni) | max\|pull\| (soglia, test) | a₁ | Protoni che escono | Esito |
|---|---|---|---|---|---|---|
| `6_ff_flat_10_19p5MeV_x_1e7` (fase 6, riuso) | +X | 0.299 / 0.790 / 0.975 | 1.44 (3.291, 7) | −0.01435 | 36727 | PASS **(a)** |
| `7_ff_flat_10_19p5MeV_z_1e7` | +Z | 0.770 / 0.250 / 0.520 | 1.69 (3.291, 7) | −0.01403 | 12503 | PASS **(a)** |
| `7_ff_flat_10_19p5MeV_obl_1e7` | (3, 2, 1)/√14 | 0.415 / 0.184 / 0.730 | 1.31 (3.317, 8) | −0.01398 | 12328 | PASS **(a)** |

- Si applicano le regole `c1` della fase 6 (§1 di `6_predictions.md`). Gli oracoli
  deterministici sono a zero o entro le soglie esatte: contatori, bilancio dei depositi, kB,
  `observed_non_hp`, offset HP − fisico.
- **Corda.** Lungo +Z la corda è unica: 20 mm entro 10⁻⁹ mm. Lungo l'obliqua la corda varia da
  evento a evento, con un massimo di 37.42 mm, sotto la diagonale del blocco (66.33 mm). Il
  pull della corda media entra fra i test: 8 pull e soglia 3.317 dal `test_plan` della
  previsione.
- **Precisazione (a).** Le attese di `riptide_sim_predict` per +Z e per l'obliqua non erano
  committate prima dei run. Sono deterministiche, ottenute dal config e dai dati di riferimento
  senza leggere il run: per l'obliqua sono state generate prima della simulazione, per +Z dopo.
  Sono copiate in `predictions/7_ff_*.json`.

### 6.2 Selezione e righe

Catena `single_np,has_proton,contained`, conversione HP, smearing D31 (σ_Ep 0.10, σ_θ 0.08 rad).

| Direzione | Selezionati | Ripetizione a N = 1000 |
|---|---|---|
| +X | 799076 | M = 240 |
| +Z | 193863 | M = 153 di 240 |
| obliqua | 222456 | M = 182 di 240 |

- **Precisazione (a).** La §4.1 prevede la ripetizione (1000, 240). Lungo +Z e l'obliqua la
  sezione trasversale è minore, quindi il pool non basta: la riga M = 40 e la ripetizione
  richiedono 2.8·10⁵ eventi disgiunti. Il split prende tutte le esperienze disponibili e il
  record lo registra (`experiments` < `experiments_planned`). Il toy gemello ha lo stesso pool
  e quindi gli stessi M.
- Si applicano le soglie con le M disponibili. Per N = 1000 conta la ripetizione, come nella
  §4.1.

### 6.3 Esiti per voce

Ordine: bias → risoluzione → pull → copertura → contrazione. "id" è il record `stats/7b_<id>`.
Senza accettanza e con `--acceptance 6acc_s20260908`.

**Bias vettoriale a N = 1000** (z = componente / errore, nella base (e_θ, e_φ) alla verità):

| id | M | z(e_θ), z(e_φ) | norma | σ_c | Soglia | Esito |
|---|---|---|---|---|---|---|
| `x` | 240 | −0.73, +0.19 | 2.5·10⁻⁴ rad | 7.71·10⁻³ | \|z\| ≤ 2 | PASS **(a)** |
| `x_acc` | 240 | −0.75, +0.17 | 2.6·10⁻⁴ | 7.71·10⁻³ | \|z\| ≤ 2 | PASS **(a)** |
| `z` | 153 | +1.10, +0.74 | 5.9·10⁻⁴ | 7.50·10⁻³ | \|z\| ≤ 2 | PASS **(a)** |
| `z_acc` | 153 | +1.13, +0.90 | 6.4·10⁻⁴ | 7.50·10⁻³ | \|z\| ≤ 2 | PASS **(a)** |
| `obl` | 182 | +1.30, −0.17 | 5.0·10⁻⁴ | 7.59·10⁻³ | norma < 0.5 σ_c | PASS **(a)**, 0.066 σ_c |
| `obl_acc` | 182 | +1.31, −0.05 | 5.0·10⁻⁴ | 7.59·10⁻³ | norma < 0.5 σ_c | PASS **(a)**, 0.066 σ_c |

Il bias dell'obliqua è sotto il 7% della risoluzione per esperienza: l'ipotesi (d) della §4
regge. Anche espresso come pull per componente, l'obliqua è entro 2σ.

**Risoluzione G4 / toy (≤ 1.15), pull rms (1 ± 0.1), copertura HPD 68/90/95 (3σ binomiali)**, a
ogni N. L'errore sul rapporto delle risoluzioni è ≈ r·√(1/2M + 1/2M) **(b)**.

| id | Rapporto rms, intervallo su N | Rapporto σ dichiarate | Pull rms, intervallo | Copertura | Esito (25 voci) |
|---|---|---|---|---|---|
| `x` | 0.92–1.12 | 0.99–1.00 | 0.955–1.020 | 15/15 | 25/25 **(a)** |
| `x_acc` | 0.97–1.07 | 1.00–1.01 | 0.954–1.016 | 15/15 | 25/25 **(a)** |
| `z` | 0.96–1.09 | 0.95–0.98 | 0.993–1.049 | 15/15 | 25/25 **(a)** |
| `z_acc` | 1.01–1.09 | 0.95–0.99 | 0.997–1.052 | 15/15 | 25/25 **(a)** |
| `obl` | 0.91–1.04 | 0.97–0.99 | 0.859–1.020 | 15/15 | 24/25 **(a)** |
| `obl_acc` | 0.86–1.04 | 0.97–1.00 | 0.845–1.017 | 15/15 | 24/25 **(a)** |

- Il caso più vicino alla soglia di risoluzione è `x` a N = 30: 1.116 ± 0.064.
- `truth_in_cap` vale 1 in tutte le esperienze di tutte le righe **(a)**.
- Il rapporto delle σ dichiarate resta fra 0.95 e 1.01: Geant4 e toy dichiarano la stessa
  risoluzione. Lungo +Z quella di Geant4 è più stretta del 2–5%, pur con errori realizzati
  uguali entro la statistica (riportato, senza soglia).

**Contrazione e coda a N = 10** (riportate, senza soglia per Geant4 nella §4):

| id | Pendenza 30–1000, con la ripetizione | Toy | rms·√N a N = 1000 (toy) | Mediana di pull·√2 a N = 10 (toy) | Coda > 3 (toy) |
|---|---|---|---|---|---|
| `x` | −0.514 | −0.470 | 0.235 (0.240) | 1.115 (1.148) | 0.0125 (0.0075) |
| `x_acc` | −0.512 | −0.486 | 0.236 (0.238) | 1.114 (1.138) | 0.0125 (0.0050) |
| `z` | −0.489 | −0.488 | 0.249 (0.241) | 1.191 (1.126) | 0.0150 (0.0125) |
| `z_acc` | −0.490 | −0.507 | 0.249 (0.229) | 1.191 (1.071) | 0.0150 (0.0100) |
| `obl` | −0.517 | −0.521 | 0.228 (0.229) | 1.165 (1.079) | 0.0175 (0.0100) |
| `obl_acc` | −0.519 | −0.512 | 0.228 (0.219) | 1.163 (1.198) | 0.0150 (0.0025) |

- Con la ripetizione la pendenza di Geant4 è −0.5 entro ±0.02 in tutte e sei le combinazioni;
  quella del toy va da −0.470 a −0.521 **(c)**.
- rms·√N di Geant4 cade nella banda 0.22–0.27 rad della 7a (§2) a ogni N, tranne la riga
  dell'obliqua a N = 300 (0.202–0.206), che è quella del fallimento (§6.4) **(c)**.
- La riga N = 1000, M = 40 è riportata solo nei record e non entra negli esiti.

### 6.4 Il fallimento: pull rms dell'obliqua a N = 300

| id | Pull rms | Toy | Copertura 68/90/95 | rms·√N (toy) | Rapporto σ |
|---|---|---|---|---|---|
| `obl` | 0.859 | 0.929 | 0.78 / 0.96 / 0.99 | 0.206 (0.225) | 0.985 |
| `obl_acc` | 0.845 | 0.971 | 0.80 / 0.96 / 0.99 | 0.202 (0.234) | 0.989 |

- Le due righe usano le stesse 100 esperienze sugli stessi eventi, quindi sono un solo
  fallimento, non due indipendenti **(a)**.
- Con M = 100 l'errore statistico del pull rms è ≈ 1/√(2M) = 0.071 **(b)**. Gli scarti sono
  quindi −2.0σ e −2.2σ.
- La soglia ±0.1 vale 1.41σ a M = 100, per cui una riga corretta fallisce con probabilità
  ≈ 16% **(b)**. Sulle 30 voci di pull di Geant4, quelle indipendenti sono circa 15, perché le
  righe con e senza accettanza condividono le esperienze. Sommando le probabilità di falso
  allarme riga per riga, con M da 100 a 400, si attendono ≈ 0.8 fallimenti, e almeno uno con
  probabilità ≈ 60% **(b)**.
- Il quadro è coerente con un insieme di errori realizzati più piccoli del solito, a σ
  dichiarata invariata:
  - rms·√N scende a 0.202–0.206, contro 0.225–0.228 a N = 100 e 1000 nella stessa direzione;
  - la copertura sale su tutti e tre i livelli, pur entro 3σ;
  - la σ dichiarata è quella del toy entro l'1.5%.
- **Ipotesi (d):** fluttuazione statistica dello split, non un difetto del modello. Non è
  verificata. La verifica naturale è la stessa riga con un `--split-seed` diverso, oppure su un
  run indipendente: è una diagnosi, da decidere a parte. Le soglie non cambiano.

### 6.5 Tempi e risorse

| Voce | Misurato |
|---|---|
| Simulazione 10⁷ eventi, 8 thread (+Z) | 156 s |
| Pilot `--pilot` su +X (799 076 eventi) | 31 s, RSS 1.4 GB |
| Chiamata completa, 4 thread (somma delle righe) | +X 211–227 s, +Z 158 s, obliqua 176 s |
| Dodici chiamate in sequenza | ≈ 38 min |

- L'RSS di 1.4 GB supera la stima "< 1 GB" del piano (non una soglia delle previsioni): è
  l'intero pool di 10⁷ eventi letto in memoria prima dei tagli. Resta sotto il tetto di 5 GB.
- Le chiamate toy durano come quelle Geant4: il tempo è nell'analisi delle righe, non nella
  lettura.

### 6.6 Diagnosi 7c: l'obliqua con un altro split

Build `release` a `7d4e4c2`, stesso run `7_ff_flat_10_19p5MeV_obl_1e7` e stesso `--seed`
20260907; cambia solo `--split-seed` (20260908). Record `stats/7c_obl{,_acc}{,_toy}`,
valutazione in `validation/7c_eval/` sulla radice dei dati (`run_7c.sh`, `eval_7c.py`,
`eval_7c.txt`).

**Esito per voce.** 62/62 **(a)**. Ogni chiamata Geant4 ha 31 voci: le 25 della tabella della
§6.3, più `truth_in_cap` per le cinque righe e il bias. Le M sono quelle della 7b: (10, 400),
(30, 300), (100, 200), (300, 100), ripetizione (1000, 182).

**Il criterio della §4.2**, a N = 300, M = 100:

| Record | Split | Pull rms | rms·√N (rad) | Copertura 68/90/95 |
|---|---|---|---|---|
| `7b_obl` | 20260907 | 0.859 | 0.206 | 0.78 / 0.96 / 0.99 |
| `7c_obl` | 20260908 | 0.974 | 0.234 | 0.68 / 0.89 / 0.95 |
| `7b_obl_acc` | 20260907 | 0.845 | 0.202 | 0.80 / 0.96 / 0.99 |
| `7c_obl_acc` | 20260908 | 0.971 | 0.233 | 0.68 / 0.89 / 0.96 |

- Il pull rms a N = 300 cade in 1 ± 0.1 in tutte e due le chiamate, e rms·√N torna nella banda
  della 7a. L'esito è "fluttuazione": l'ipotesi (d) della §6.4 passa a **(c)**.
- La σ dichiarata resta quella della 7b entro l'1% (rapporto G4 / toy 0.988 / 0.993): cambia
  l'insieme degli errori realizzati, non il modello.
- Limite dichiarato nella §4.2: il pool è lo stesso. Un difetto legato a questo campione di
  eventi non si vedrebbe con questo test.

**Altre voci, riportate.**

| id | Rapporto rms, intervallo su N | Pull rms, intervallo | Pendenza con la ripetizione (toy) |
|---|---|---|---|
| `7c_obl` | 0.925–1.053 | 0.906–1.020 | −0.537 (−0.505) |
| `7c_obl_acc` | 0.886–1.051 | 0.906–1.016 | −0.537 (−0.507) |

- Il pull rms più basso è ora quello della ripetizione a N = 1000: 0.906, dentro la soglia di
  0.094. Nella 7b la stessa riga dava 0.951 **(a)**.
- **Bias a N = 1000 (d).** La norma vale 8.5·10⁻⁴ / 7.9·10⁻⁴ rad, cioè 0.11 / 0.10 σ_c, sotto la
  soglia di 0.5 σ_c. La componente e_θ è positiva in tutti e due gli split: z = +1.30 nella 7b e
  +2.19 / +2.05 nella 7c. Il toy gemello non la mostra (norma ≤ 4.1·10⁻⁴, componenti negative).
  - Gli split usano lo stesso pool, quindi i due valori non sono indipendenti.
  - Se l'effetto è reale, è il candidato atteso da D45: la parte della selezione che dipende
    dalla direzione. Vale circa un decimo della risoluzione per esperienza a N = 1000.
  - Non è un fallimento e non apre una diagnosi adesso. Va rivisto se una fase successiva (Caso
    C, N più grandi, fase 9) renderà rilevante un bias di questa scala.
- Tempi: quattro chiamate in sequenza, 4 thread, nello stesso ordine di grandezza della 7b
  (≈ 3 min ciascuna); il tetto di 5 GB non è stato raggiunto.
