# Risultati della fase 7 — Caso B

Data: 2026-11-25. Previsioni e soglie: `docs/validation/7_predictions.md`, versione 1 in
`09634a8`, committata **prima** del codice del Caso B e dei run. Dopo i run le soglie non sono
state toccate. Questo file copre per ora la parte 7a (port contro toy); la 7b (Geant4) seguirà.

**Esito 7a.**
- La checklist passa **70/76** soglie lette alla lettera.
- Nessuno dei sei fallimenti indica un difetto del port; per ognuno la causa è misurata (§4):
  - due vengono dall'arrotondamento della tabella del toy;
  - due dal punto N = 1000 con M = 40, che fallisce allo stesso modo anche nel toy;
  - due da soglie che il testo lascia ambigue, lette nel modo più severo.
- Il port del Caso B è validato contro il toy; si passa alla 7b.

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
