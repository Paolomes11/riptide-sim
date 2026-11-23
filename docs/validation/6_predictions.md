# Previsioni della fase 6 — scritte prima dei run

Versione 1 (2026-11-23). Va committata prima dei run di simulazione della fase 6 e prima di
qualsiasi esecuzione di `riptide_case_a_geant4` sui run Geant4. I risultati (`6_results.md`)
citeranno l'hash di quel commit. Dopo i run le soglie non si toccano.

La fase 6 applica il Caso A (D32, D33), validato sul toy nella fase 5, a eventi Geant4 lungo +X
(D36). `riptide_case_a_geant4` legge un run, applica la catena `single_np,has_proton,contained`,
la conversione HP (D24) e lo smearing D31 (σ_Ep = 0.10 MeV, σ_θ = 0.08 rad). Poi divide gli
eventi selezionati in esperienze disgiunte (D37) e scrive `stats/<id>/case_a.json`. Con `--toy`
la stessa app genera gli eventi dal modello del toy con la stessa verità (D38). La differenza
Geant4 − toy, voce per voce, è il costo del realismo del rivelatore.

Fonti dei numeri:
- le attese di `riptide_sim_predict` per i quattro config nuovi, generate dalla build `release`
  di `6055ab1` e copiate in `predictions/6_ff_*.json`;
- per il 14 MeV +X, il run `3_ff_mono_14MeV_x_1e7` con le sue verifiche `c1` (`3_results.md`) e
  la selezione della fase 4 (`4_results.md`);
- la calibrazione della fase 5 (`5_results.md`) e il controllo di coerenza dell'app: `--toy` con
  verità `flat:10:20` e 20 000 eventi ridà la checklist della §1 di `5_results.md` numero per
  numero (c);
- le ipotesi D39.

Etichette: **(a)** fatto stabilito (sorgenti, dati, costruzione del codice); **(b)** stima con
ipotesi esplicite; **(c)** verificato con Monte Carlo; **(d)** ipotesi da mettere alla prova.

## Run previsti

Build `release` del commit di questo file, Geant4 11.4.1, `QGSP_BIC_HP`, G4NDL 4.7.1, blocco di
BC-408 da 60 × 20 × 20 mm, direzione +X (D27, D36). Seme del progetto `20260907`. Ogni comando
gira sotto `scripts/run_capped.sh -m 5G`, uno alla volta, dopo `scripts/data_root.sh check`.

### Simulazione e verifiche `c1`

| Run (`run_id`) | Config | Eventi | Thread |
|---|---|---|---|
| `6_ff_flat_10_19p5MeV_x_1e7` | `6_ff_flat_10_19p5MeV_x.json` | 10⁷ | 8 |
| `6_ff_mono_11MeV_x_1e6` | `6_ff_mono_11MeV_x.json` | 10⁶ | 8 |
| `6_ff_mono_17MeV_x_1e6` | `6_ff_mono_17MeV_x.json` | 10⁶ | 8 |
| `6_ff_mono_19MeV_x_1e6` | `6_ff_mono_19MeV_x.json` | 10⁶ | 8 |
| `3_ff_mono_14MeV_x_1e7` | esiste già (fase 3), non si rilancia | 10⁷ | 8 |

```bash
scripts/run_capped.sh -m 5G -- build/release/sim/riptide_sim --config config/<c>.json \
  --events <n> --threads 8 --run-id <run_id>
scripts/run_capped.sh -m 5G -- build/release/sim/riptide_sim_check --run-id <run_id> \
  --check-id c1 --bins-u 50 --bins-z 50
```

Le attese si rigenerano con:

```bash
build/release/sim/riptide_sim_predict --config config/<c>.json --events <n> \
  --bins-u 50 --bins-z 50
```

### Caso A, Geant4 e toy gemello

Per ogni run, nell'ordine:

1. il pilota;
2. il run completo;
3. il toy gemello.

`--events` del gemello è `source.events_selected` del `case_a.json` Geant4 completo (D38). La
verità del gemello è quella del run: `mono:<E>`, oppure `flat:10:19.5` per il piatto.

```bash
B=build/release/analysis/apps/riptide_case_a_geant4
scripts/run_capped.sh -m 5G -- $B --stats-id 6_<run_id>_pilot --run-id <run_id> --pilot
scripts/run_capped.sh -m 5G -- $B --stats-id 6_<run_id> --run-id <run_id>
scripts/run_capped.sh -m 5G -- $B --stats-id 6_<run_id>_toy --toy --truth <spec> \
  --events <events_selected>
```

Il 14 MeV ha `run.json` v1, quindi i suoi due comandi Geant4 prendono anche
`--hp-from-check c1` (a). Seme dello smearing e del toy, e seme dello split: i default
`20260907`.

**Criterio pilota → completo.** Il completo parte se il pilota finisce senza errori, con RSS
< 3 GB e tempo < 60 s. Altrimenti si torna al piano.

## 1. Verifiche `c1` dei run nuovi

Valgono le stesse regole della fase 3 per la direzione +X (`3_predictions.md`), con le attese dei
nuovi file di previsione. Un run entra nell'analisi solo se le passa.

- **Oracoli deterministici, §2 della fase 3 (soglie esatte):**
  - contatori `must_be_zero` a zero;
  - corda di 60 mm entro 10⁻⁹ mm;
  - bilancio dei depositi entro [−10⁻⁹, 0.01] MeV;
  - kB 0.155;
  - `model_mix.observed_non_hp` = 0 (D25; `overlaps_model_window` è `false` per tutti e
    quattro);
  - offset HP − fisico nell'intervallo di `prediction.Ep_physical_offset_relative`, con margine
    relativo 10⁻⁹:
    - piatto: [5.33500, 5.50179]·10⁻⁴;
    - 11 MeV: [5.38277, 5.47756]·10⁻⁴;
    - 17 MeV: [5.34896, 5.49469]·10⁻⁴;
    - 19 MeV: [5.33779, 5.50037]·10⁻⁴.
- **Test statistici contro `jeff33`, §3 della fase 3:** 3 χ² e 7 pull per run (a₁…a₆ e
  accettanza), con α di famiglia 0.01, α per test 10⁻³ e soglia |pull| < 3.291 (`test_plan`).
- **Senza soglia, riportati in `6_results.md`:**
  - il confronto con PSTAR (le bande della fase 3 sono state fissate solo per 14 MeV e per il
    piatto isotropo);
  - i protoni che escono.

Attese principali (a), da `libraries[jeff33]`:

| Run | P(prima interazione n–p) | a₁ atteso | σ(a₁) |
|---|---|---|---|
| piatto 10–19.5 | 0.10117 | −1.448·10⁻² | 5.8·10⁻⁴ |
| 11 MeV | 0.12563 | −1.255·10⁻² | 1.63·10⁻³ |
| 14 MeV (fase 3) | 0.10438 | −1.48·10⁻² | — |
| 17 MeV | 0.08685 | −1.588·10⁻² | 1.97·10⁻³ |
| 19 MeV | 0.07744 | −1.609·10⁻² | 2.09·10⁻³ |

L'accettanza geometrica vale 0.636620 per tutti (D27).

## 2. Selezione e numero di esperienze

**Eventi selezionati (b).** Per il 14 MeV il numero è noto (a): 826 859 su 10⁷. Diviso per gli
eventi con prima interazione n–p attesi (1.0438·10⁶), dà s = 0.7922. s raccoglie le perdite di
`single_np` e di `contained`. Se s resta costante con E_n:

| Run | Eventi selezionati attesi | Banda |
|---|---|---|
| piatto 10–19.5 | 8.01·10⁵ | ±10% |
| 11 MeV | 9.95·10⁴ | ±10% |
| 17 MeV | 6.88·10⁴ | ±10% |
| 19 MeV | 6.13·10⁴ | ±10% |

La banda di ±10% copre la variazione di s con E_n (D39). Il verso non si prevede, perché due
effetti si oppongono:
- a energia alta escono più protoni, quindi `contained` toglie di più;
- a energia bassa il neutrone diffuso interagisce di più, quindi `single_np` toglie di più.

Un valore fuori banda non invalida il run, ma va spiegato in `6_results.md`.

**M per riga (a).** Le righe sono quelle di D37:
- contrazione: N/M = 10/400, 100/200 e 1000/40;
- prior su θ_p: 1000/20;
- prior su E_n: 1/400, 3/400, 10/400, 30/200 e 100/200;
- checklist del piatto: 20 000 eventi.

Ogni riga ha il suo split. La più esigente vuole 1000 · 40 = 40 000 eventi, meno del minimo
atteso (6.1·10⁴, anche al bordo basso della banda). Quindi **nessuna riga si riduce**: in ogni
`case_a.json` completo vale `experiments` = `experiments_planned` e `skipped` = `null`. Nel
pilota M è diviso per 10, con lo stesso N.

## 3. Prove multiple e regole di confronto

**Confronto Geant4 − toy.** Per ogni voce Δ = X_G4 − X_toy, con errore combinato
σ_c = √(σ_G4² + σ_toy²). I due campioni sono indipendenti. Gli errori di una voce sono:
- `bias_mean`: `bias_error`;
- `std_sqrt_N`: `std_sqrt_N_error`;
- `shift_mean`: `shift_error`;
- `pull_mean`: `pull_width`/√M;
- `pull_width`: `pull_width`/√(2M);
- `coverage`: √(p(1−p)/M), con p nominale;
- bias per bin: `rms_error`/√`counts`, conservativo come nella fase 5.

Le voci della famiglia «Geant4 = toy» sono K = 117:
- 5 bin del bias del piatto;
- per ognuno dei quattro mono, 28 voci:
  - contrazione: 3 righe × 7 voci (`bias_mean`, `std_sqrt_N`, `pull_mean`, `pull_width` e tre
    `coverage`);
  - prior su θ_p: `flat.bias_mean` e `physical.bias_mean`;
  - prior su E_n: 5 righe di `shift_mean`.

Con α di famiglia 0.01 e Bonferroni la soglia è **|Δ| ≤ 3.93·σ_c**. Le voci della stessa riga
sono correlate, quindi la correzione è conservativa (b).

**Lettura.**
- Una voce entro soglia vuol dire: costo del realismo non risolto con questa statistica.
- Una voce oltre soglia è un **costo misurato**, non un errore del codice. Va riportata con
  grandezza e verso, e confrontata con le stime (b) delle sezioni seguenti.
  - Se il verso è quello previsto e la grandezza è entro un fattore 2, è spiegata.
  - Altrimenti è un'anomalia da indagare prima della chiusura della fase.

**Calibrazione assoluta.** È una seconda famiglia, separata, per Geant4 e per il toy. Per ogni
mono, alle righe di contrazione N = 100 e 1000:
- `pull_mean` = 0 ± z·`pull_width`/√M;
- le tre `coverage` entro z errori binomiali del nominale.

Le voci sono 8 per mono e per sorgente, cioè K = 32 per sorgente. Con α di famiglia 0.01 la
soglia è z = 3.60. Atteso (b): a N ≥ 100 il prior su E_n non conta, anche a 11 e 19 MeV, perché
la larghezza combinata (≲ 0.1 MeV) sta lontano dai bordi 10 e 20. Quindi il toy passa tutte le
32 voci. Una voce del toy che fallisce indica un difetto dell'app, non del rivelatore. Una voce
Geant4 che fallisce con il toy che passa rende il toy **circolare** su quella voce (D38).

**Fluttuazioni.** Come nella fase 5: un fallimento isolato fra la soglia e 4σ_c è un candidato
fluttuazione. Come diagnosi, e non in sostituzione, si ripete il solo `case_a` di quel run con
`--seed 20260908 --split-seed 20260908`, Geant4 e gemello, e si riportano entrambi gli esiti.

## 4. Checklist sul piatto (`checklist`)

Sono 20 000 eventi selezionati a caso, con verità per evento. Il posterior per evento usa il
prior piatto su [10, 20] e sin 2θ_p.

**Spettro selezionato non uniforme (b).** Nel piatto, gli eventi selezionati hanno E_n
distribuita come P_np(E)·s(E), non uniforme (D39). P_np scende da 0.1256 a 11 MeV a 0.0774 a
19 MeV, cioè d ln P_np/dE ≈ −0.060/MeV (a, dalle attese). Con s costante, i 5 bin uguali su
[10, 19.5] contengono circa 24.8 / 22.1 / 19.7 / 17.6 / 15.7% degli eventi, contro il 20% del
toy.

**Bias per bin (d): Geant4 = toy** entro la soglia della §3. Per singolo evento l'effetto
dell'accettanza è ≈ g·σ₁², dove:
- g = d ln s/dE, con |g| ≲ 0.0125/MeV per la banda della §2;
- σ₁ ≈ 1.9 MeV è la larghezza di un evento.

Quindi l'effetto è ≲ 0.05 MeV, contro una soglia di ~0.17 MeV per bin (σ_c ≈ 0.042 MeV, con
`rms_error` ≈ 1.9 MeV e ~4000 eventi per bin). I bin dei due campioni coincidono entro 10⁻³ MeV,
perché entrambi coprono le verità [10, 19.5] (a).

**Bias globale, previsione con formula fissata (b).** Il bias medio del piatto Geant4 si
confronta con il toy ripesato:

  Δ_atteso = Σ_k (f_k^G4 − f_k^toy) · b_k^toy,

dove:
- f_k è la frazione di eventi nel bin k (`counts`/Σ`counts`) di ciascun campione;
- b_k^toy è il bias per bin del toy.

Soglia: |(`bias_mean`_G4 − `bias_mean`_toy) − Δ_atteso| ≤ 3.93·σ_c + 0.03 MeV. Lo 0.03 MeV copre
la pendenza del bias dentro un bin. Stima prima del run, con i bias per bin del toy
`flat:10:20`: **Δ_atteso ≈ +0.19 MeV**, e da +0.15 a +0.23 MeV per |g| ≤ 0.0125/MeV. È positivo
perché il piatto Geant4 ha più eventi nel bin basso, dove il prior spinge le stime verso l'alto.
È ~15 volte `bias_error`: il bias globale del piatto Geant4 **non** è zero, e non va letto come
un difetto.

**Pull, larghezza e coverage globali:** senza soglia. La verità non viene dal prior, quindi la
calibrazione esatta della fase 5 non vale (D39), e il toy gemello non è un riferimento diretto.
Si riportano Geant4 e toy affiancati. Atteso (b): la coverage di Geant4 è entro ~0.01 da quella
del toy.

## 5. Contrazione (`contraction`) sui quattro mono

**Geant4 = toy** su tutte le 21 voci per mono (d).

**Bias da accettanza (b).** Il modello non ha il fattore di selezione A(E) = s(E). Con N eventi,
il massimo del combinato si sposta di ≈ g·σ₁², dove σ₁ = `std_sqrt_N` ≈ 0.93–0.98 MeV (fase 5).
Lo spostamento non dipende da N, mentre il pull cresce come √N.
- Con |g| ≲ 0.0125/MeV: |Δ `bias_mean`| ≲ 0.011 MeV.
- A N = 1000: Δ `pull_mean` ≲ 0.4, contro una soglia di 3.93·√2/√40 ≈ 0.88.

Quindi **nessuna voce della contrazione oltre soglia** (d). Il bias a N = 1000 si confronta dopo
il run con la formula fissata qui:
- s(E) misurato = `events_selected` / (eventi · P_np) sui quattro mono, con il 14 MeV da 10⁷;
- g(E) da differenze finite centrate (a 11 e 19 MeV, quelle unilaterali);
- Δ `bias_mean`_atteso = g · `std_sqrt_N`².

Il confronto è riportato, ma non è una soglia, perché l'errore di Δ a N = 1000 (~0.007 MeV) è
dello stesso ordine dell'effetto.

**Anisotropia in CM (b, d).** κ = 3⟨P₁⟩ ≈ −0.044, dai valori di a₁ della §1. La riponderazione
in cos θ_CM sposta il combinato di ≲ |κ|·0.1 MeV ≈ 0.005 MeV. Il verso non si prevede: la
riponderazione è antisimmetrica attorno a θ_p = 45°, e il suo effetto su E_n dipende dalla forma
della verosimiglianza. È sotto la sensibilità.

**Contrazione 1/√N (d):** per Geant4 come per il toy, `std_sqrt_N` decresce da N = 10 a 100 e
cambia meno del 10% da 100 a 1000.

**Bordi del prior.** A 11 e 19 MeV e N = 10, il prior attira il combinato verso il centro. È lo
stesso nei due gemelli, quindi si cancella in Δ (a, stessa verità e stesso prior).

## 6. Prior su θ_p (`theta_prior`)

Sono 20 esperimenti di N = 1000 per ogni mono, gli stessi eventi con sin 2θ_p e con il prior
piatto. Il riferimento di D37 è il 14 MeV; gli altri tre si riportano con le stesse soglie.
- **`flat.bias_mean` e `physical.bias_mean`: Geant4 = toy** (d). Lo spostamento da accettanza
  (≲ 0.011 MeV) sta sotto σ_c ≈ √2·0.008 MeV · 3.93 ≈ 0.045 MeV.
- Atteso dal toy della fase 5 a 15 MeV (c): il prior piatto dà un bias di circa −0.10 MeV,
  quello fisico circa 0. Il gemello misura i valori a 11, 14, 17 e 19 MeV.

## 7. Prior su E_n (`energy_prior`)

Il prior è piatto su [10, 20] contro [10, 30], alle righe N = 1, 3, 10, 30 e 100, con la verità
fissa del mono. Nella fase 5 la verità era uniforme, quindi quei `shift_mean` non sono un
riferimento.

- **`shift_mean`: Geant4 = toy** su ogni riga (d). La selezione cambia la verosimiglianza di
  pochi centesimi di MeV, mentre lo spostamento a N piccolo vale decimi di MeV.
- A 19 MeV lo spostamento a N = 1 è il più grande dei quattro mono (b): il prior largo lascia
  andare il posterior sopra 20 MeV. A 11 MeV è il più piccolo.
- `shift_mean` decresce con N per ogni mono e per entrambe le sorgenti (d).
- Il bordo del prior a 20 MeV, contro la sorgente piatta a 19.5, non entra qui, perché i mono
  hanno una sola energia (D39).

## 8. Tempi, memoria e disco (b)

| Passo | Tempo | RSS | Disco |
|---|---|---|---|
| sim piatto 10⁷, 8 thread | ~150–170 s (fase 3) | < 0.2 GB | ~1.2 GB |
| sim mono 10⁶, 8 thread | ~20 s | < 0.2 GB | ~0.12 GB |
| `c1` su 10⁷ / 10⁶ | ~45 s / ~5 s | < 0.2 GB | trascurabile |
| `case_a` piatto, Geant4 | ~15 s (lettura di 10⁷ ~9 s, fase 4, e 20 000 posterior) | ~1.4 GB (come `prepare` in fase 4) | < 1 MB |
| `case_a` mono 10⁶, Geant4 | ~40 s (lavoro per riga come la calibrazione della fase 5, 39.6 s) | < 0.3 GB | < 1 MB |
| `case_a` 14 MeV 10⁷, Geant4 | ~50 s | ~1.4 GB | < 1 MB |
| toy gemelli | come il Geant4 corrispondente, senza lettura | < 0.1 GB | < 1 MB |

Il disco totale è ~1.6 GB, contro 1.6 TB liberi sull'SSD. I tempi entrano in `6_results.md`
senza soglia.

## 9. Riepilogo delle previsioni

| Previsione | Etichetta |
|---|---|
| I quattro run nuovi passano `c1` | (d) |
| Eventi selezionati entro ±10% della §2; nessuna riga ridotta | (b) |
| Bias per bin del piatto: Geant4 = toy | (d) |
| Bias globale del piatto: Geant4 − toy ≈ +0.19 MeV, secondo la formula della §4 | (b) |
| Contrazione, θ_p e prior su E_n: le 112 voci Geant4 = toy | (d) |
| Calibrazione assoluta a N ≥ 100: toy e Geant4 passano | (b) per il toy, (d) per Geant4 |
| Costo del realismo nel Caso A lungo +X sotto la sensibilità di questi campioni, tranne la forma dello spettro selezionato del piatto | (d) |
