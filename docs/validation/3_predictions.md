# Previsioni della fase 3 — scritte prima dei run

Versione 2 (2026-11-21), scritta anche questa prima dei run. Rispetto alla versione 1 (`8429cc3`)
cambia solo il riferimento della luce per bin (§5), che ora è nel repo; le soglie restano quelle.

La geometria (blocco 60 × 20 × 20 mm) e la sorgente a campo lontano sono
descritte in D27. Le verifiche per evento sono quelle di `checks_version` 10
(`docs/architecture/data_layout.md`). Questo file va committato prima dei tre run di validazione;
i risultati (`3_results.md`) citeranno l'hash di quel commit. Dopo i run le soglie non si toccano.

Le fonti dei numeri sono tre:
- le attese di `riptide_sim_predict` per le tre config a 10⁷ eventi, generate dal commit
  `a15d899` e copiate in `predictions/3_ff_*_1e7.json`;
- i run da 10⁷ eventi della fase 2 (`2_results.md`), per ciò che la geometria non cambia (luce
  contro PSTAR per bin di E_p, tassi per protone);
- tre piloti da 10⁵ eventi a 1 thread, uno per config, con la build `release` di `a15d899`
  (`pilot3_*_1e5_1t`). I piloti non sono run di validazione: servono a stimare ciò che non ha
  un'attesa calcolabile (protoni che escono, classi dell'evento, range contro PSTAR), e i loro
  test non entrano nei risultati. Sono una misura (c) di questo codice su un campione con lo
  stesso seme dei run di validazione ma 100 volte più piccolo, quindi non indipendente dai primi
  10⁵ eventi di quei run.

Etichette: **(a)** fatto stabilito (sorgenti, dati, costruzione del codice); **(b)** stima con
ipotesi esplicite; **(c)** verificato con Monte Carlo; **(d)** ipotesi da mettere alla prova.

## Run previsti

Build `release`, Geant4 11.4.1, `QGSP_BIC_HP`, G4NDL 4.7.1, blocco di BC-408 da 60 × 20 × 20 mm
centrato nell'origine con il lato lungo su X (D27), seme del progetto `20260907`. Tutte e tre le
config hanno `config_version` 3, schema 4, sorgente `far_field` con margine 10 mm, kB = 0.155
mm/MeV (D18) e funzione di passo del protone (0.2, 0.1 mm) (D21). La sorgente piatta si ferma a
19.5 MeV, sotto la finestra HP/CHIPS (D25). Ogni run gira sotto `scripts/run_capped.sh -m 5G`,
uno alla volta, dopo `scripts/data_root.sh check`.

| Run (`run_id`) | Config | Eventi | Thread | Ruolo |
|---|---|---|---|---|
| `3_ff_mono_14MeV_x_1e7` | `3_ff_mono_14MeV_x.json` | 10⁷ | 8 | direzione fissa +X (Caso A) |
| `3_ff_mono_14MeV_isotropic_1e7` | `3_ff_mono_14MeV_isotropic.json` | 10⁷ | 8 | direzioni isotrope, mono |
| `3_ff_flat_10_19p5MeV_isotropic_1e7` | `3_ff_flat_10_19p5MeV_isotropic.json` | 10⁷ | 8 | isotropa, sorgente di produzione |

```bash
scripts/run_capped.sh -m 5G -- build/release/sim/riptide_sim --config config/<c>.json \
  --events 10000000 --threads 8 --run-id <c>_1e7
scripts/run_capped.sh -m 5G -- build/release/sim/riptide_sim_check --run-id <c>_1e7 \
  --check-id c1 --bins-u 50 --bins-z 50
```

Le attese si rigenerano con:

```bash
build/release/sim/riptide_sim_predict --config config/<c>.json --events 10000000 \
  --bins-u 50 --bins-z 50
```

## 1. Geometria della sorgente

Accettanza e corda media sono in forma chiusa (D27), quindi non dipendono da nessuna libreria
(a):

| | +X | isotropa |
|---|---|---|
| raggio del disco | 14.142 mm (√(10² + 10²)) | 33.166 mm (mezza diagonale) |
| accettanza attesa = area proiettata / area del disco | 400 / 628.32 = 0.636620 | 1400 / 3455.75 = 0.405122 |
| σ binomiale dell'accettanza a 10⁷ | 1.52·10⁻⁴ | 1.55·10⁻⁴ |
| raggi che colpiscono il blocco, attesi | 6.366·10⁶ | 4.051·10⁶ |
| corda media attesa | 60 mm (una sola corda) | 4V/S = 17.142857 mm (Cauchy) |
| deviazione standard della corda | 0 | 10.49 mm (campione da 10⁸ raggi, (c)) |
| errore sulla corda media a 10⁷ | — | 10.49 / √(4.05·10⁶) = 0.0052 mm |

L'accettanza entra nei test come `geometric_acceptance.pull` e la corda media come
`mean_chord_pull` (§3). Nella direzione fissa la corda ha una sola lunghezza, quindi il pull della
corda è `null` e al suo posto c'è un oracolo esatto (§2).

## 2. Oracoli deterministici (soglie esatte)

Valgono per tutti e tre i run.

| Verifica | Soglia | Etichetta |
|---|---|---|
| `results.counters.must_be_zero.total` (compresi i nuovi `start_not_on_source` e `interaction_off_chord`) | = 0 | (d), solo un bug li rende ≠ 0 |
| `results.deposits.counters.must_be_zero.total` | = 0 | (d), come sopra |
| `deposits.skipped`, `deposits.pstar.skipped` | `null` | (a), schema 4 |
| Σ `deposits.event_classes` | = `deposits.events` = `results.events` = 10⁷ | (a) |
| `geometric_acceptance.rays` | = `results.events` | (a) |
| `geometric_acceptance.hits` ≤ `rays`; `event_classes.no_interaction` ≥ `rays` − `hits` | vero | (a): un raggio che manca il blocco non interagisce |
| solo +X: `chord_mm.min`, `chord_mm.max` | entro 10⁻⁹ mm da 60 mm, e `chord_mm.sd` < 10⁻⁹ mm | (a); pilota: 60 − 7·10⁻¹⁵ mm |
| solo isotropa: `chord_mm.max` | ≤ 2·R_d = 66.33 mm (diagonale del blocco) | (a); pilota: 64.4 / 64.6 mm |
| `deposits.contained_balance_MeV.max` | ≤ 0.01 MeV (slack di D19) | (d); piloti: ≤ 5.6·10⁻⁴ MeV |
| `deposits.contained_balance_MeV.min` | ≥ −10⁻⁹ MeV (protoni senza interazioni adroniche né rinculi nucleari) | (d); piloti: ≥ −7.1·10⁻¹⁵ MeV |
| `deposits.pstar.samples` + `below_table` | = `contained_balance_MeV.samples` | (a), stessi protoni |
| `deposits.pstar.sha256` | `eec195fe…081b` di D20 | (a) |
| `deposits.setup.birks_kB_mm_per_MeV` | 0.155 | (a) |
| `Ep_residuals.max_abs` | < 10⁻⁹·E_n,max | (d); piloti: ≤ 2.1·10⁻¹⁴ MeV |
| `model_mix.observed_non_hp` | = 0: nessuna delle tre sorgenti raggiunge la finestra HP/CHIPS | (a), D25 |
| offset HP − fisico, `rounding_free_relative_{min,max}` | dentro l'intervallo atteso sotto, con margine relativo 10⁻⁹ | (a), D14 |
| `hp_elastic_H1.awr` | 0.999167 | (a) |

Intervallo atteso dell'offset relativo (da `prediction.Ep_physical_offset_relative`):
- mono 14 MeV (+X e isotropa): [5.36581·10⁻⁴, 5.48614·10⁻⁴];
- piatto 10–19.5 MeV: [5.33500·10⁻⁴, 5.50179·10⁻⁴].

I tre piloti ci cadono dentro.

## 3. Test statistici contro `jeff33`

La prima interazione dipende dalla geometria solo attraverso la corda dell'evento. Le attese per
evento usano la corda L_i dell'evento (`checks_version` 10), e la previsione le media sulla
distribuzione delle corde della sorgente (`prediction.chords`) (a). Ipotesi (d): la simulazione
campiona esattamente le attese di `jeff33` lungo la corda di ogni evento.

Piano dei test (Bonferroni, α di famiglia 0.01 per run, dal campo `test_plan`; i pull di
Legendre sono correlati e la correzione resta conservativa):

| Run | χ² | pull | α per test | soglia \|pull\| |
|---|---|---|---|---|
| +X | 3 (u, profondità s = z/L_i, frazioni) | 7 (a₁…a₆, accettanza) | 1.000·10⁻³ | 3.291 |
| isotropa mono | 3 (come sopra) | 8 (a₁…a₆, accettanza, corda media) | 9.09·10⁻⁴ | 3.317 |
| isotropa piatta | 3 (come sopra) | 8 (come sopra) | 9.09·10⁻⁴ | 3.317 |

Nessun run ha il test `model_mix` (D25).

**Esito:** un χ² passa se `p_value` > α per test, un pull se |pull| < soglia. Un run passa se
passano tutti i suoi test e gli oracoli della §2. Con un sistema corretto la probabilità di
almeno un falso allarme sui tre run è circa 3% (b).

### Attese principali (jeff33, 10⁷ eventi) — (b)

Ipotesi:
- prima interazione da sezioni d'urto `jeff33` lungo la corda di ogni evento;
- un raggio che manca il blocco conta come «nessuna interazione»;
- nessun contributo del deuterio (§6);
- quadratura in E_n trascurabile (piatto: 1.6·10⁻⁴σ sulle frazioni, misurata raddoppiando i
  nodi).

| | +X | isotropa mono | isotropa piatta |
|---|---|---|---|
| Σ_tot [1/mm] | 0.0097908 | 0.0097908 | 0.0101486 |
| P(n–p) | 0.104381 | 0.022454 | 0.021910 |
| P(n–C el.) | 0.112058 | 0.024106 | 0.026153 |
| P(n–C non el.) | 0.066366 | 0.014277 | 0.014736 |
| P(altro) | 1.70·10⁻⁵ | 3.65·10⁻⁶ | 5.37·10⁻⁶ |
| P(nessuna), mancati compresi | 0.717178 | 0.939160 | 0.937196 |
| σ relativo P(n–p) | 9.3·10⁻⁴ | 2.1·10⁻³ | 2.1·10⁻³ |
| σ relativo P(n–C non el.) | 1.19·10⁻³ | 2.6·10⁻³ | 2.6·10⁻³ |
| eventi n–p HP | 1.0437·10⁶ | 2.245·10⁵ | 2.191·10⁵ |
| eventi n–p non HP | 0 | 0 | 0 |
| conteggio u minimo per bin | 20 103 | 4 325 | 4 235 |

Lungo +X ogni neutrone che entra attraversa 60 mm, e 1 − e^(−60 Σ) = 0.4442 con 0.6366
dell'accettanza dà la probabilità d'interazione 0.2828 (a, verifica a mano della tabella). Con
la sorgente isotropa l'interazione è più rara: la corda media è 17.1 mm e l'accettanza 0.405.

**Forma di u.** u dipende solo da E_n e dalla cinematica n–p, non dalla geometria (a). Le
probabilità per bin sono quindi quelle della fase 1b alla stessa energia (mono 14 MeV) o quelle
del nuovo spettro 10–19.5 MeV; cambia solo il numero di eventi.

**Profondità.** La profondità s = z/L_i ∈ [0, 1] ha per ogni evento l'attesa di un esponenziale
troncato su [0, L_i]. Con la sorgente isotropa le corde vanno da circa 0 a 65 mm, quindi la
forma attesa è la media di esponenziali quasi piatti (Σ L_i ≲ 0.65). Il χ² non ha un valore
previsto a parte: il test è la sua p-value.

Momenti di Legendre ⟨P_l(μ)⟩ = a_l sugli eventi HP, con l'errore atteso a 10⁷ eventi:

| l | mono 14 (+X): a_l ± σ | mono 14 (isotropa): σ | piatto 10–19.5: a_l ± σ |
|---|---|---|---|
| 1 | −1.4798·10⁻² ± 5.7·10⁻⁴ | 1.22·10⁻³ | −1.4467·10⁻² ± 1.23·10⁻³ |
| 2 | 1.893·10⁻³ ± 4.4·10⁻⁴ | 9.5·10⁻⁴ | 2.431·10⁻³ ± 9.6·10⁻⁴ |
| 3 | −6.61·10⁻⁴ ± 3.7·10⁻⁴ | 8.0·10⁻⁴ | −7.23·10⁻⁴ ± 8.1·10⁻⁴ |
| 4 | 2.52·10⁻⁴ ± 3.3·10⁻⁴ | 7.0·10⁻⁴ | 3.20·10⁻⁴ ± 7.1·10⁻⁴ |
| 5 | −1.0·10⁻⁵ ± 3.0·10⁻⁴ | 6.4·10⁻⁴ | −1.5·10⁻⁵ ± 6.5·10⁻⁴ |
| 6 | 1.4·10⁻⁶ ± 2.7·10⁻⁴ | 5.9·10⁻⁴ | 2.2·10⁻⁶ ± 5.9·10⁻⁴ |

Gli a_l a 14 MeV sono quelli della fase 1b; cambia solo l'errore. Conseguenza (b): a₁ è misurato
a circa 26σ lungo +X e 12σ con la sorgente isotropa; a₂–a₆ sono compatibili con zero entro circa
4σ (+X) e 2.5σ (isotropa).

**Accettanza e corda media.** Pull attesi distribuiti come N(0, 1) (b). Il campione da 10⁸ raggi
della previsione ha dato pull −0.24 sull'accettanza isotropa e 0.93 lungo +X. I piloti da 10⁵ eventi
hanno dato −0.64 (+X), +0.82 e −0.40 (isotropa) sull'accettanza, e −0.43 e +0.12 sulla corda
media (c).

## 4. Confronto con `endfb80` — attese, non soglie

Come nella fase 1b, `endfb80` non è un criterio di accettazione (D13). Differenze attese in
unità dell'errore statistico del run (`shift_from_first_in_sigma`, (b)):

| Run | max \|Δ frazioni\| | n–C el. | n–p | max \|Δa_l\| | non centralità del χ² di u |
|---|---|---|---|---|---|
| +X | 11.05σ (n–C non el.) | +6.8σ | +0.85σ | 0.91σ (a₃) | 1.43 |
| isotropa mono | 5.02σ (n–C non el.) | +3.0σ | +0.34σ | 0.42σ (a₃) | 0.31 |
| isotropa piatta | 4.65σ (n–C non el.) | +3.0σ | +0.28σ | 0.45σ (a₃) | 0.36 |

Previsioni (b):
- il χ² delle frazioni contro `endfb80` fallisce in tutti e tre i run. La non centralità è circa
  170 lungo +X e circa 35 nei due isotropi, su 3 gradi di libertà. Viene dalle valutazioni del
  carbonio e non è un errore;
- su n–p e sugli a_l le due librerie restano entro 1σ: anche qui la simulazione non distingue
  le due revisioni di H-1.

## 5. Confronto con PSTAR (soglie a banda)

La luce e il range per protone dipendono solo da E_p e dal materiale, non dalla geometria. La
geometria decide però quali protoni restano contenuti, e così cambia lo spettro di E_p su cui si
fa la media (a).

**Luce, `pstar.light_MeV.relative_residual`:** stesse bande della fase 2 (b). Il residuo medio
ha un errore statistico di circa 0.0014%/√(N/10⁵) e i piloti ci stanno dentro (c):

| Run | Soglia | Pilota (c) | Fase 2, 10⁷ (c) |
|---|---|---|---|
| +X | [−0.070%, −0.055%] | −0.0618% | mono −0.0620% |
| isotropa mono | [−0.070%, −0.055%] | −0.0618% | mono −0.0620% |
| isotropa piatta | [−0.070%, −0.052%] | −0.0600% | piatto 10–20 −0.0600% |

**Luce per bin di E_p.** Il residuo per bin, (`Evis_mean` − `L_mean`) / `L_mean` in
`pstar_comparison.csv`, è una funzione liscia di E_p (fase 2, (c)). Soglia per i bin con
`Ep_low` ≥ 2 MeV e almeno 100 protoni: entro **0.005% assoluto** dal riferimento. Il
riferimento è committato in `predictions/`: colonne `Ep_low`, `Ep_high`, `reference_residual`,
un bin per riga, compresi tutti i bin del riferimento con almeno 2 protoni. La regola dei 100
protoni si applica al run.
- +X e isotropa mono: `predictions/3_pstar_light_reference_mono_14MeV.csv`. È il residuo dello
  stesso bin in `validation/2_mono_14MeV_1e7/c2/pstar_comparison.csv` sull'SSD. I bordi
  coincidono (50 bin su [0, 14] MeV) (a).
- isotropa piatta: `predictions/3_pstar_light_reference_flat_10_19p5MeV.csv`. È l'interpolazione
  lineare, al centro del bin, dei residui di `validation/2_flat_10_20MeV_1e7/c2/
  pstar_comparison.csv` presi ai centri dei loro bin. I bin sono larghi 0.39 MeV invece di 0.40,
  quindi i bordi non coincidono. Il primo bin, sotto il primo centro, non ha riferimento.

I bin si accoppiano per `Ep_low` uguale entro 10⁻⁶ MeV.

La banda di 0.005% è la stessa della fase 2. Prove fatte prima di fissarla (c):
- i piloti, su tutti i bin da 2 MeV in su, distano dal riferimento al più 0.0010% (+X), 0.0008%
  (isotropa mono) e 0.0019% (isotropa piatta);
- l'interpolazione del piatto applicata ai run mono e di bordo della fase 2 dà al più 0.0014%.
  È la differenza dovuta alla forma dello spettro dentro il bin.

**Range, `pstar.range_mm.relative_residual`:** qui la banda della fase 2 non si eredita. I
protoni che escono sono più energetici della media, quindi la loro esclusione cambia lo spettro
dei protoni confrontati. Il residuo per bin sale verso +8% sotto 3 MeV, quindi la media non è più
quella della fase 2 (−0.106% mono, −0.103% piatto). Le bande sono il valore del pilota ± 4 errori
del pilota (c). A 10⁷ eventi l'errore sarà circa 0.001–0.003%, e la banda è larga perché il
pilota ha pochi protoni:

| Run | Pilota (c) | Soglia |
|---|---|---|
| +X | −0.079% ± 0.011% | [−0.125%, −0.035%] |
| isotropa mono | −0.138% ± 0.024% | [−0.235%, −0.041%] |
| isotropa piatta | −0.017% ± 0.027% | [−0.125%, +0.091%] |

Previsione (d), non soglia: il residuo dei tre run è meno negativo di quello della fase 2, perché
pesa di più la parte sotto 3 MeV. Il pilota isotropo mono va nel verso opposto, ma solo di 1.3σ:
la previsione resta.

## 6. Protoni che escono e altre attese (non soglie)

**Protoni che escono** (`protons_exiting` / `np_events`, protone della prima n–p). Nel blocco la
faccia più vicina dista al più 10 mm, mentre nel cubo della fase 2 distava 50 mm. La frazione
quindi sale di quasi un ordine di grandezza (d, confermata dai piloti (c)):

| | Fase 2, mono 1e7 (c) | +X pilota | isotropa mono pilota | isotropa piatta pilota |
|---|---|---|---|---|
| frazione che esce | 0.47% (13.1 su 2778 per 10⁴) | 3.35% ± 0.17% (372 / 11 095) | 3.98% ± 0.40% (94 / 2 363) | 4.85% ± 0.45% (113 / 2 329) |

Previsioni per i run da 10⁷ eventi, dove l'errore scende a circa 0.02% (+X) e 0.04% (isotropa)
(b):
- ogni frazione cade entro 3 errori del pilota: [2.8%, 3.9%] (+X), [2.8%, 5.2%] (isotropa mono),
  [3.5%, 6.2%] (isotropa piatta);
- il piatto supera il mono isotropo, perché arriva fino a 19.5 MeV e quindi a range più lunghi;
- isotropa mono contro +X: nel pilota la differenza è di 1.5σ, quindi qui non si prevede il
  verso.

Questa è la grandezza nuova della fase 3: l'effetto dei bordi sui protoni (D27). Nei risultati si
riporta anche la distribuzione per E_p, da `proton_exits.csv`.

**Classi dell'evento (D19)**, dai piloti (c; errore binomiale su 10⁵ eventi):

| Classe | +X | isotropa mono | isotropa piatta |
|---|---|---|---|
| `no_interaction` (mancati compresi) | 71.72% ± 0.14 | 93.72% ± 0.08 | 93.67% ± 0.08 |
| `nc_elastic_only` | 10.10% ± 0.10 | 2.37% ± 0.05 | 2.43% ± 0.05 |
| `single_np` | 8.52% ± 0.09 | 1.84% ± 0.04 | 1.85% ± 0.04 |
| `nc_inelastic_no_np` | 7.08% ± 0.08 | 1.55% ± 0.04 | 1.58% ± 0.04 |
| `mixed` | 1.69% ± 0.04 | 0.33% ± 0.02 | 0.29% ± 0.02 |
| `multiple_np` | 0.89% ± 0.03 | 0.19% ± 0.01 | 0.19% ± 0.01 |
| `other` | 0.003% ± 0.002 | 0 | 0 |

Rispetto alla prima interazione, il blocco sottile riduce le interazioni multiple. Nella fase 2
`mixed` + `multiple_np` erano il 18% degli eventi (cubo da 100 mm); qui `single_np` vale circa
4.6 volte (`mixed` + `multiple_np`) lungo +X, contro 0.55 nella fase 2 (b).

Altre attese, scalate dai piloti a 10⁷ eventi (c, errore di Poisson del pilota):

| Grandezza | +X | isotropa mono | isotropa piatta |
|---|---|---|---|
| `np_events` | 1.11·10⁶ | 2.36·10⁵ | 2.33·10⁵ |
| protoni contenuti confrontati con PSTAR | 1.07·10⁶ | 2.26·10⁵ | 2.20·10⁵ |
| `protons_hadronic` | 5 100 ± 700 | 900 ± 300 | 1 500 ± 400 |
| `events_secondary_n_hadronic` | 1.38·10⁵ | 3.0·10⁴ | 2.8·10⁴ |
| `Edep_overflow` | 7 300 ± 850 | 900 ± 300 | < 300 (0 nel pilota) |
| `Evis_overflow` | 0 | 0 | 0 |
| `expected_rare.np_first_without_proton` | circa 120 | circa 26 | circa 25 |

`np_first_without_proton` (n–d con deutone di rinculo, D17) è stimato come n_H2/n_H1 = 1.15·10⁻⁴
per gli eventi n–p, come nella fase 1b (b). Il pilota +X ne ha 1.

**Rinculi nucleari, `protons_nuclear_recoil`** (D22). Il tasso della fase 2 è circa 1.5·10⁻⁵ per
protone contenuto a 14 MeV e 2·10⁻⁵ sul piatto. Con i protoni contenuti della tabella sopra si
attendono circa 16 (+X), 3 (isotropa mono) e 4 (isotropa piatta) casi (b). I piloti ne hanno 0.
Non è una soglia: il conteggio è piccolo e si riporta così com'è.

## 7. Residui noti non modellati

Sono gli stessi della §4 di `1b_predictions.md` e della §6 di `2_predictions.md`, ma con errori
diversi. Ce n'è uno nuovo, il campionamento delle corde.

| Effetto | Grandezza | Rispetto all'errore a 10⁷ | Etichetta |
|---|---|---|---|
| Distribuzione delle corde della previsione da 10⁸ raggi (isotropa) | il secondo campione (seme 20260908) sposta le frazioni di 0.21σ, gli a_l e u di nulla | ≤ 0.21σ sulle frazioni, 0 lungo +X | (c), `chord_sample_shift` |
| C-13 contato due volte nel modello HP (D15) | circa 2.9·10⁻⁴ relativo sul non elastico del carbonio | 0.24σ (+X), 0.11σ (isotropa) | (b) |
| n–d nella verità come n–p (D17) | n_H2/n_H1 = 1.15·10⁻⁴ → circa 120 / 26 / 25 eventi | 0.12σ (+X), 0.05σ (isotropa) sul conteggio n–p | (b) |
| Offset HP della massa di H-1 (D14) | E_p +5.4·10⁻⁴ relativo | modellato nelle verifiche | (a) |
| Stopping nucleare escluso dalla luce, S di Bethe di Geant4 contro PSTAR (D21) | come nella fase 2 | dentro le bande della §5 | (a) / (d) |

## 8. Prestazioni attese (dai piloti, 2026-11-21)

| Grandezza | +X | isotropa mono | isotropa piatta |
|---|---|---|---|
| tempo per 10⁵ eventi a 1 thread (inizializzazione compresa) | 8.4 s | 7.1 s | 7.4 s |
| RSS massimo, sim a 1 thread | circa 91 MB | circa 91 MB | circa 91 MB |
| byte per evento su disco | 110 | 112 | 119 |
| `riptide_sim_check` su 10⁵ eventi | 0.67 s, 65 MB | 0.55 s, 65 MB | 0.58 s, 65 MB |

Molti eventi mancano il blocco e costano solo il volo nel vuoto, per cui un evento richiede meno
tempo che nella fase 2: il pilota v3 della fase 2 faceva 7.8·10³ eventi/s, qui 1.2–1.4·10⁴.
Stime per un run da 10⁷ eventi a 8 thread (b), con la scalatura di 8 thread della fase 2 (4.6·10³
eventi/s per thread contro 7.8·10³ a 1 thread):
- tempo: 2–6 minuti;
- RSS: sotto 0.3 GB;
- spazio su disco: 1.1–1.2 GB per run, circa 3.5 GB in tutto (sull'SSD ne restano 1.7 TB);
- verifica: circa un minuto.

`kBytesPerEvent` (256) resta valido.

## 9. Se una soglia fallisce

Vale la procedura della §6 di `1b_predictions.md` e della §8 di `2_predictions.md`:
- le soglie non si cambiano e non si ripete il run con un altro seme;
- la causa si cerca sul run esistente;
- un oracolo esatto che fallisce è un bug finché non si dimostra il contrario;
- una convenzione non modellata porta a una decisione, a `checks_version` + 1 e a un'altra
  verifica sugli stessi eventi;
- nei risultati si riportano entrambi gli esiti.

In più, per le verifiche nuove della fase 3:
- se fallisce il pull dell'accettanza o della corda media, si controllano prima i raggi uno per
  uno con `rayBoxChord` contro il punto di partenza registrato, poi il campionamento del disco;
- se fallisce il χ² della profondità ma non quello delle frazioni, l'errore è nella corda per
  evento e non nelle sezioni d'urto;
- se il range contro PSTAR esce dalla banda del pilota, nei risultati lo si scrive come
  previsione fallita, studiandolo per bin di E_p.
