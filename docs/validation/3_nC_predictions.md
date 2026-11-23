# Previsioni — pull del n–C non elastico lungo +X (run gemello)

Versione 1 (2026-11-23), scritta **prima** del run. Riguarda il punto aperto della fase 3
(`3_results.md`, §2): nel run `3_ff_mono_14MeV_x_1e7` la frazione di prime interazioni n–C non
elastiche è a +3.73σ da `jeff33`. Le soglie e le regole di lettura qui sotto non si toccano dopo
il run.

Etichette: **(a)** fatto stabilito; **(b)** stima con ipotesi esplicite; **(c)** verificato con
Monte Carlo; **(d)** ipotesi da mettere alla prova.

## 1. Cosa si sa già, dagli stessi eventi

**Frazioni condizionate all'interazione (c).** Fra le prime interazioni n–p, n–C elastiche e
n–C non elastiche, le frazioni dipendono solo dai rapporti delle sezioni d'urto a E_n, non dalla
geometria. Le attese `jeff33` sono quindi le stesse per tutti i run a 14 MeV: 0.36909 / 0.39624 /
0.23467.

| Campione a 14 MeV (indipendenti) | Interazioni | n–C non el., frazione | Pull |
|---|---|---|---|
| `3_ff_mono_14MeV_x_1e7` (+X) | 2 831 300 | 0.23552 | **+3.37** |
| `3_ff_mono_14MeV_isotropic_1e7` | 608 232 | 0.23411 | −1.03 |
| `2_mono_14MeV_1e7` (fascio a matita, cubo da 100 mm) | 6 241 491 | 0.23433 | −2.03 |

`1b_mono_14MeV_1e7` ha gli stessi eventi di `2_mono_14MeV_1e7` e non conta come campione
indipendente. L'eccesso c'è solo lungo +X.

**Diagnostico per profondità e per facce (c).** È una macro RDataFrame fuori dal repo, in sola
lettura su `sim/3_ff_mono_14MeV_x_1e7`. Calcola la frazione condizionata del n–C non elastico
in 6 bin di profondità lungo il fascio (0–60 mm) e in 5 bin di distanza dalla faccia più vicina
(0–10 mm). Il riferimento è la frazione media del run.
- **profondità:** χ² = 2.00 su 5 gradi di libertà, nessun andamento;
- **facce:** χ² = 8.22 su 4 gradi di libertà (p ≈ 0.08), con pull +2.02 / −1.55 / +0.04 /
  −1.20 / −0.54 dalla faccia verso il centro, senza un andamento monotono.

L'eccesso è distribuito uniformemente nel blocco. Le due cause indicate in `3_results.md` §2
(profondità dei vertici, vicinanza alle facce) non trovano conferma **(c)**. Resta l'ipotesi
della fluttuazione statistica **(d)**: un |pull| ≥ 3.37 ha p ≈ 7.5·10⁻⁴ a due code, e su circa
10 pull per categoria guardati nei run da 10⁷ la probabilità che uno arrivi così lontano è circa
10⁻².

## 2. Run

| Run (`run_id`) | Config | Eventi | Thread | Seme |
|---|---|---|---|---|
| `3_ff_mono_14MeV_x_1e7_s20260908` | `config/3_ff_mono_14MeV_x.json` | 10⁷ | 8 | `20260908` |

- Build `release` all'HEAD di questo branch. Rispetto a `ded6c88`, la build della fase 3, in
  `sim/` sono cambiati solo `run.json` v2 (D28) e il reader: geometria, fisica e generatore sono
  identici **(a)**.
- `20260908` è il seme diagnostico già previsto nella fase 5 (`5_results.md`), dove poi non è
  servito: con un altro seme gli eventi sono indipendenti da quelli del run originale.
- Comando, dopo `scripts/data_root.sh check`:
  `scripts/run_capped.sh -m 5G -- build/release/sim/riptide_sim --config
  config/3_ff_mono_14MeV_x.json --events 10000000 --threads 8 --seed 20260908 --run-id
  3_ff_mono_14MeV_x_1e7_s20260908`;
- poi `riptide_sim_check --check-id c1 --bins-u 50 --bins-z 50`, come nella fase 3.
- **Stima (b), dal run originale:** circa 2.8 minuti e meno di 0.2 GB di RSS per la
  simulazione, circa 50 s per le verifiche, 1.1 GB su disco.

## 3. Previsioni e regole di lettura

- **Verifiche `c1` (b):** passano tutte le soglie della fase 3 (`3_predictions.md` v2), con gli
  stessi valori attesi entro la statistica:
  - circa 6.37·10⁶ raggi sul blocco;
  - circa 2.83·10⁶ interazioni;
  - n–C non elastiche attese 663 889 ± 787, dalla stessa attesa del run originale.
- **Test principale (d → c).** È il pull del n–C non elastico contro `jeff33` nel run gemello.
  - |pull| < 3.0: fluttuazione. Il punto aperto si chiude e si riporta il pull combinato dei due
    run, (3.73 + x)/√2.
  - pull ≥ +3.0, con lo stesso segno: effetto sistematico della configurazione +X. Prima di
    toccare il codice si apre una decisione (D36) e si cerca nel generatore a campo lontano e nel
    conteggio della prima interazione.
  - pull ≤ −3.0: due anomalie opposte, da studiare con un terzo seme.
- Si riportano anche la frazione condizionata (§1) e il diagnostico per profondità e per facce
  sul run gemello. Sono attese, non soglie.
