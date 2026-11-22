# Previsioni della fase 4 — scritte prima dei run

Versione 1 (2026-11-22). Va committata prima di qualsiasi run di `riptide_analysis_prepare` sui
run della fase 3. I risultati (`4_results.md`) citeranno l'hash di quel commit. Dopo i run le
soglie non si toccano.

La fase 4 non simula niente: legge in sola lettura i tre run da 10⁷ eventi della fase 3. Poi
applica la selezione (D30), la conversione HP (D28) e la risposta del rivelatore (D31). Le chiavi
di `analysis.json` e il formato di `cutflow.csv` sono descritti in
`docs/architecture/data_layout.md`.

Le fonti dei numeri sono due:
- `validation/<run>/c1/checks.json` dei tre run, prodotti nella fase 3 (`3_results.md`);
- le attese di `riptide_sim_predict` in `predictions/3_ff_*_1e7.json`, per la banda analitica
  dell'offset HP.

Etichette: **(a)** fatto stabilito (sorgenti, dati, costruzione del codice); **(b)** stima con
ipotesi esplicite; **(c)** verificato con Monte Carlo; **(d)** ipotesi da mettere alla prova.

## Run previsti

Si usa la build `release` del commit di questo file. `RIPTIDE_BUILD_ROOT_IO` è ON per default; la
build `release` va riconfigurata, perché è stata configurata prima che l'opzione esistesse. I run della
fase 3 hanno `run_record_version` 1, quindi la massa HP viene da `--hp-from-check c1` (D28).
Ogni comando gira sotto `scripts/run_capped.sh -m 5G`, uno alla volta, dopo
`scripts/data_root.sh check`. Parametri del toy: σ_Ep = 0.10 MeV, σ_θ = 0.08 rad.
Il seme è quello del progetto, `20260907`, e la selezione è quella di default per la fase 5,
`single_np,contained`.

```bash
scripts/run_capped.sh -m 5G -- build/release/analysis/apps/riptide_analysis_prepare \
  --run-id <run> --analysis-id a1 --seed 20260907 --sigma-Ep 0.10 --sigma-theta 0.08 \
  --cuts single_np,contained --hp-from-check c1
```

| Run (`run_id`) | `analysis_id` | Ruolo |
|---|---|---|
| `3_ff_mono_14MeV_x_1e7` | `p1e5` | `--max-events 100000`, misura dell'RSS prima dei run completi |
| `3_ff_mono_14MeV_x_1e7` | `a1` | campione del Caso A |
| `3_ff_mono_14MeV_isotropic_1e7` | `a1` | direzioni isotrope, mono |
| `3_ff_flat_10_19p5MeV_isotropic_1e7` | `a1` | sorgente di produzione |
| `3_ff_mono_14MeV_x_1e7` | `repro_a`, `repro_b` | `--max-events 1000000`, stesso seme, riproducibilità |

**Memoria (b).** Il lettore tiene le 13 colonne di un run: 3 intere da 4 byte e 10 doppie da
8 byte. Per 10⁷ eventi sono ~0.9 GB, più una copia temporanea per colonna durante il riordino
per `event_id`, più ROOT. La stima dà un RSS < 2.5 GB. Il pilota da 10⁵ eventi misura l'RSS. Se
l'estrapolazione lineare a 10⁷ supera 4 GB, i run completi non partono e si torna al piano.

## 1. Selezione

**Conteggi di `single_np` (a).** La classe è un campo dell'evento, e `riptide_sim_check` la
conta sugli stessi file (`deposits.event_classes.single_np`). Il taglio `single_np` del cut-flow
deve quindi coincidere con quel conteggio **esattamente**. Anche gli eventi letti devono valere
esattamente 10⁷.

| Run | `all` | `single_np` |
|---|---|---|
| +X, 14 MeV | 10 000 000 | 852 271 |
| isotropa, 14 MeV | 10 000 000 | 186 115 |
| isotropa, 10–19.5 MeV | 10 000 000 | 181 232 |

**Protoni contenuti (d).** `riptide_sim_check` non conta i protoni che escono dagli eventi
`single_np`, ma conta quelli di tutti gli n–p validi. Per questi ultimi, la frazione che esce è
`protons_exiting` / (`protons_exiting` + `contained_balance_MeV.samples`):

| Run | protoni che escono / tutti i protoni n–p |
|---|---|
| +X, 14 MeV | 3.19% |
| isotropa, 14 MeV | 4.39% |
| isotropa, 10–19.5 MeV | 4.78% |

In un evento `single_np` il protone nasce nello stesso volume e con la stessa cinematica.
L'ipotesi è che la sua frazione di uscita sia la stessa entro un punto percentuale. Il rapporto
`contained` / `single_np` atteso è quindi:

| Run | banda |
|---|---|
| +X, 14 MeV | [0.958, 0.978] |
| isotropa, 14 MeV | [0.946, 0.966] |
| isotropa, 10–19.5 MeV | [0.942, 0.962] |

Un rapporto fuori banda non blocca la fase. Va però spiegato prima della fase 5, perché cambia il
campione del Caso A.

## 2. Conversione HP

Le soglie seguenti valgono sugli eventi selezionati, da `summaries.hp_conversion`.

- **Eventi senza offset (a).** `events_without_offset` è il numero di eventi selezionati con
  `truth_Ep_phys` = 0, cioè con `truth_valid_np` = 0. Non può superare il contatore
  `np_first_without_proton` dello stesso `checks.json`: 108 / 21 / 23.
- **Banda dell'offset (a).**
  - Un evento `single_np` ha un solo urto, che è n–p. Quindi E_n è l'energia della sorgente e
    θ_p è misurato dalla direzione del fascio.
  - `truth_Ep_phys` è l'inversa esatta della cinematica HP. Perciò l'offset relativo deve cadere
    nella banda analitica di `riptide_sim_predict` (`Ep_physical_offset_relative`).
  - La tolleranza è 10⁻⁹ in assoluto, molto sopra l'arrotondamento in doppia precisione.

  | Run | min | max |
  |---|---|---|
  | +X e isotropa, 14 MeV | 5.365814·10⁻⁴ | 5.486139·10⁻⁴ |
  | isotropa, 10–19.5 MeV | 5.335002·10⁻⁴ | 5.501793·10⁻⁴ |
- **Cinematica fisica (a/c).** `lab_angle_relative_max` < 6·10⁻⁶ in ogni run (D31). Sul run
  d'integrazione da 1000 eventi si è misurato 4.0·10⁻⁶.

## 3. Risposta del rivelatore

Le soglie seguenti si leggono da `summaries.resolution`.

- **Pull (b).** Con i semi per evento e un campionatore gaussiano corretto, ogni pull è N(0, 1).
  Le soglie sono:
  - |`pull_mean`| < 3.3 e |`pull_variance`| < 3.3, per `Ep_MeV` e per `theta_p_rad`;
  - |`p_dir.pull_mean`| < 3.3.

  Sono 15 confronti. Con 3.3σ la probabilità di almeno un falso fallimento è ~1.5%.
- **Larghezze (b).** `std` deve coincidere con σ entro 4·σ/√(2N), con N preso dal taglio finale:
  `Ep_MeV.std` con 0.10 MeV, `theta_p_rad.std` con 0.08 rad.
- **Versori (a).** `p_dir.norm_deviation_max` < 10⁻¹².
- **Estremi (b).** Fra ~8·10⁵ gaussiane il massimo atteso di |x|/σ è circa 4.9. La soglia è
  `max_abs` < 6.5 σ.

## 4. Riproducibilità

- **Stesso seme (a).** `repro_a` e `repro_b` sono costruiti uguali, quindi i loro `analysis.json`
  devono essere identici byte per byte, a meno di `analysis_id`. Anche i `cutflow.csv` devono
  essere identici.
- **Sottoinsieme (a).** Il motore è per evento (D31), quindi non dipende da quanti eventi si
  leggono. I primi 10⁶ eventi di `repro_a` danno un cut-flow con tutti i conteggi ≤ quelli di `a1`
  dello stesso run.

## Cosa non si verifica qui

Ci sono tre cose che questa fase non verifica:
- il Caso A stesso (fase 5);
- l'effetto degli eventi a più urti, dei protoni che escono e della soglia `evis_min`. Lo si
  misurerà come sistematica (D30);
- il pull aperto di +3.73σ della fase 3 (n–C non elastico, +X). Non entra nella selezione, che
  usa la classe vera; resta da indagare prima della fase 6.
