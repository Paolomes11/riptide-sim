# Risultati della fase 4 — I/O dell'analisi

Data: 2026-11-22. Previsioni e soglie: `docs/validation/4_predictions.md` versione 1, committata
in `d19fd67` **prima** dei run. Dopo i run le soglie non sono state toccate.

**Esito.** Il primo giro (`a1`, build di `d19fd67`) ha fallito una soglia:
`p_dir.pull_mean` era `null` nei tre run. La causa è un difetto del codice, non della fisica: gli
eventi selezionati senza protone producevano una direzione ricostruita NaN (§3). Dopo la
correzione (`6f8817e`) il secondo giro (`a2`) passa **tutte** le soglie di §1–§4. `a2` coincide
con `a1` in tutto tranne il blocco `p_dir`, che `a1` calcolava male.

Etichette: **(a)** fatto stabilito; **(b)** stima con ipotesi esplicite; **(c)** verificato con
Monte Carlo; **(d)** ipotesi da mettere alla prova.

## Run

Build `release` con `RIPTIDE_BUILD_ROOT_IO` ON. `a1` è stato eseguito a `d19fd67`, `a2` e
`repro_c`/`repro_d` a `6f8817e` (`git describe` v0.5.0-10-g6f8817e). Argomenti comuni:
- seme 20260907;
- σ_Ep = 0.10 MeV, σ_θ = 0.08 rad;
- `--cuts single_np,contained`;
- `--hp-from-check c1`.

Ogni comando è passato da `scripts/run_capped.sh -m 5G`, uno alla volta, dopo
`scripts/data_root.sh check`. Gli output stanno sull'SSD, in `analysis/<run>/<analysis_id>/`.

| Run | `analysis_id` | Eventi letti | Selezionati | RSS massimo | Tempo |
|---|---|---|---|---|---|
| `3_ff_mono_14MeV_x_1e7` | `p1e5` | 10⁵ | 8 256 | 313 MB | 2.9 s |
| `3_ff_mono_14MeV_x_1e7` | `repro_a`, `repro_b` | 10⁶ | 82 919 | 413 MB | 2.1 s |
| `3_ff_mono_14MeV_x_1e7` | `a1` | 10⁷ | 826 859 | 1.37 GB | 8.6 s |
| `3_ff_mono_14MeV_isotropic_1e7` | `a1` | 10⁷ | 178 067 | 1.36 GB | 9.2 s |
| `3_ff_flat_10_19p5MeV_isotropic_1e7` | `a1` | 10⁷ | 172 774 | 1.36 GB | 9.6 s |
| `3_ff_mono_14MeV_x_1e7` | `a2` | 10⁷ | 826 859 | 1.36 GB | 9.2 s |
| `3_ff_mono_14MeV_isotropic_1e7` | `a2` | 10⁷ | 178 067 | 1.36 GB | 8.9 s |
| `3_ff_flat_10_19p5MeV_isotropic_1e7` | `a2` | 10⁷ | 172 774 | 1.36 GB | 9.8 s |
| `3_ff_mono_14MeV_x_1e7` | `repro_c`, `repro_d` | 10⁶ | 82 919 | 415 MB | 2.0–2.2 s |

L'RSS e il tempo trascorso vengono da `/usr/bin/time -v` (a). Per `repro_a` e `repro_b` è
riportato un solo valore.

**Memoria (b → a).** La previsione era un RSS < 2.5 GB. Il massimo misurato è 1.37 GB, quindi
passa. Fra 10⁵ e 10⁷ eventi l'RSS cresce di circa 110 byte per evento, sopra una base di circa
0.3 GB. L'estrapolazione del pilota restava ben sotto il limite di 4 GB.

## 1. Selezione

| Run | `all` | `single_np` | atteso | `contained` | `contained` / `single_np` | banda (d) |
|---|---|---|---|---|---|---|
| +X, 14 MeV | 10 000 000 | 852 271 | 852 271 | 826 859 | 0.97018 | [0.958, 0.978] |
| isotropa, 14 MeV | 10 000 000 | 186 115 | 186 115 | 178 067 | 0.95676 | [0.946, 0.966] |
| isotropa, 10–19.5 MeV | 10 000 000 | 181 232 | 181 232 | 172 774 | 0.95333 | [0.942, 0.962] |

- **Conteggi (a):** gli eventi letti sono esattamente 10⁷ e `single_np` coincide esattamente con
  `deposits.event_classes.single_np` di `c1` ✅.
- **Protoni contenuti (d):** nei tre run il rapporto cade nella banda ✅. La frazione di uscita
  `1 − contained/single_np` vale 2.98% / 4.32% / 4.67%. È circa 0.1–0.2 punti sotto quella di
  tutti gli n–p validi (3.19% / 4.39% / 4.78%), cioè dentro il punto percentuale dell'ipotesi.

## 2. Conversione HP

Questi valori sono uguali in `a1` e in `a2`.

| Run | `events_without_offset` | tetto (a) | offset min | offset max | banda (a) | `lab_angle_relative_max` |
|---|---|---|---|---|---|---|
| +X, 14 MeV | 89 | 108 | 5.365815·10⁻⁴ | 5.445855·10⁻⁴ | [5.365814, 5.486139]·10⁻⁴ | 4.03·10⁻⁶ |
| isotropa, 14 MeV | 17 | 21 | 5.365815·10⁻⁴ | 5.445854·10⁻⁴ | [5.365814, 5.486139]·10⁻⁴ | 4.03·10⁻⁶ |
| isotropa, 10–19.5 MeV | 19 | 23 | 5.335497·10⁻⁴ | 5.445860·10⁻⁴ | [5.335002, 5.501793]·10⁻⁴ | 5.58·10⁻⁶ |

- **Eventi senza offset (a):** sono sotto il tetto `np_first_without_proton` ✅. Sono eventi
  `single_np` con `truth_valid_np` = 0, quindi senza protone tracciato.
- **Banda dell'offset (a):** min e max cadono nella banda analitica, con una tolleranza di 10⁻⁹
  ✅. Nei run a 14 MeV il minimo coincide con il bordo inferiore alle cifre riportate; il
  massimo resta sotto il bordo superiore.
- **Cinematica fisica (a/c):** `lab_angle_relative_max` < 6·10⁻⁶ in ogni run ✅.

## 3. Risposta del rivelatore

### Il difetto di `a1` (a)

In `a1` `p_dir.pull_mean` era `null` nei tre run, quindi la soglia |`p_dir.pull_mean`| < 3.3 è
**fallita** ❌. La causa:
- per gli eventi con `truth_valid_np` = 0 (gli stessi 89 / 17 / 19 della §2) `truth_p_dir`
  vale (0, 0, 0);
- `smearDirection` normalizzava quel vettore nullo, e `reco_p_dir` diventava NaN;
- nella media il NaN si propagava al pull, che nlohmann scrive come `null`;
- `std::max(x, NaN)` restituisce `x`, quindi `angle_max` e `norm_deviation_max` nascondevano il
  NaN. Anche il valore di `norm_deviation_max` di `a1` non era una verifica vera.

La correzione (`6f8817e`, addendum a D31 in `decisions.md`):
- una traccia nulla significa «nessuna traccia»: `reco_p_dir` = (0, 0, 0) e nessuna estrazione
  per la direzione;
- una traccia che non è né nulla né un versore è un errore;
- i riassunti saltano gli eventi senza traccia e li contano (`events_without_track`, record
  versione 2). Un residuo non finito fa fallire il run invece di sparire in un massimo;
- la selezione (D30) e i conteggi non cambiano. Il seme per evento fa sì che l'Ep e il θ_p degli
  eventi senza traccia, e tutto quello degli altri eventi, restino identici bit per bit.

Il confronto fra `a1` e `a2` con `jq` (a) mostra che i due `analysis.json` sono uguali dopo aver
tolto `analysis_id`, `analysis_record_version`, `git_describe` e `p_dir`. Anche i `cutflow.csv`
sono identici. Le soglie non sono state cambiate.

### Risultati di `a2`

| Grandezza | +X, 14 MeV | isotropa, 14 MeV | isotropa, 10–19.5 MeV | soglia |
|---|---|---|---|---|
| `Ep_MeV.pull_mean` | 0.839 | 0.085 | 0.056 | < 3.3 |
| `Ep_MeV.pull_variance` | −1.016 | −0.374 | 0.555 | < 3.3 |
| `Ep_MeV.std` (MeV) | 0.099921 | 0.099937 | 0.100094 | 0.10 ± 4σ/√(2N) |
| `Ep_MeV.max_abs` / σ | 4.81 | 4.73 | 4.73 | < 6.5 |
| `theta_p_rad.pull_mean` | −1.358 | −0.039 | −0.564 | < 3.3 |
| `theta_p_rad.pull_variance` | −0.318 | −0.302 | 1.091 | < 3.3 |
| `theta_p_rad.std` (rad) | 0.079980 | 0.079960 | 0.080148 | 0.08 ± 4σ/√(2N) |
| `theta_p_rad.max_abs` / σ | 5.12 | 4.86 | 4.58 | < 6.5 |
| `p_dir.pull_mean` | −0.355 | −1.918 | 0.550 | < 3.3 |
| `p_dir.scaled_angle_squared_mean` | 0.99961 | 0.99546 | 1.00132 | — |
| `p_dir.events_without_track` | 89 | 17 | 19 | = §2 |
| `p_dir.norm_deviation_max` | 4.4·10⁻¹⁶ | 6.7·10⁻¹⁶ | 6.7·10⁻¹⁶ | < 10⁻¹² |

- **Pull (b):** i 15 pull hanno |pull| < 3.3 ✅. Il più grande è −1.92, sulla direzione nel run
  isotropo a 14 MeV.
- **Larghezze (b):** tutte entro 4σ/√(2N) ✅, con N = 826 859 / 178 067 / 172 774. La deviazione
  più grande, relativa al margine, è quella di θ_p nel piatto: 1.5·10⁻⁴ rad su 5.4·10⁻⁴.
- **Versori (a):** `norm_deviation_max` < 10⁻¹² ✅, ora su tutti gli eventi con traccia.
- **Estremi (b):** `max_abs` < 6.5σ ✅. Il massimo è 5.12σ su θ_p lungo +X, contro un'attesa
  di ~4.9 per ~8·10⁵ gaussiane.
- `events_without_track` coincide con `events_without_offset` della §2 (a), come deve: sono gli
  stessi eventi con `truth_valid_np` = 0.

## 4. Riproducibilità

- **Stesso seme (a):** `repro_a` ≡ `repro_b` (build di `d19fd67`) e `repro_c` ≡ `repro_d` (build
  di `6f8817e`). In ogni coppia i due `analysis.json` sono identici a meno di `analysis_id`, e i
  due `cutflow.csv` sono identici byte per byte ✅.
- **Sottoinsieme (a):** sui primi 10⁶ eventi il cut-flow è 1 000 000 / 85 448 / 82 919, con ogni
  conteggio ≤ quello di `a1` (10⁷ / 852 271 / 826 859) ✅.

## Punti aperti

- **Eventi con `truth_Ep_phys` = 0 nel campione del Caso A (a):** sono 89 / 17 / 19 eventi
  `single_np` contenuti senza protone tracciato. La selezione li tiene (D30) e il loro `reco_Ep` è
  solo rumore attorno a 0. Va deciso nella fase 5 come entrano nella verosimiglianza del Caso A:
  esclusi da un taglio sulla traccia, o tenuti come parte della risposta.
- **Warning della build `release` (a):** con `-O2` GCC dà `-Wnull-dereference` e `-Warray-bounds`
  dentro header esterni (libstdc++ `streambuf`, nlohmann/json, Geant4 `G4NavigationHistory.icc`).
  Erano già noti dalla fase 2 e non vengono dal codice del progetto. Le build `debug` e `ci`
  sono senza warning.
- Restano aperti come prima: il pull di +3.73σ della fase 3 (n–C non elastico, +X), da indagare
  prima della fase 6, e la misura degli eventi a più urti, dei protoni che escono e della soglia
  `evis_min` come sistematica (D30).
