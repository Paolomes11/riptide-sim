# Roadmap

Una fase alla volta: si passa alla successiva solo con tutti i test verdi e la fase documentata in
`status.md`. Ogni risultato porta un'etichetta (a)–(d) (vedi `toy_model_findings.md`). Le
previsioni numeriche si scrivono **prima** di ogni run.

A ogni chiusura di fase si aggiorna anche il documento del progetto, `docs/riptide_report.pdf`
(teoria, scelte e risultati in forma leggibile), con la stessa versione e data di `status.md`.

**Fuori scope (deciso):** ottica di lettura (fotoni di scintillazione, lenti, fotocatodo). Il
progetto ricostruisce il neutrone incidente da ciò che interagisce nello scintillatore. Niente
`G4OpticalPhysics`, niente trasporto di fotoni.

## Fase 0 — Infrastruttura ✅ completata (v0.1.0, 2026-11-16)

- Superbuild CMake (C++17, Ninja, preset `debug` / `release` / `ci`), warning severi,
  `-ffast-math` vietato nelle statistiche (`#error`), sanitizer in CI.
- `common/`: unità (identiche a CLHEP), schema degli eventi (bozza v0), radice dei dati
  (`data_paths.hpp`: niente scritture sul disco di sistema, controllo dello spazio).
- Scheletri compilabili di `sim/` (stampa la versione di Geant4) e `analysis/stats`.
- Docs, CI GitHub Actions (senza Geant4/ROOT), wrapper con tetto di memoria.
- Primo commit locale; repository GitHub non ancora creato (solo su richiesta esplicita).

## Fase 1 — Simulazione minima con verità MC

### 1a — Simulazione funzionante ✅ completata (v0.2.0, 2026-11-17)

- Mondo nel vuoto + blocco di BC-408 (`G4_PLASTIC_SC_VINYLTOLUENE`, dimensioni dal config),
  `QGSP_BIC_HP`, run MT.
- Fascio di neutroni monoenergetico o piatto su [E_min, E_max], direzione fissa.
- CLI (`--config`, `--data-dir`, `--events`, `--threads`, `--seed`, `--run-id`), configurazione
  JSON versionata in `config/`, `run.json` scritto accanto ai dati (config, seme, git, versioni).
- Semi per evento da (event_id, seme del run): output indipendente dal numero di thread (D11).
- Ntupla di verità secondo `event_schema.hpp` v1, un file per thread.
- Test d'integrazione `sim_integration`: run a 1 e a 2 thread identici colonna per colonna.

### 1b — Verifiche fisiche e prestazioni ✅ completata (v0.3.0, 2026-11-19)

- Verifiche **non circolari**:
  - spettro E_p del primo scattering n–p piatto su [0, E_n] (a);
  - E_p(θ_p) contro la formula relativistica (a);
  - distribuzione angolare in CM contro JEFF-3.3 (implementazione di Geant4) ed ENDF/B-VIII.0;
    per H-1 le due librerie sono la stessa analisi LANL, quindi il confronto n–p indipendente
    resta da fare (D13);
  - frazione di prime interazioni n–p / n–C contro le sezioni d'urto;
  - profondità del primo vertice contro exp(−Σ_tot z);
  - comportamento a 19.5–20 MeV (bordo di validità di `G4NeutronHP`);
  - convenzioni del modello HP verificate sui dati che Geant4 usa davvero: massa del bersaglio di
    H-1 (D14), cattura esclusa dall'inelastico del carbonio e composizione isotopica del
    carbonio (D15, oracolo G4NDL ↔ JEFF-3.3); sopra 19.5 MeV elastico HP e CHIPS sovrapposti
    (D16); break-up n–d su H-2 nelle attese, da G4NDL (D17).
- Prestazioni: eventi/s per thread, byte/evento misurato (`kBytesPerEvent` per
  `estimateRunBytes`).
- Risultati: `validation/1b_results.md` (tre run da 10⁷, tutte le verifiche passate).

## Fase 2 — Deposito di energia e quenching di Birks ✅ completata (v0.4.0, 2026-11-21)

- Energia depositata e luce visibile con `G4EmSaturation` (kB da letteratura citata, mai
  inventato): kB = 0.155 mm/MeV, intervallo sistematico [0.145, 0.165] (D18).
- Luce dei protoni continua, L(K_pre) − L(K_pre − E_dep) dallo S elettronico di Geant4,
  indipendente dalla funzione di passo (D21).
- Traccia del protone: inizio, fine, lunghezza, uscita dal volume, interazioni adroniche e
  rinculi nucleari da scattering Coulombiano (D22).
- Classificazione dell'evento: singolo n–p, n–p multiplo, n–C elastico/inelastico, misto (D19).
- Verifiche dei depositi (schema 3, `checks_version` 8, poi 9): oracoli di bilancio, confronto di range
  e luce con PSTAR per bin di E_p (D20).
- Lettura delle righe in ordine di `event_id`: le verifiche non dipendono dalla divisione degli
  eventi fra i thread (D23, `checks_version` 9).
- Run di validazione: tre da 10⁷ (mono 14 MeV, piatto 10–20 MeV, bordo 19.5–20 MeV) e tre di
  sensibilità a kB, con previsioni in `validation/2_predictions.md` v2.
- Risultati: `validation/2_results.md` (tutte le soglie passate nelle verifiche `c2`).

## Fase 3 — Geometria reale ✅ completata (v0.5.0, 2026-11-22)

- Solo il blocco scintillatore, senza componenti ottici: 60 × 20 × 20 mm di BC-408, lato lungo
  su X, l'unico volume 3D del progetto di ottimizzazione (D27).
- Sorgente a campo lontano (D27): direzione fissa (Caso A) o isotropa (Casi B–C), punto di
  partenza uniforme su un disco ⟂ Ω_n. Il fascio a matita resta per la regressione.
- Config v3, schema 4 (`truth_n_start_{x,y,z}`), verifiche per evento con la corda L_i,
  accettanza geometrica e corda media di Cauchy (`checks_version` 10).
- Run di validazione: tre da 10⁷ eventi (mono 14 MeV lungo +X, mono 14 MeV isotropo, piatto
  10–19.5 MeV isotropo), con previsioni in `validation/3_predictions.md`.
- Risultati: `validation/3_results.md` (tutte le soglie passate in `c1`), con l'effetto dei
  bordi sui protoni che escono.

## Fase 4 — I/O dell'analisi ✅ completata (v0.6.0, 2026-11-22)

- `analysis/io`: `RDataFrame` → buffer struct-of-arrays (nessun ROOT a valle).
- Selezione dei candidati a singolo scattering n–p.
- Risposta del rivelatore (σ_Ep, σ_θ) applicata **in analisi**, con seme registrato: la stessa
  simulazione serve a più risoluzioni.
- `truth_Ep` convertito all'E_p fisico con `physicalProtonEnergyFromHp` prima di σ_Ep; la massa
  HP del bersaglio viene dai dati della run, mai scritta a mano (D24, D28).
- Risultati: `validation/4_results.md` (tutte le soglie passate in `a2`, dopo la correzione degli
  eventi senza traccia del protone), con D28–D31.

## Fase 5 — Port del motore statistico (Caso A) ✅ completata (v0.7.0, 2026-11-22)

- `analysis/stats`: griglie, prior, logsumexp, modello del Caso A (D33), combinazione,
  validazione, generatore di pseudo-esperienze (D35).
- Oracoli R1/R2/R3 e test di limite come test Catch2.
- Sorgente con E_n ≤ 19.5 MeV nei run delle fasi 5–6 (D25). Dominio del prior: parametro, default
  [10, 20] MeV (D32).
- Taglio `has_proton` per gli eventi `single_np` senza protone tracciato (D34).
- Calibrazione C++ (`riptide_stats_calibrate`) confrontata con i numeri del toy.
- Target di prestazione più stretti del toy (log-verosimiglianza N = 1000 ≤ 0.10 s, contro
  0.71 s), griglie contigue, parallelismo sugli eventi solo dopo aver profilato (D35).

## Fase 6 — Caso A su dati Geant4 ✅ completata (v0.8.0, 2026-11-25)

- Checklist completa: bias → risoluzione → pull → coverage → contrazione 1/√N → robustezza al
  prior.
- **Confronto toy vs Geant4 sullo stesso motore**: misura della circolarità del toy.

## Fase 7 — Caso B 🔄 in corso (design: D41–D45, `validation/7_predictions.md`)

Direzione Ω_n su griglia di Fibonacci con raffinamento locale; campo lontano.

## Fase 8 — Caso C

Due stadi (Ω_n come B, poi griglia 2D su μ_E, σ_E); mai griglia 4D. Verifica della troncatura al
dominio.

## Fase 9 — Realismo

σ_np(E) nella verosimiglianza (spettro rivelato ≠ incidente), soglia su E_p, eventi n–C nel
modello, sorgenti multiple e fondo (mistura con frazione f), spettri non gaussiani. Sistematiche:
bordo 19.5–20 MeV del modello HP (D25), eccesso di `ionIoni` sui rinculi C12 (D26) e doppio
conteggio del C-13 (D15).

## Fase 10 — Produzione batch

Run su cluster (CNAF) se i volumi lo richiedono; unione dei file per thread solo se serve.
