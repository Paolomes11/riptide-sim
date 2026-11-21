# Roadmap

Una fase alla volta: si passa alla successiva solo con tutti i test verdi e la fase documentata in
`status.md`. Ogni risultato porta un'etichetta (a)–(d) (vedi `toy_model_findings.md`). Le
previsioni numeriche si scrivono **prima** di ogni run.

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

## Fase 2 — Deposito di energia e quenching di Birks 🔄 in corso (branch `feat/sim-deposits`)

- Energia depositata e luce visibile con `G4EmSaturation` (kB da letteratura citata, mai
  inventato): kB = 0.155 mm/MeV, intervallo sistematico [0.145, 0.165] (D18).
- Luce dei protoni continua, L(K_pre) − L(K_pre − E_dep) dallo S elettronico di Geant4,
  indipendente dalla funzione di passo (D21).
- Traccia del protone: inizio, fine, lunghezza, uscita dal volume, interazioni adroniche e
  rinculi nucleari da scattering Coulombiano (D22).
- Classificazione dell'evento: singolo n–p, n–p multiplo, n–C elastico/inelastico, misto (D19).
- Verifiche dei depositi (schema 3, `checks_version` 8): oracoli di bilancio, confronto di range
  e luce con PSTAR per bin di E_p (D20).
- Run di validazione: tre da 10⁷ (mono 14 MeV, piatto 10–20 MeV, bordo 19.5–20 MeV) e tre di
  sensibilità a kB, con previsioni in `validation/2_predictions.md` v2.

## Fase 3 — Geometria reale

- Solo il blocco scintillatore (dimensioni e posizione dal progetto di ottimizzazione, nessun
  componente ottico), sorgente a distanza finita.
- **Punto aperto:** lì lo scintillatore non è un volume GDML e i parametri differiscono fra config
  e test ⇒ prima di questa fase servono le dimensioni definitive.

## Fase 4 — I/O dell'analisi

- `analysis/io`: `RDataFrame` → buffer struct-of-arrays (nessun ROOT a valle).
- Selezione dei candidati a singolo scattering n–p.
- Risposta del rivelatore (σ_Ep, σ_θ) applicata **in analisi**, con seme registrato: la stessa
  simulazione serve a più risoluzioni.

## Fase 5 — Port del motore statistico (Caso A)

- `analysis/stats`: kinematics, forward_model, posterior_A, combine, validate.
- Oracoli R1/R2/R3 e test di limite come test Catch2.
- Target di prestazione più stretti del toy (log-verosimiglianza N = 1000 ≪ 0.71 s), griglie
  contigue, parallelismo sugli eventi (TBB) solo dopo aver profilato.

## Fase 6 — Caso A su dati Geant4

- Checklist completa: bias → risoluzione → pull → coverage → contrazione 1/√N → robustezza al
  prior.
- **Confronto toy vs Geant4 sullo stesso motore**: misura della circolarità del toy.

## Fase 7 — Caso B

Direzione Ω_n su griglia di Fibonacci con raffinamento locale; campo lontano.

## Fase 8 — Caso C

Due stadi (Ω_n come B, poi griglia 2D su μ_E, σ_E); mai griglia 4D. Verifica della troncatura al
dominio.

## Fase 9 — Realismo

σ_np(E) nella verosimiglianza (spettro rivelato ≠ incidente), soglia su E_p, eventi n–C nel
modello, sorgenti multiple e fondo (mistura con frazione f), spettri non gaussiani.

## Fase 10 — Produzione batch

Run su cluster (CNAF) se i volumi lo richiedono; unione dei file per thread solo se serve.
