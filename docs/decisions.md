# Registro delle decisioni

Formato: **Decisione** · **Motivo** · **Alternative** · **Conseguenze**. Le decisioni non si
cancellano: se cambiano, se ne aggiunge una nuova che rimanda alla vecchia.

## D1 — C++17 (2026-11-16)

- **Decisione:** tutto il progetto in C++17, estensioni GNU disattivate.
- **Motivo:** ROOT 6.36 sulla macchina di sviluppo è compilato in C++17; mescolare standard con
  ROOT causa problemi di ABI e di dizionari. Geant4 11.4 accetta C++17.
- **Alternative:** C++20 (concepts, `std::span`): richiederebbe ricompilare ROOT.
- **Conseguenze:** niente `std::span`/`std::format`; si usano puntatore + lunghezza o viste
  minimali scritte a mano.

## D2 — Output in TTree, lettura con RDataFrame (2026-11-16)

- **Decisione:** la simulazione scrive TTree con `G4AnalysisManager`, un file per thread senza
  unione; l'analisi legge con `RDataFrame`.
- **Motivo:** in Geant4 11.4 la scrittura RNTuple non è supportata da `G4AnalysisManager`;
  `RDataFrame` legge già anche RNTuple, quindi la migrazione futura tocca solo la sim. Il file per
  thread evita la serializzazione in fase di unione.
- **Alternative:** unione automatica dei file (costo in memoria e tempo); HDF5 (fuori
  dall'ecosistema Geant4/ROOT).
- **Conseguenze:** l'ordine degli eventi non è deterministico in MT ⇒ ogni evento porta
  `event_id` e i due semi del suo RNG (`rng_seed_0/1`), sufficienti a ripeterlo.

## D3 — Libreria statistica senza ROOT né Geant4 (2026-11-16)

- **Decisione:** `analysis/stats` è C++ puro (solo libreria standard); I/O ROOT confinato in
  `analysis/io`.
- **Motivo:** test veloci e portabili, CI senza installare Geant4/ROOT, motore riusabile su dati
  non ROOT.
- **Alternative:** usare RooFit/Minuit direttamente.
- **Conseguenze:** gli adattatori ROOT → struct-of-arrays stanno in `analysis/io`; un eventuale
  minimizzatore esterno entra dietro un'interfaccia.

## D4 — Nessun `-ffast-math` (2026-11-16)

- **Decisione:** vietato `-ffast-math` / `-ffinite-math-only`; `numerics.hpp` dà `#error` se
  attivi.
- **Motivo:** le log-verosimiglianze valgono legittimamente −∞ fuori dal supporto; fast-math
  permette al compilatore di assumere assenza di NaN/Inf e di riordinare le somme (risultati non
  riproducibili).
- **Alternative:** fast-math solo su funzioni selezionate.
- **Conseguenze:** le ottimizzazioni passano da layout dei dati, vettorizzazione e parallelismo;
  `-march=native` resta opzionale (`RIPTIDE_ENABLE_NATIVE`).

## D5 — Dati fuori dal disco di sistema (2026-11-16)

- **Decisione:** radice dei dati = CLI `--data-dir` > `$RIPTIDE_DATA_DIR` >
  `/mnt/external_ssd/riptide-sim`; la radice deve esistere, essere scrivibile, **non** stare sul
  dispositivo di `/` e avere spazio per la stima del run + 2 GiB. Nessun fallback, la radice non
  viene mai creata dal programma.
- **Motivo:** il disco di sistema è quasi pieno (97%); un run che lo riempie blocca la macchina.
  Creare la radice mascherebbe un SSD non montato (si scriverebbe nel mount point vuoto, su `/`).
- **Alternative:** cartella `data/` nel repo (ignorata da git).
- **Conseguenze:** senza SSD montato i run falliscono con un messaggio esplicito; per test brevi
  esiste l'opzione esplicita `allow_system_disk`.

## D6 — Niente ottica (2026-11-16)

- **Decisione:** nessun trasporto di fotoni, lenti o fotocatodo.
- **Motivo:** l'obiettivo è ricostruire il neutrone da ciò che interagisce nello scintillatore;
  l'ottica appartiene a un altro progetto e complicherebbe solo l'analisi.
- **Alternative:** catena completa fino al fotocatodo.
- **Conseguenze:** il realismo del rivelatore passa da deposito + quenching di Birks e da una
  risposta (σ_Ep, σ_θ) applicata in analisi.

## D7 — `truth_valid_np` nello schema (2026-11-16)

- **Decisione:** colonna intera `truth_valid_np` = 1 se le colonne del protone di rinculo sono
  riempite (prima interazione n–p), 0 altrimenti.
- **Motivo:** i valori sentinella (es. E_p = −1) si confondono con dati veri nei tagli e negli
  istogrammi.
- **Alternative:** NaN come sentinella (incompatibile con il divieto di trattare NaN come dato).
- **Conseguenze:** ogni selezione dell'analisi parte da `truth_valid_np == 1`.

## D8 — CLI e configurazione rimandate alla fase 1 (2026-11-16)

- **Decisione:** il parser della riga di comando (lyra) e il lettore JSON (nlohmann/json) entrano
  in `external/` solo nella fase 1, quando servono davvero.
- **Motivo:** nessuna dipendenza senza un uso che la giustifichi.
- **Conseguenze:** in fase 0 `riptide_sim` accetta solo `--version`.

## D9 — Risposta del rivelatore in analisi (2026-11-16)

- **Decisione:** la simulazione salva verità e depositi; lo smearing (σ_Ep, σ_θ) si applica in
  analisi con seme registrato.
- **Motivo:** la stessa simulazione (costosa) serve a più ipotesi di risoluzione.
- **Conseguenze:** l'analisi deve registrare il seme dello smearing nei suoi metadati.
