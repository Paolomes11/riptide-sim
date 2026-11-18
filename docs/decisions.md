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

## D10 — Mondo nel vuoto (2026-11-17)

- **Decisione:** il volume mondo è di `G4_Galactic`; l'unico materiale denso è il blocco di
  scintillatore.
- **Motivo:** il neutrone non deve interagire prima del blocco, così la verità della prima
  interazione è sempre quella nello scintillatore e lo spettro incidente è quello della sorgente.
- **Alternative:** aria (interazioni rare ma presenti, che sporcano la verità senza servire alla
  ricostruzione).
- **Conseguenze:** la geometria reale (fase 3) riapre la questione: supporti e materiali attorno
  allo scintillatore andranno introdotti in modo esplicito.

## D11 — Semi per evento (2026-11-17)

- **Decisione:** all'inizio di `GeneratePrimaries` il motore del thread viene riseminato con
  (event_id, seme del run) tramite `MixMaxRng::setSeeds` (stream derivati con
  `seed_uniquestream`). Il seme del run è al più 2³¹−1, perché `event_id` e i due semi stanno in
  colonne Int32 dell'ntupla. Prima di riseminare si invalida la cache delle sezioni d'urto
  adroniche del neutrone con una valutazione fittizia a 0.5·E_n,min.
- **Motivo:** ogni evento si riproduce da solo e l'output non dipende dal numero di thread. La
  sola risemina non bastava: le sezioni d'urto elastiche HP integrano il Doppler con un Monte
  Carlo che consuma numeri casuali, e `G4CrossSectionDataStore` riusa l'ultimo valore se energia,
  materiale e particella non cambiano. Il numero di estrazioni di un evento dipendeva quindi
  dall'evento precedente sullo stesso thread: senza invalidazione circa il 3% degli eventi
  differiva fra run a 1 e a 2 thread.
- **Alternative:** i semi per thread di Geant4 (riproducibili solo a parità di thread);
  `G4ParticleHPManager::SetNeglectDoppler` (scartata perché cambia la fisica).
- **Conseguenze:** ogni cambio di fisica, geometria o versione di Geant4 deve ripassare il test
  d'integrazione `sim_integration`, che confronta colonna per colonna un run a 1 thread e uno a 2.

## D12 — Versioni delle dipendenze (2026-11-17)

- **Decisione:** Geant4 11.4.1 MT con `QGSP_BIC_HP` (dati G4NDL); nlohmann/json 3.11.3 e Lyra
  1.6.1 copiati in `external/` (header-only, con sorgente, sha256 e licenza in ciascun `VERSION`);
  Catch2 3 dal sistema (oggi 3.4.0).
- **Motivo:** build riproducibile senza rete; le versioni finiscono in `run.json`
  (`geant4_version`, `git_describe`), quindi ogni run dichiara con cosa è stato prodotto.
- **Alternative:** `FetchContent` al configure (richiede rete e sposta la verifica degli hash).
- **Conseguenze:** aggiornare una dipendenza è una decisione esplicita: nuovo `VERSION`, `NOTICE`
  aggiornato, `sim_integration` e verifiche della fase 1b rilanciati.

## D13 — Valutazioni di riferimento per la fase 1b (2026-11-17)

- **Decisione:** le verifiche fisiche della fase 1b confrontano gli eventi Geant4 con due
  valutazioni: **JEFF-3.3** come verifica dell'implementazione e **ENDF/B-VIII.0** come secondo
  confronto, indipendente per il carbonio ma non per H-1 (vedi Motivo). Servono i file di H-1,
  C-12 e C-13: da MF3 le sezioni d'urto, e da MF4/MT2 di H-1 i coefficienti di Legendre della
  distribuzione angolare elastica nel CM. I file stanno in
  `<radice>/reference/{endfb80,jeff33}/` con URL e sha256 in `docs/architecture/data_layout.md`,
  e non entrano nel repository. JEFF-3.3 non ha un file C-12, ma solo C-nat (C-0).
- **Motivo:**
  - G4NDL 4.7.1 è G4NDL 4.6 con modifiche ad argon e scattering termico, e G4NDL 4.6 è ricavato
    da JEFF-3.3 (README dei dati G4NDL). Un accordo con JEFF-3.3 dice quindi che Geant4 legge e
    campiona correttamente i suoi dati, non che i dati siano giusti.
  - **ENDF/B-VIII.0 non è indipendente per H-1.** In entrambe le librerie H-1 è la valutazione
    LANL di G.M. Hale dall'analisi R-matrix del sistema N-N, con lo stesso testo descrittivo:
    in JEFF-3.3 è la versione EVAL-OCT05, in ENDF/B-VIII.0 la EVAL-JUL16. Lo scarto fra le due
    misura quindi la revisione di una stessa analisi, e non stima l'incertezza sui dati di H-1:
    gli errori comuni all'analisi non vi compaiono. Valori a confronto: σ_el differisce dello 0.27%
    a 10 MeV e < 0.2% a 15–20 MeV; a_1 vale −1.06e-2 (ENDF/B-VIII.0) contro −1.15e-2
    (JEFF-3.3) a 10 MeV, e circa −1.6e-2 per entrambe a 20 MeV.
  - Per il carbonio le valutazioni sono distinte: in JEFF-3.3 C-nat è di LANL/ORNL (Chadwick,
    Young, Fu, EVAL-JUN96) e C-13 viene da TENDL-2015; in ENDF/B-VIII.0 C-12 e C-13 sono
    valutazioni LANL del 2015 (Hale e altri).
  - Per il carbonio G4NDL ha `6_nat_Carbon` e `6_13_Carbon`, ma non `6_12`. Quando manca il file
    dell'isotopo, `G4ParticleHPNames::GetName` prova prima il file `Z_nat`, quindi a C-12 Geant4
    applica i dati del carbonio naturale, segnalandolo come abbondanza naturale. Il confronto per
    il carbonio si fa perciò sulla miscela del materiale (n_C12 σ_C12 + n_C13 σ_C13), non
    isotopo per isotopo.
  - ENDF/B-VIII.0 viene da NNDC (BNL). IAEA-NDS e OECD-NEA rispondono 403 agli accessi
    automatici, quindi JEFF-3.3 è scaricato a mano dal browser da IAEA-NDS.
- **Formato verificato** (manuale ENDF-102, §4.1–4.2):
  - in MF4 vale f(μ, E) = Σ_l (2l+1)/2 a_l(E) P_l(μ), con a_0 = 1 implicito, quindi
    ⟨P_l(μ)⟩ = a_l;
  - per H-1 in ENDF/B-VIII.0 si ha LTT = 1 (solo Legendre), LCT = 2 (CM), 153 energie,
    interpolazione lineare-lineare e NL = 6;
  - μ è il coseno CM del neutrone e vale E_p = E_p,max(1 − μ)/2 per masse qualsiasi
    (`np_kinematics.hpp`).
- **Alternative:** solo JEFF-3.3 (verifica circolare dei dati); solo ENDF/B-VIII.0 (non separa
  un errore di implementazione da una differenza fra valutazioni); risultati sperimentali (meno
  omogenei, rimandati alla validazione del rivelatore).
- **Conseguenze:** un cambio di versione di G4NDL richiede di ricontrollare la provenienza; se
  smette di derivare da JEFF-3.3, va cambiato il riferimento "di implementazione".
  Per H-1 la fase 1b verifica quindi l'implementazione (JEFF-3.3) e mostra di quanto si spostano
  i numeri fra due revisioni della stessa analisi (ENDF/B-VIII.0). Un confronto davvero
  indipendente per n-p (altra valutazione o analisi in onde parziali dei dati sperimentali) è
  rimandato a una fase successiva, e i risultati della 1b non vanno presentati come tale.

## D14 — Massa del bersaglio HP per H-1 nelle verifiche (2026-11-18)

- **Fatti** (verificati sui sorgenti di Geant4 11.4.1 e su G4NDL 4.7.1):
  - `G4ParticleHPElasticFS` costruisce il bersaglio con massa AWR · m_n, dove AWR è letto dalla
    prima riga del file dello stato finale. Per H-1 (`Elastic/FS/1_1_Hydrogen.z`) la riga è
    `1 0.999167 2` (rappresentazione, AWR, angoli nel CM), quindi M_T = 938.7827 MeV. È la
    massa dell'atomo di idrogeno, elettrone compreso, e non quella del protone.
  - Il rinculo è creato come `G4Proton` con l'impulso della cinematica a massa M_T, quindi la sua
    energia cinetica è √(p² + m_p²) − m_p. Ne segue un E_p più alto della cinematica fisica
    n–p di (5.37–5.49)·10⁻⁴ in relativo (≈ m_e/m_p) su tutto l'intervallo angolare.
  - A 10–20 MeV il bersaglio è a riposo: il moto termico (SVT) si applica solo sotto circa
    400 kT. Quindi u > 1 e θ_p > π/2 sono cinematicamente impossibili.
  - Questo modello riproduce i residui dello smoke test entro 1.2·10⁻¹⁴ MeV.
- **Decisione:** le verifiche della fase 1b modellano esattamente questa convenzione
  (`sim/validation/hp_recoil.hpp`). u = (1 − μ)/2 viene da E_p con la massa HP, μ da θ_p con la
  stessa massa, e il residuo di E_p è calcolato contro la cinematica HP, quindi è un oracolo a
  livello di arrotondamento. M_T è letto dall'intestazione di G4NDL a ogni verifica, non scritto
  a mano, e `checks.json` registra file, sha256 e AWR. Lo scarto dalla cinematica fisica è
  riportato a parte (`Ep_physical_offset`). u > 1 e θ_p > π/2 passano fra i contatori che
  devono essere zero.
- **Motivo:** con la normalizzazione a m_p, u arriva a 1.00054. Gli eventi oltre 1 sarebbero
  esclusi dall'ultimo bin, e alla statistica di produzione questo darebbe un difetto di circa
  l'1% contro un errore statistico di circa 0.3%: un falso fallimento del χ². Il residuo di E_p,
  inoltre, misurerebbe la convenzione e non un errore.
- **Alternative:** tolleranza sul residuo e su u (nasconde l'effetto invece di spiegarlo);
  correggere Geant4 (fuori perimetro, e la simulazione resta quella di riferimento).
- **Conseguenze:**
  - `common/np_kinematics.hpp` resta la cinematica fisica, usata dall'analisi.
  - **Aperto, per la fase 2:** la verità di Geant4 ha E_p più alto di circa 5.4·10⁻⁴ in
    relativo, circa 7.6 keV a E_p = 14 MeV. Un E_n ricostruito con la cinematica fisica avrebbe
    un bias relativo dello stesso ordine. È piccolo rispetto a σ_Ep, ma visibile nei test di bias
    ad alta statistica. Va deciso prima dell'analisi se correggerlo, se modellarlo o se
    riportarlo come sistematica nota.
  - Un cambio di G4NDL che modifichi l'AWR di H-1 si vede in `checks.json` e nello smoke test,
    che controlla AWR = 0.999167.
