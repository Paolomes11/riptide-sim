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

## D15 — Carbonio HP: inelastico senza cattura e composizione isotopica (2026-11-18)

- **Fatti** (verificati sui sorgenti di Geant4 11.4.1, su G4NDL 4.7.1 e su JEFF-3.3; oracolo
  permanente in `sim/tests/hp_data_reference_test.cpp`):
  - La sezione d'urto inelastica HP del carbonio esclude la cattura. `Inelastic/CrossSection/
    6_nat_Carbon` parte dalla soglia (4.812 MeV) e coincide con JEFF-3.3 C-0 MT3 − MT102.
    L'integrale su 10–20 MeV differisce di +2.5·10⁻⁷, contro −2.67·10⁻⁴ rispetto a MT3, cioè
    proprio −σ_cap/MT3. Su alcuni nodi lo scarto arriva a 6.2·10⁻⁴, con origine non verificata.
    L'elastico e la cattura di H-1, C-0 e C-13 sono copiati da JEFF-3.3.
  - G4NDL non ha `6_12_Carbon`: `G4ParticleHPNames::GetName` ripiega su `6_nat_Carbon`. L'elemento
    è Σ abbondanza · σ(file), quindi σ_C = 0.9893 · σ(C-0) + 0.0107 · σ(C-13). Il C-13 entra due
    volte, una dentro C-0 e una come isotopo a sé.
  - L'inelastico di `6_13_Carbon` è la somma dei canali parziali di JEFF-3.3 C-13 entro 2.7·10⁻⁴,
    ma non è MT3 − MT102: ne differisce fino a 1.84·10⁻², perché in JEFF-3.3 C-13 MT4 non è la
    somma di MT51–91.
  - Nelle tabelle di G4NDL un'energia ripetuta segna una discontinuità, come in ENDF: in C-13
    l'elastico e la cattura ripetono 6.9501 MeV, la cattura anche 30 MeV. Il lettore
    (`sim/validation/hp_data.hpp`) le accetta e vi prende il valore sopra; tre energie uguali o
    una decrescente sono un errore.
  - L'oracolo legge G4NDL da `G4NEUTRONHPDATA` e JEFF-3.3 da `<root>/reference/jeff33`, e viene
    saltato se manca l'uno o l'altro.
- **Decisione:** l'attesa `jeff33` è composta come la compone il modello HP: H-1, C-0 con
  densità n_C12 e C-13 con densità n_C13. Verifica quindi l'implementazione, non la fisica di C-0.
  `endfb80` resta H-1 + C-12 + C-13. `checks_version` passa a 3.
- **Motivo:** con C-0 a tutta la densità del carbonio, il non elastico atteso sul carbonio era
  basso di circa 5.3·10⁻³ in relativo a 14 MeV (0.724 b per C-13 contro 0.483 b per C-0). A
  10⁷ eventi l'errore statistico è circa 8·10⁻⁴: un falso fallimento di circa 6σ, dello stesso
  tipo di D14. Sullo smoke test la frazione attesa nc_inelastic sale di 4.5·10⁻³ in relativo.
- **Alternative:** modellare l'inelastico di C-13 come somma dei parziali, che richiederebbe una
  tabella non elastica nel modello dei nuclidi delle verifiche (rimandato: l'effetto è sotto la
  statistica di produzione); confrontare con C-0 a tutta densità e una tolleranza (nasconde
  l'effetto).
- **Conseguenze:**
  - Resta un residuo noto e non modellato di circa 0.0107 × 0.0130 / 0.486 ≈ 2.9·10⁻⁴ sul non
    elastico del carbonio vicino a 14 MeV, da riportare nelle previsioni della fase 1b.
  - **Per la fase 2:** il doppio conteggio del C-13 è un artefatto di Geant4 (circa +0.5% sul non
    elastico del carbonio) e va tenuto fra le sistematiche del fondo da carbonio.

## D16 — Bordo a 20 MeV: elastico HP e CHIPS sovrapposti in 19.5–20 MeV (2026-11-19)

- **Fatti** (verificati sui sorgenti di Geant4 11.4.1):
  - `G4ParticleHPElastic` vale da 0 a 20 MeV (`SetMaxEnergy(20. * MeV)`).
    `G4HadronElasticPhysicsHP::ConstructProcess` abbassa il minimo del modello elastico del
    neutrone già registrato (`G4ChipsElasticModel`, da `G4HadronElasticPhysics`) a 19.5 MeV e
    poi registra HP. In 19.5–20 MeV i due modelli si sovrappongono.
  - Con due modelli applicabili, `G4EnergyRangeManager::GetHadronicInteraction` ne sceglie uno a
    caso, evento per evento, con una rampa lineare: P(CHIPS) = (E_n − 19.5 MeV)/(0.5 MeV),
    P(HP) = 1 − P(CHIPS). Sotto 19.5 MeV c'è solo HP.
  - Le sezioni d'urto restano quelle HP fino a 20 MeV (`G4ParticleHPElasticData` è aggiunto per
    ultimo e prevale). La sovrapposizione cambia quindi solo lo stato finale, non le frazioni
    della prima interazione né la profondità.
  - L'inelastico segue la stessa logica: BIC da 19.9 MeV (`G4HadronPhysicsQGSP_BIC_HP`), HP fino
    a 20 MeV (`G4NeutronPHPBuilder`). Cambia solo lo stato finale non elastico, che le verifiche
    non usano.
  - `G4ChipsElasticModel` usa la cinematica di `G4HadronElastic::ApplyYourself`: proiettile con
    la massa PDG del neutrone, bersaglio con `G4NucleiProperties::GetNuclearMass(1, 1)`, cioè la
    massa del protone, e rinculo `G4Proton`. Le masse CLHEP coincidono con quelle di
    `riptide::kinematics`, quindi un evento CHIPS segue la cinematica fisica n–p, non quella HP
    di D14. La distribuzione angolare è quella di CHIPS, non JEFF-3.3.
- **Decisione:** `riptide_sim_check` classifica ogni evento n–p dalla cinematica, senza nuove
  colonne nello schema. È HP se E_p coincide con la cinematica HP entro 10⁻¹¹·E_n, non HP se
  coincide con la cinematica fisica, ambiguo se con entrambe (trattato come HP, in
  `expected_rare`). Se non coincide con nessuna, o se è non HP sotto 19.5 MeV, va fra i contatori
  `must_be_zero`. Spettro u, Legendre e residui di E_p usano solo gli eventi HP, che seguono
  ancora JEFF-3.3. Gli eventi non HP formano un oracolo a sé (`model_mix`): conteggio contro
  Σ p(E_n) con pull binomiale e χ² su 10 bin di E_n nella finestra, più il loro residuo contro la
  cinematica fisica. `checks_version` passa a 4.
- **Motivo:** la tolleranza viene dallo smoke test a 14 MeV (residuo HP massimo 1.24·10⁻¹⁴ MeV),
  ben sotto lo scarto fra le due cinematiche (circa 5.4·10⁻⁴ in relativo). Senza
  classificazione, nel run piatto 10–20 MeV poco meno del 2.5% degli eventi n–p (media della
  rampa sulla sorgente, pesata con σ_np decrescente) sarebbe stato CHIPS:
  l'oracolo dei residui sarebbe fallito per costruzione e spettro u e Legendre avrebbero
  mescolato due modelli angolari.
- **Alternative:** tagliare E_n > 19.5 MeV nelle verifiche (perde il controllo del bordo);
  salvare il modello in una nuova colonna (cambio di schema, non necessario perché la cinematica
  lo distingue); tolleranza larga sui residui (nasconderebbe errori veri).
- **Conseguenze:**
  - Il run di bordo 19.5–20 MeV verifica la rampa. Gli eventi CHIPS ambigui attesi sono
    trascurabili (circa 0.06 su 1.5·10⁶ eventi n–p), ma un evento ambiguo abbasserebbe
    `Ep_physical_offset.relative_min`.
  - Lo scarto relativo dalla cinematica fisica dipende da E_n: (5.37–5.49)·10⁻⁴ di D14 vale a
    14 MeV e scende a circa 5.33·10⁻⁴ vicino a 20 MeV. Le previsioni lo calcolano per ogni run.
  - **Per la fase 2:** sopra 19.5 MeV lo stato finale elastico di Geant4 non è più solo JEFF-3.3.
    Il bordo della sorgente a 20 MeV va trattato come sistematica o evitato.

## D17 — Esito c1 della fase 1b: H-2 nelle attese e offset al netto dell'arrotondamento (2026-11-19)

- **Fatti** (run `1b_*_1e7`, verifiche `c1` con `checks_version` 4, soglie delle previsioni
  committate prima dei run):
  - piatto 10–20 MeV e bordo 19.5–20 MeV: PASS; mono 14 MeV: FAIL su due verifiche;
  - **offset relativo massimo, mono:** osservato 5.4861409·10⁻⁴ contro il massimo previsto
    5.4861387·10⁻⁴ (margine 10⁻⁹ relativo). L'evento che lo dà ha E_p ≈ 9 eV e residuo
    E_p − E_p,HP = 2·10⁻¹⁵ MeV: un arrotondamento, che diviso per pochi eV sposta lo scarto
    relativo di ~2·10⁻¹⁰. La sola parte HP dello stesso evento vale 5.4861386·10⁻⁴, dentro il
    range. L'errore è nella formulazione dell'oracolo, non nella simulazione;
  - **eccesso di "other":** pull di +4.8 (mono), +2.9 (piatto), +2.5 (bordo). Nel mono:
    395 osservati contro 310.4 attesi, eccesso 85 ± 18 (c). Il break-up n–d su H-2 (Z = 1,
    anelastico) finisce in "other" nella simulazione, ma le attese contenevano solo H-1 e C: le
    previsioni avevano stimato l'effetto del deuterio solo sul conteggio n–p. Con le sezioni
    d'urto di G4NDL a 14 MeV (σ_inel(H-2) = 0.1659 b, σ_el(H-1) = 0.6871 b) e n(H-2)/n(H-1) =
    1.15·10⁻⁴ l'aggiunta attesa è di circa 64 eventi.
- **Decisione:**
  - H-2 entra nelle attese con le sezioni d'urto elastica, anelastica e di cattura di **G4NDL
    4.7.1**, cioè gli stessi dati che usa HP, in entrambe le librerie di riferimento. Per H-1
    sarebbe circolare; per una componente di 1.15·10⁻⁴ basta dichiararlo: la verifica di H-2 **non
    è indipendente**. Le categorie seguono la simulazione: l'elastico su H-2 conta come n–p (Z = 1)
    ma non ha protone di rinculo, break-up e cattura vanno in "other". Il peso n–p HP (spettro u,
    Legendre, `model_mix` previsto) usa solo l'elastico su H-1. Le tabelle G4NDL hanno
    discontinuità vere (punto doppio a 14.25 MeV in elastico e cattura), riportate nel totale.
  - L'oracolo dell'offset confronta col range previsto gli estremi `rounding_free_relative_*`:
    a ogni evento si toglie R/E_p,fis, con R = `Ep_residuals.max_abs` del run (già vincolato a
    < 10⁻⁹·E_n,max dal suo oracolo). Il confronto è una condizione necessaria sull'offset esatto e
    non può più fallire per arrotondamento; `relative_min/max` restano riportati.
  - `checks_version` passa a 5. Soglie, α, margine e range previsti **non cambiano**. Come previsto
    al §6 delle previsioni, `riptide_sim_check` si rilancia (`c2`) sugli **stessi eventi** e si
    riportano entrambi gli esiti.
- **Previsioni del rilancio** (scritte prima di `c2`; JSON in
  `docs/validation/predictions/*_1e7_v5.json`, stessi comandi delle previsioni originali):
  - attese di "other" (jeff33): mono 310.4 → 374.4, piatto 500.1 → 562.3, bordo 776.2 → 839.9;
    le altre attese cambiano di ≤ 10⁻⁴ in relativo, i coefficienti di Legendre attesi di
    ≤ 6·10⁻⁹, piano dei test e range dell'offset sono identici (a);
  - pull delle frazioni attesi in `c2`, stimati dai conteggi di `c1` con le nuove attese (d):

    | Run | unknown | n–p | n–C el. | n–C inel. | other | χ² frazioni |
    |---|---|---|---|---|---|---|
    | mono | +1.0 | +0.5 | +0.2 | −2.3 | +1.1 | ≈ 6/4, p ≈ 0.2 |
    | piatto | 0.0 | −1.2 | +0.9 | +0.3 | +0.1 | ≈ 2/4, p ≈ 0.7 |
    | bordo | 0.0 | +0.3 | −0.3 | 0.0 | +0.2 | ≈ 0.2/4, p ≈ 0.99 |

  - spettro u, profondità e Legendre: stessi esiti di `c1` entro ~10⁻³ in χ² (le attese cambiano
    solo al livello di 10⁻⁴); `model_mix` invariato;
  - offset: `rounding_free_relative_max` del mono sotto 5.4861387·10⁻⁴; tutti e tre i run PASS.
- **Alternative:** un'altra valutazione di H-2 (indipendente, ma senza un file pronto nel formato
  già letto e per un effetto di 10⁻⁴); allargare il margine dell'offset (nasconderebbe errori veri
  a E_p alti); tagliare gli eventi a E_p piccolo (perde proprio il bordo dello spettro).
- **Conseguenze:** l'esito della fase 1b è quello di `c2`, con `c1` riportato accanto. Il
  confronto n–p indipendente resta quello di D13; H-2 è controllato solo per coerenza con G4NDL.
