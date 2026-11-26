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
  - **Violazione del bilancio energetico (2026-11-20).** L'offset non è solo uno spostamento di
    E_p: ogni urto n–p HP crea l'energia ε = E_p(m_p) − E_p(M_T), che finisce nel deposito.
    Lungo una termalizzazione la somma è Σε ≈ 5.4·10⁻⁴ · ΣE_p, e ΣE_p arriva a circa E_n.
    Questo spiega la banda di `det_Edep` a circa +7.5 keV sopra E_n, vista nel run
    d'integrazione a 14 MeV (28 eventi su 1000 sopra E_n).
    - **Metodo.** Un diagnostico fuori dal repo replica `riptide_sim` con le stesse azioni,
      avvolte. Lo stesso seme, lo stesso config e lo stesso binario debug riproducono il run
      d'integrazione bit per bit su tutte le colonne dei 1000 eventi. Per ogni step nello
      scintillatore il diagnostico calcola dK = K_post + ΣK_sec + E_dep − K_pre. Per gli step
      adronici calcola anche Q dalle masse nucleari di Geant4 (`G4NucleiProperties`) e ε = dK − Q.
    - **Bilancio per evento (c).** E_n + ΣdK − E_uscita − `det_Edep` = 0 entro 4·10⁻¹⁵ MeV su
      tutti i 1000 eventi, e il deposito ricostruito coincide con quello di `riptide_sim`.
    - **Per urto n–p (c).** Sopra 400 kT (1017 urti), ε misurato e ε di D14 coincidono entro
      4.5·10⁻¹³ MeV. ε/E_p vale (5.36–5.44)·10⁻⁴ per E_p > 1 keV. Sotto 1 keV i valori sono
      dominati dall'arrotondamento: un urto a E_p = 28 eV dà 5.57·10⁻⁴. Sui 505 urti termici
      Σε = 2·10⁻⁵ MeV, trascurabile.
    - **I 27 eventi della banda (c), un solo run.**
      - **Fatti:** eccesso fra 2.2 e 11.8 keV, media 6.9 keV, sd 1.8 keV; Σε_np fra 5.9 e
        7.6 keV, media 7.43 keV, sd 0.31 keV.
      - **Scomposizione:** eccesso = Σε_np + Σε_C − E_uscita, più (in un solo evento)
        l'eccesso degli ioni descritto sotto. Il termine dominante è Σε_np in tutti i 27
        eventi.
      - **Il carbonio:** l'elastico HP sul C dà anch'esso ε > 0, ma piccolo: fino a 1.1 keV
        per evento, 0.1 MeV in totale su 649 urti.
      - **La dispersione** viene da E_uscita. In 26 eventi il neutrone, rallentato da 2–47 urti
        n–p, esce dal blocco senza essere catturato, con un'energia fra 5 meV e 5.3 keV.
      - **I due eventi con cattura nel blocco.**
        - Evento 958, nella banda: il γ esce tutto.
        - Evento 56, a 15.6 MeV: dei 2.224 MeV di Q restano 1.613 MeV, perché il γ esce con
          0.611 MeV.
    - **Previsione sbagliata:** la cattura n–p conserva l'energia. Mi aspettavo, dalla lettura
      di `G4NeutronHPCaptureFS`, circa +1.3 keV per cattura dal rinculo del deuterone. Invece
      K_d (1.31–1.33 keV) ed E_γ sommano a Q = m_n + M(¹H) − M(²H) entro 10⁻⁷ MeV, cioè entro
      l'energia del neutrone termico.
    - **Altro trovato (c):** a fine percorso `ionIoni` deposita più del K iniziale dello step.
      - **Dove:** in 56 rinculi C12 (su 486 sopra 0.1 MeV) e in 4 Be9, con K_pre fra 0.11 e
        0.40 MeV.
      - **Quanto:** l'eccesso relativo va da 0.2% a 15%, mediana 4.4%. In totale fa 0.74 MeV su
        492 MeV di rinculi C12.
      - **Perché:** la causa non è verificata sui sorgenti (d). È dello stesso tipo del
        G4NuclearStopping del protone in D19, ma relativamente più grande. Pesa poco sulla luce,
        perché i rinculi di C sono fortemente spenti da Birks, ma entra in `det_Edep`.
      - **Le conversioni γ → e⁺e⁻ e le annichilazioni** danno ±2m_e: per gli step non adronici
        il diagnostico non calcola Q, quindi non sono violazioni.
    - **Conseguenza:** `det_Edep` > E_n non è un errore del codice, ma il segno della convenzione
      HP. Le verifiche dei depositi devono tollerare Σε ≲ 6·10⁻⁴ · E_n senza un γ di cattura
      depositato.

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

## D18 — Costante di Birks per lo scintillatore PVT (2026-11-19)

- **Fatti:**
  - in Geant4 11.4.1 `G4EmSaturation` applica Birks passo per passo, con kB letto da
    `G4IonisParamMat::GetBirksConstant()`;
  - i valori di serie esistono solo per pochi materiali (fra cui `G4_POLYSTYRENE`), non per
    `G4_PLASTIC_SC_VINYLTOLUENE`, per cui senza un valore esplicito la luce visibile coincide con
    l'energia depositata.
- **Decisione:**
  - **kB = 0.155 mm/MeV**, con incertezza ±0.005 mm/MeV (intervallo di credibilità al 68%, che
    comprende statistica e sistematica).
  - **Fonte:** T. Pöschl, D. Greenwald, M. J. Losekamm, S. Paul, *Measurement of ionization
    quenching in plastic scintillators*, Nucl. Instrum. Methods Phys. Res. A 988 (2021) 164865,
    doi:10.1016/j.nima.2020.164865 (arXiv:2007.08366), Tabella 1, riga "Birks", colonna BC-408.
  - **Misura:** BC-408, a base di poliviniltoluene, con protoni da 30–100 MeV che si fermano in
    fibre scintillanti; Bragg curve fino a circa 500 MeV/cm; stopping power CSDA da NIST.
  - Gli autori dichiarano accordo con il valore BC-408 di M. Almurayshid et al., J. Appl. Clin.
    Med. Phys. 18 (5) (2017) 210–219, doi:10.1002/acm2.12143, valore non letto direttamente qui.
  - BC-408 ed EJ-200 sono formulazioni commerciali equivalenti; il materiale Geant4 resta
    `G4_PLASTIC_SC_VINYLTOLUENE`.
  - kB entra come parametro **obbligatorio** del config (`config_version` 2), con la fonte
    registrata in `run.json`; kB = 0 resta ammesso come oracolo (luce visibile = deposito).
- **Limiti noti** (dalla stessa fonte):
  - per BC-408 Birks è il modello preferito fra quelli confrontati, ma nessun modello descrive bene
    l'intero intervallo di dE/dx; le misure sullo stesso scintillatore spesso non concordano;
  - il fit copre dE/dx fino a circa 500 MeV/cm, mentre i protoni di rinculo a fine percorso
    superano questo valore: lì kB è un'estrapolazione;
  - un solo kB non descrive il quenching di ioni più pesanti (carbonio di rinculo, alfa).
- **Sistematica** (fissata il 2026-11-20, prima delle previsioni della fase 2):
  - **intervallo kB ∈ [0.145, 0.165] mm/MeV**, ±2σ del valore Birks di Pöschl (Tabella 1,
    BC-408, 0.155 ± 0.005, 68% di credibilità con statistica e sistematica insieme, fra cui un
    prior del 4% sullo stopping NIST). Viene **solo** da Pöschl;
  - nella stessa Tabella 1 il modello di Chou dà per BC-408 kB = 0.151 ± 0.04 mm/MeV con C
    compatibile con 0: dentro l'intervallo, ma con un errore otto volte più largo;
  - T. A. Laplace et al. (2020), letto il 2026-11-20 (arXiv:2009.07217): misura la resa luminosa
    dei protoni in EJ-200, EJ-204 ed EJ-208 relativa a un elettrone da 477 keV, **senza** fit di
    Birks né valori di kB. Un fit proprio su quei dati sonderebbe proprio la zona oltre
    500 MeV/cm; non è fatto;
  - non letti: Almurayshid 2017 (citato da Pöschl come in accordo, valore non riportato) e
    Hirschberg 1992;
  - punto aperto: l'intervallo non copre l'estrapolazione oltre 500 MeV/cm, dove arriva la fine
    del range di ogni protone.
- **Alternative:**
  - il kB di `G4_POLYSTYRENE` (materiale diverso);
  - valori da misure con elettroni o gamma, che non sondano l'alto dE/dx dei protoni;
  - modelli a più parametri (Chou, Voltz), che `G4EmSaturation` non implementa.
- **Conseguenze:**
  - la luce visibile simulata dipende dal passo del protone: la funzione di passo diventa un
    parametro registrato, con studio di convergenza nella fase 2 (superato da D21: per i protoni la
    luce non dipende più dal passo);
  - i PDF delle fonti restano fuori dal repo.

## D19 — Schema v2: depositi, traccia del protone e classe dell'evento (2026-11-19)

- **Fatti** (sorgente di Geant4 11.4.1):
  - `G4EmSaturation` riempie tabelle statiche condivise e non è thread-safe da costruire;
    `G4EmParameters::GetEmSaturation()` la crea al primo uso e attiva Birks;
  - i modelli HP anelastici uccidono il neutrone primario ed emettono il neutrone uscente come
    traccia nuova;
  - `G4NuclearStopping::AlongStepDoIt` limita la perdita del passo all'energia cinetica *prima*
    del passo, non a quella rimasta dopo la ionizzazione, e `G4ParticleChangeForLoss` porta poi
    l'energia a 0 depositando comunque l'intera perdita. Sull'ultimo passo di un protone che si
    ferma, la traccia può quindi depositare un po' più dell'energia con cui è nata: 1.9·10⁻⁴ MeV
    su un protone da 11.6 MeV in un run di prova a 14 MeV (perdita nucleare del passo
    1.3·10⁻³ MeV, energia residua 1.1·10⁻³ MeV).
- **Decisione:**
  - `G4EmSaturation` si crea una volta sul master in `DetectorConstruction::Construct()`, dopo
    materiali e kB; i worker la usano in sola lettura per `VisibleEnergyDepositionAtAStep`.
  - Il protone di rinculo è il secondario protone del passo del **primo** n–p (D7), marcato con
    un `G4VUserTrackInformation`; il ruolo si legge dal marcatore, mai confrontando energie.
  - `det_Edep` e `det_Evis` sommano i passi di **tutte** le particelle nello scintillatore;
    `truth_p_*` solo quelli della traccia marcata (i suoi secondari no).
  - Le colonne del vertice n–p e del protone valgono 0 se `truth_valid_np` = 0, per costruzione.
  - La classe dell'evento viene dalle interazioni adroniche del neutrone **primario** nello
    scintillatore, con precedenza dall'alto:

    | Codice | Classe | Condizione |
    |---|---|---|
    | 6 | altro | almeno un'interazione che non è n–p, n–C elastico o n–C anelastico |
    | 5 | misto | almeno un n–p e almeno un n–C |
    | 1 | singolo n–p | esattamente un n–p, nient'altro |
    | 2 | n–p multiplo | due o più n–p, nient'altro |
    | 4 | n–C anelastico senza n–p | almeno un n–C anelastico (n–C elastici ammessi) |
    | 3 | solo n–C elastico | solo n–C elastici |
    | 0 | nessuna | il primario attraversa senza interagire |

    I codici sono Int32 dello schema e non si riusano mai. Come in D17, n–p comprende lo
    scattering su H-2 (Z = 1).
  - Nei confronti fra `truth_p_Edep` e `truth_Ep` si ammette un margine di 0.01 MeV per l'effetto
    di `G4NuclearStopping` (molto sopra l'eccesso visto, molto sotto l'errore di MeV di una
    traccia sbagliata). La distribuzione dell'eccesso si misura nelle verifiche della fase 2.
- **Limite noto:** le interazioni dei neutroni secondari (per esempio un n–p dopo un n–C
  anelastico) non entrano nella classe; `truth_n_secondary_n_hadronic` dice quanto spesso
  accade.
- **Alternative:**
  - riconoscere il protone confrontando energia e direzione dei secondari, fragile;
  - creare `G4EmSaturation` in ogni worker, che non è sicuro con le tabelle statiche;
  - seguire nella classe anche i neutroni secondari, che richiede di marcare la catena delle
    tracce: rimandato finché il contatore non mostra che serve.
- **Conseguenze:** lo schema passa alla versione 2 (18 colonne in coda); il lettore accetta
  ancora la versione 1, quindi `riptide_sim_check` legge i run della fase 1b.

## D20 — Riferimento PSTAR per range e luce dei protoni (2026-11-20)

- **Dati:** tabella PSTAR del NIST per il materiale 216, "PLASTIC SCINTILLATOR
  (vinyltoluene-based)": stopping power elettronico, nucleare e totale, range CSDA e proiettato,
  132 energie da 1 keV a 10 GeV. Scaricata il 2026-11-20 da
  https://physics.nist.gov/cgi-bin/Star/apdata.pl (formato testo) e salvata in
  `reference/pstar/`, sola lettura, con origine e sha256 in `architecture/data_layout.md`.
- **Materiale a confronto:**

  | | PSTAR 216 | `G4_PLASTIC_SC_VINYLTOLUENE` (Geant4 11.4.1) |
  |---|---|---|
  | densità | 1.032 g/cm³ | 1.032 g/cm³ |
  | I | 64.7 eV | 64.7 eV |
  | composizione | H 0.085000, C 0.915000 in peso | C9H10 per numero di atomi |

  C9H10 con le masse atomiche standard dà circa 0.0853 di H in peso. Il sorgente di
  `G4NistMaterialBuilder.cc` lascia commentate le frazioni NIST e annota che non corrispondono al
  PVT puro: la differenza è nota a Geant4 e qui si dichiara, senza stimarne l'effetto a occhio.
- **Fatti** (sorgente di Geant4 11.4.1, `QGSP_BIC_HP` con `G4EmStandardPhysics_option4`):
  - per il protone `G4hIonisation` usa `G4BraggModel` sotto 2 MeV e `G4BetheBlochModel` sopra;
  - `G4BraggModel` legge `G4PSTARStopping`, e per questo materiale la tabella (60 punti da 1 keV
    a 2 MeV) coincide **valore per valore** con lo stopping elettronico della tabella scaricata;
  - `G4BetheBlochModel` calcola la formula di Bethe con correzione di densità, di shell e le
    correzioni di ordine superiore di `G4EmCorrections` (Barkas, Bloch, Mott), con la stessa I;
    i dati ICRU90 non si applicano (coprono solo acqua, aria e grafite);
  - nella tabella di dE/dx, `G4EmModelManager::FillDEDXVector` moltiplica il modello alto per
    1 + (S_Bragg(E_lim)/S_BB(E_lim) − 1)·E_lim/E, con E_lim = 2 MeV: sopra 2 MeV Geant4 resta agganciato a
    PSTAR, con un peso che decresce come 1/E;
  - lo stopping nucleare viene da un processo separato (`G4NuclearStopping`), non da PSTAR.
- **Decisione:**
  - PSTAR è il riferimento (c) per il range CSDA e per la luce attesa da Birks dei protoni
    contenuti senza interazioni adroniche:
    - ⟨`truth_p_length`⟩ contro R_CSDA(E_p);
    - ⟨`truth_p_Evis`⟩ contro L(E_p) = ∫₀^E_p dE / (1 + kB·S(E)), con S totale e kB di D18;
  - in entrambi si usa lo stopping **totale** (elettronico + nucleare), perché la simulazione
    deposita anche la perdita nucleare;
  - tra i punti della tabella si interpola in log-log;
  - le soglie di accordo si fissano nelle previsioni della fase 2, prima dei run.
- **Limite di indipendenza:**
  - sotto 2 MeV il confronto verifica solo l'implementazione (interpolazione, passo, raccordo fra
    i modelli), non la fisica: Geant4 usa la stessa tabella;
  - sopra 2 MeV il confronto è indipendente solo in parte: stessa I, correzioni implementate in
    modo diverso, aggancio a PSTAR a 2 MeV;
  - un accordo dice quindi che la simulazione trasporta il protone come previsto da questi dati,
    non che i dati descrivono bene il PVT reale.
- **Alternative:**
  - tabelle ICRU49 lette direttamente (stessa origine di PSTAR);
  - SRIM, con un modello diverso, ma senza una tabella pubblica con origine fissata come PSTAR;
  - il solo `G4EmCalculator`, che confronterebbe Geant4 con sé stesso.
- **Conseguenze:**
  - nuovo lettore della tabella in `sim/validation` (senza Geant4, compilato anche in `ci`);
  - la sezione `pstar` delle verifiche arriva con `checks_version` 7.

## D21 — Luce dei protoni continua e funzione di passo (2026-11-20)

- **Problema:** con `G4EmSaturation` Birks si applica step per step, eloss/(1 + kB·eloss/ℓ). La
  luce di un protone che si ferma dipende allora dalla funzione di passo. Il confronto PSTAR
  (D20) sullo smoke test dava +0.57% ± 0.04% di luce sopra L(E_p; kB).
- **Scan** (fuori dal repo, mono 14 MeV, 10⁴ eventi, 1 thread, seme 20260907, circa 2770 protoni
  contenuti senza interazioni adroniche per run; config uguali a `config/mono_14MeV.json` salvo
  `proton_step_function`). Residuo relativo medio di `truth_p_Evis` contro L di PSTAR e di
  `truth_p_length` contro R_CSDA, (c):

  | (dRoverRange, finalRange) | luce, Birks per step | luce, L continua (questa decisione) | range |
  |---|---|---|---|
  | (0.2, 0.1 mm) | +0.511% ± 0.012% | −0.062% ± 0.001% | −0.101% ± 0.029% |
  | (0.1, 0.1 mm) | +0.471% | −0.062% | −0.100% |
  | (0.05, 0.1 mm) | +0.447% | −0.062% | −0.103% |
  | (0.02, 0.1 mm) | +0.442% | −0.062% | −0.099% |
  | (0.01, 0.1 mm) | +0.439% | −0.062% | −0.099% |
  | (0.01, 0.01 mm) | −0.167% | −0.062% | −0.183% |
  | (0.01, 1 µm) | −0.392% ± 0.016% | −0.062% ± 0.001% | +0.065% ± 0.032% |
  | piatto 10–20 MeV, (0.2, 0.1 mm) | +0.495% | −0.060% | −0.066% |
  | piatto 10–20 MeV, (0.01, 1 µm) | −0.377% | −0.060% | +0.108% |

  Il range non dipende dalla definizione della luce: i due scan danno gli stessi valori, perché il
  trasporto è lo stesso evento per evento.
- **Diagnosi** (c), con le fluttuazioni di perdita spente in un binario di prova fuori dal repo:
  (0.2, 0.1 mm) +0.481%, (0.01, 1 µm) +0.068%. Gli effetti sono due, di segno opposto, entrambi
  da concavità di S/(1 + kB·S):
  - step lunghi: la media di S sullo step dà più luce dell'integrale (+);
  - step corti con perdite fluttuanti: la luce media per step cala (−).

  Birks per step con fluttuazioni non ha quindi un limite per step → 0.
- **Previsioni scritte prima dello scan, ed esito:**
  - P1 (luce monotona verso 0 al calare di dRoverRange, |·| < 0.05% al passo più fine):
    **fallita**, plateau +0.44% e cambio di segno con finalRange;
  - P2 (range insensibile al passo entro 2·10⁻⁴): **fallita**, si sposta fino a circa 0.25%,
    soprattutto sotto 2 MeV (+0.78% → +0.34% fra il passo più grosso e il più fine);
  - P3 (finalRange sotto 0.1 mm cambia la luce di meno di 0.05%): **fallita**, −0.83 punti;
  - P4 (eccesso a fine range di `ionIoni` invariato dalla config MuHad): **passata**; gli ioni
    hanno già finalRange 1 µm da `G4EmStandardPhysics_option4`;
  - P5 (costo circa 1/dRoverRange): **fallita**, 2.0–2.8 s per 10⁴ eventi su tutta la griglia,
    RSS < 92 MiB.
- **Decisione:**
  - per i **protoni** la luce di uno step è L(K_pre) − L(K_pre − E_dep), con
    L(K) = ∫₀^K dE / (1 + kB·S(E)) e S lo stopping elettronico medio delle tabelle di
    ionizzazione del protone in Geant4 (`G4EmCalculator::GetDEDX`), su una griglia logaritmica
    da 1 keV a 100 MeV con 50 punti per decade, interpolata log-log, con Simpson a 8 pannelli per
    intervallo (`sim/truth`, senza Geant4, con test). La luce di un protone che si ferma è L(E_p)
    qualunque siano gli step, che è la definizione con cui kB è misurato (D18: fit su stopping
    NIST); il deposito tiene le sue fluttuazioni;
  - un deposito oltre K_pre (stopping nucleare a fine range, deposito locale di un'interazione
    adronica) ha l'eccesso quenchato come le energie più basse della tabella; con kB = 0 la luce
    è esattamente il deposito;
  - le **altre particelle** restano su `G4EmSaturation`: per gli elettroni kB·S è piccolo, per
    gli ioni un solo kB non vale comunque (D18);
  - la funzione di passo resta **(0.2, 0.1 mm)**: la luce non ne dipende più e il costo del passo
    più fine non serve a nulla.
- **Residuo che resta** (−0.062%, uguale per tutte le config): differenza fra lo S di Geant4 e lo
  S totale di PSTAR del checker.
  - Sotto 2 MeV +0.13%: lo stopping nucleare (`G4NuclearStopping`, attivo sotto 1 MeV) deposita ma
    non è nella tabella del protone. Su PSTAR stesso L_elec/L_tot − 1 vale +0.17% a 1 MeV e
    +0.03% a 13 MeV (c, calcolo fuori dal repo).
  - Sopra 2 MeV circa −0.07%, quindi circa −0.1% al netto della parte nucleare: lo S di Bethe di
    Geant4 un po' sopra PSTAR (d), da non correggere; vale 1/8 della σ_Ep relativa attesa a
    10 MeV.
- **Alternative:**
  - Birks per step con passo più fine: escluso, non converge;
  - stopping totale nella tabella: richiederebbe `ComputeNuclearDEDX`, che spegne le
    fluttuazioni del modello condiviso e cambierebbe la simulazione;
  - luce calcolata in analisi da `truth_Ep`: perde la luce dei protoni che escono o interagiscono.
- **Conseguenze:**
  - stesso schema (v2) e stesse colonne: cambia il significato di `det_Evis` e `truth_p_Evis` per
    i run dal commit 3af4c3a in poi. I run precedenti con kB > 0 erano solo di prova;
  - la funzione di passo resta un parametro registrato, ma non è più una sistematica della luce dei
    protoni;
  - D18, "Conseguenze": la frase sulla dipendenza dal passo è superata da questa decisione.

## D22 — Rinculi nucleari del protone nello scattering Coulombiano (2026-11-20)

- **Problema:** nel pilota della fase 2 (`pilot2_mono_14MeV_4e5_8t`) un protone contenuto e senza
  interazioni adroniche ha `truth_p_Edep` − `truth_Ep` = −2.027 MeV (evento 246185, E_p =
  13.617 MeV). La previsione del bilancio (minimo ≥ −10⁻⁹ MeV) è fallita. Il run ripetuto con la
  traccia degli step dimostra la causa: uno step `CoulombScat` (scattering Coulombiano singolo,
  `G4eCoulombScatteringModel`, processo EM e non adronico) ha creato un rinculo di C12 da
  2.027 MeV. Sopra la soglia di rinculo il modello crea il nucleo bersaglio come secondario
  (`G4IonTable::GetIon`). Quell'energia lascia la traccia del protone, che per D19 non comprende i
  secondari, e finisce in `det_Edep` come deposito di uno ione.
- **Frequenza (d):** 1 caso su 137 806 protoni contenuti nei due piloti a 14 MeV; nessun altro
  caso fra −10⁻⁹ e −2 MeV.
- **Decisione:**
  - nuova colonna `truth_p_n_nuclear_recoil` (schema 3): il numero di step non adronici del
    protone seguito che hanno fra i secondari una particella con numero barionico > 0. Per un
    rinculo su H il modello crea un G4Proton, per un rinculo su C uno ione: hanno tutti e due
    numero barionico 1 o più. Elettroni e fotoni hanno 0. Uno step adronico è già contato in
    `truth_p_n_hadronic`, quindi i due contatori non contano mai lo stesso step;
  - le verifiche dei depositi (`checks_version` 8) escludono questi protoni da
    `contained_balance_MeV` e dal confronto PSTAR, e li contano a parte in
    `deposits.protons_nuclear_recoil`. Il loro range e la loro luce non sono quelli di un protone
    che rallenta solo per ionizzazione;
  - le verifiche dei depositi richiedono lo schema 3. Un file v2 darebbe 0 alla colonna nuova e
    rimetterebbe questi protoni nel bilancio. Gli unici run v2 sono i due piloti, che non sono run
    di validazione.
- **Alternative:**
  - sommare l'energia del rinculo al deposito del protone: cambierebbe la definizione di D19, e la
    luce dello ione segue `G4EmSaturation`, non la luce continua dei protoni (D21);
  - spegnere lo scattering Coulombiano singolo: cambierebbe la fisica di `G4EmStandardPhysics_option4`;
  - allargare la tolleranza del bilancio: nasconderebbe il caso invece di contarlo.
- **Conseguenze:**
  - `docs/validation/2_predictions.md` passa alla versione 2, corretta dopo il pilota e detta tale:
    l'oracolo del bilancio vale sui protoni senza rinculi nucleari, e il loro conteggio ha una
    previsione propria;
  - i piloti vanno ripetuti con lo schema 3 prima della produzione;
  - `det_Edep` e `det_Evis` non cambiano: il rinculo deposita comunque nello scintillatore.

## D23 — Lettura delle righe in ordine di `event_id` (2026-11-21)

- **Problema:** nei run di sensibilità a kB (`2_mono_14MeV_kB0p{145,155,165}_1e6`, verifiche `c1`)
  quattro verifiche di identità della §1 di `2_predictions.md` v2 sono fallite:
  `contained_balance_MeV.mean` e `pstar.range_mm` differiscono fra run che hanno gli stessi eventi.
  Tutti gli interi sono identici, i float differiscono di al più 1.3·10⁻¹² relativo. Causa: in MT
  Geant4 assegna gli eventi ai thread dinamicamente, quindi la divisione fra i file
  `events_t<k>.root` cambia da run a run. Il checker leggeva i file uno dopo l'altro e sommava i
  double in un ordine che dipendeva da quella divisione.
- **Decisione:**
  - il reader (`forEachEventRow`) fonde i file per thread in ordine crescente di `event_id`
    (`mergeByEventId`, k-way merge). Ogni worker elabora id crescenti, quindi ogni file è già
    ordinato; il merge rifiuta (`RunReadError`) un file che non sia strettamente crescente, perché
    l'ordine tornerebbe a dipendere dalla divisione;
  - un `event_id` presente in due file (run rotto) si visita due volte, nell'ordine dei file, così
    `row_count` di `loadRun` continua a rivelarlo;
  - `checks_version` passa a 9. Nuovo oracolo nello smoke test: lo stesso run a 1 e a 2 thread dà
    CSV identici byte per byte e `checks.json` identico a meno di `run_id`.
- **Alternative:**
  - somma compensata (Kahan) o in precisione estesa: riduce la differenza ma non la annulla, e
    l'identità resterebbe un'approssimazione;
  - tolleranza relativa sulle verifiche di identità: cambierebbe una soglia dopo il run, cosa
    esclusa; e nasconderebbe un'eventuale differenza vera;
  - un solo thread per i run di sensibilità: lento, e non risolve il problema per gli altri run.
- **Conseguenze:**
  - le soglie di `2_predictions.md` v2 non cambiano. Il fallimento delle verifiche `c1` resta scritto
    in `2_results.md`; le verifiche si ripetono come `c2` sugli stessi run, senza nuovi run di
    simulazione;
  - il merge legge tutti i file in parallelo (un cursore di riga per ntupla nel reader di Geant4);
    la memoria resta quella di una riga per file.

## D24 — Offset HP di E_p corretto in analisi (2026-11-21)

- **Problema:** il modello HP di Geant4 risolve lo scattering n–p con la massa dell'atomo di H-1
  (AWR = 0.999167 in G4NDL 4.7.1) e crea un G4Proton con l'impulso che ne risulta (D14). Il
  `truth_Ep` della simulazione è così più alto della cinematica fisica di (5.33–5.49)·10⁻⁴ in
  relativo su 10–20 MeV, circa 7.6 keV a E_p = 14 MeV. I protoni veri non hanno questa
  convenzione: lasciata nei dati, darebbe un bias relativo dello stesso ordine su E_n ricostruito.
- **Decisione:**
  - in `analysis/io` (fase 4) il `truth_Ep` dei candidati si converte all'E_p fisico con
    `physicalProtonEnergyFromHp` (`common/hp_kinematics.hpp`) **prima** di applicare σ_Ep. La
    funzione ricava cos_cm esatto con la massa HP e lo rimappa con la cinematica del protone
    nudo: è la stessa conversione dell'oracolo di D14, con residuo all'arrotondamento;
  - la massa del bersaglio è un argomento: si legge dai dati nucleari usati dalla run (oggi
    `hp_elastic_H1.awr` in `checks.json`), mai scritta a mano nell'analisi. Come arriva
    all'analisi si decide nella fase 4;
  - restano non corretti: θ_p, che il modello sposta al più di 1.4·10⁻⁶ rad (trascurabile
    rispetto a σ_θ), e `det_Edep` / `det_Evis`, con l'energia creata fino a circa 6·10⁻⁴ · E_n,
    che resta una sistematica dei depositi.
- **Alternative:**
  - ricalcolare `protonEnergy(truth_En, truth_theta_p)`: più semplice, ma lascia fino a
    5.7·10⁻⁶ relativo per lo spostamento di θ_p, e costa quanto la conversione esatta;
  - modellare l'offset nella verosimiglianza: porterebbe nella statistica una convenzione di
    Geant4 che i dati veri non hanno;
  - riportarlo come sistematica: sprecherebbe un effetto noto all'arrotondamento.
- **Conseguenze:**
  - le funzioni di cinematica HP passano da `sim/validation` a `common`, header-only, perché
    l'analisi dipende solo da `common`; i loro test girano anche in `ci`;
  - i test fissano: andata e ritorno entro 10⁻¹² MeV, identità con M_T = m_p, offset tolto in
    (5.33–5.49)·10⁻⁴, scarto da `protonEnergy(En, θ_p)` sotto 6·10⁻⁶ relativo.

## D25 — Sorgente fino a 19.5 MeV nelle fasi 5–6 (2026-11-21)

- **Problema:** in 19.5–20 MeV lo stato finale elastico n–p è HP o CHIPS secondo una rampa (D16).
  Lì la distribuzione angolare non è più solo quella di JEFF-3.3, quindi il motore statistico
  vedrebbe un modello diverso sull'ultimo mezzo MeV.
- **Decisione:**
  - i config dei run delle fasi 5–6 hanno `En_max_MeV` ≤ 19.5: lo stato finale elastico è solo
    HP/JEFF-3.3;
  - il dominio del prior su E_n si decide nella fase 5;
  - nella fase 9 il bordo 19.5–20 MeV torna come sistematica, con la verifica `model_mix` già
    validata nella fase 2.
- **Alternative:**
  - tenere 20 MeV e modellare la rampa: lega l'analisi a un dettaglio di Geant4 prima che il
    Caso A sia validato;
  - trattare subito il bordo come sistematica: mescola un effetto noto con la validazione del
    motore.
- **Conseguenze:** i config della fase 2 restano invariati; i nuovi config delle fasi 5–6
  scrivono il limite esplicitamente.

## D26 — Eccesso di `ionIoni` sui rinculi C12 come sistematica documentata (2026-11-21)

- **Fatti:** a fine percorso `ionIoni` deposita più del K iniziale dello step su 56 rinculi C12
  su 486 e 4 Be9, con K fra 0.11 e 0.40 MeV; in totale 0.74 MeV su 492 MeV di rinculi C12, cioè
  0.15% del loro deposito (fase 2). La causa non è verificata; non dipende dalla funzione di
  passo del protone (D21). Sulla luce pesa meno, per il quenching dei rinculi pesanti.
- **Decisione:** nessuna modifica al codice. L'eccesso resta fra le sistematiche del fondo da
  carbonio; la causa si studia solo se nella fase 9 quel fondo conterà.
- **Alternative:**
  - studiare ora la causa nel codice di Geant4: costo alto per un effetto dello 0.15% su eventi
    n–C, che la selezione dei candidati a singolo n–p esclude;
  - correggere il deposito a mano: senza la causa sarebbe una correzione inventata.
- **Conseguenze:** gli eventi n–C non entrano nelle fasi 5–8; nella fase 9 il numero va riportato
  accanto a quello del C-13 (D15).

## D27 — Blocco 60 × 20 × 20 mm e sorgente a campo lontano (2026-11-21)

- **Problema:**
  - il cubo di 10 cm delle fasi 1–2 era provvisorio;
  - il fascio a matita illumina una sola corda del blocco. Il Caso A e i Casi B–C vogliono
    invece un neutrone da sorgente lontana: Ω_n è la stessa in tutto il blocco e l'interazione
    può avvenire in tutto il volume.
- **Decisione — geometria:**
  - blocco di BC-408 (materiale e kB di D18) di 60 × 20 × 20 mm, centrato nell'origine, con il
    lato lungo su X: è l'unico volume 3D del progetto di ottimizzazione, che non definisce
    materiale né sorgente;
  - il mondo resta nel vuoto (D10), senza supporti, perché qui interessa la statistica del
    singolo scattering.
- **Decisione — sorgente a campo lontano (`source.beam = "far_field"`):**
  - per ogni evento la direzione Ω_n è fissa da config oppure isotropa sulla sfera;
  - il punto di partenza è uniforme su un disco ⟂ Ω_n, centrato sulla retta per il centro del
    blocco, a distanza R_d + margine dal centro all'indietro lungo −Ω_n, dove R_d = √(30² + 10² +
    10²) mm = 33.17 mm è la mezza diagonale. Nessun neutrone parte dentro il blocco;
  - raggio del disco:
    - con direzione fissa: la massima distanza da Ω_n dei vertici del blocco, cioè il cerchio
      minimo centrato che contiene l'ombra del blocco (14.14 mm lungo +X);
    - con direzione isotropa: R_d per tutte le direzioni. Un raggio per direzione cambierebbe
      il peso relativo delle direzioni fra gli eventi che colpiscono il blocco, che deve
      restare proporzionale all'area proiettata, come per un flusso isotropo vero;
  - è un fascio parallelo esteso: illuminazione uniforme del blocco per ogni Ω_n e ipotesi di
    campo lontano esatta. La distanza fisica dalla sorgente non entra;
  - accettanza geometrica attesa = area proiettata / area del disco, con area proiettata
    4·(|Ω_x| h_y h_z + |Ω_y| h_x h_z + |Ω_z| h_x h_y). Vale 400 / 628.3 = 63.7% lungo +X e
    ⟨S/4⟩ / πR_d² = 1400 / 3455.8 = 40.5% isotropa. Gli eventi che mancano il blocco costano
    solo il trasporto nel vuoto;
  - nei raggi che colpiscono il blocco con direzione isotropa la corda media è 4V/S =
    17.14 mm (formula di Cauchy): un oracolo esatto del campionamento.
- **Decisione — fascio a matita (`source.beam = "pencil"`):** resta, con `position_mm` e
  `direction`, per la regressione con le fasi 1–2. I config v2 si leggono come pencil.
- **Alternative:**
  - punto unico con direzioni isotrope: campo vicino, cioè Ω_n cambia lungo il blocco e la
    maggior parte delle direzioni manca il blocco;
  - fascio a matita: non copre il volume, quindi non misura l'effetto dei bordi sui protoni che
    escono;
  - disco di raggio minimo anche per la direzione isotropa: più efficiente, ma distorce la
    distribuzione delle direzioni che interagiscono;
  - interazione forzata con pesi: più efficiente, ma porta pesi in un'analisi che oggi conta
    eventi.
- **Conseguenze:**
  - config v3 (`source.beam`, `margin_mm` obbligatorio per `far_field`);
  - schema v4 con il punto di partenza del neutrone (`truth_n_start_{x,y,z}`), perché le
    verifiche calcolano la corda per evento;
  - le verifiche della fase 2 che assumevano una sola corda passano a quantità per evento
    (`checks_version` 10). Con la corda costante del fascio a matita devono dare gli stessi
    numeri.

## D28 — Massa del bersaglio HP registrata in `run.json` (2026-11-22)

- **Problema:** D24 converte `truth_Ep` all'E_p fisico con la massa HP di H-1, che va letta dai
  dati della run e mai scritta a mano. Oggi la massa sta solo in `checks.json`
  (`hp_elastic_H1`), cioè nell'esito di una verifica, non nella provenienza del run.
- **Decisione:**
  - `run.json` passa a `run_record_version` 2. Aggiunge il blocco `hp_elastic_H1`, scritto da
    `riptide_sim` con le stesse funzioni della verifica (`readHpElasticH1`, `hpJson`): file,
    sha256, `rep_flag`, `awr`, `frame_flag`, `target_mass_MeV`. È lo stesso blocco di
    `checks.json`;
  - le colonne non cambiano, quindi `kSchemaVersion` resta 4;
  - l'analisi usa `target_mass_MeV` e rifiuta il run se non vale `awr` · m_n entro 10⁻¹² in
    relativo;
  - per i run con `run_record_version` 1 (quelli della fase 3) la massa si prende solo dal
    `checks.json` di una verifica dello stesso run, indicata in modo esplicito
    (`--hp-from-check <check_id>`), con lo stesso `run_id`. Senza nessuna delle due fonti il
    run è rifiutato. Non esiste un valore di default.
- **Alternative:**
  - leggere G4NDL al momento dell'analisi: lega l'analisi all'installazione di Geant4 e può
    leggere una libreria diversa da quella usata dal run;
  - leggere sempre `checks.json`: rende l'analisi dipendente da una verifica che non è
    obbligatoria.
- **Conseguenze:** `readRunRecord` di `sim/reader` accetta le versioni 1 e 2; il test
  d'integrazione controlla `hp_elastic_H1.awr` = 0.999167 in `run.json`.

## D29 — Struttura di `analysis/io` (2026-11-22)

- **Problema:** `sim_reader` legge i file con il lettore ROOT di Geant4, ma l'analisi può
  dipendere solo da `common` (e da ROOT). La CI non ha né Geant4 né ROOT, ma le parti pure
  dell'I/O devono essere testate anche lì.
- **Decisione:**
  - `riptide::io_core`, senza ROOT, compilata anche in `ci`: lettura di `run.json`, selezione e
    cut-flow, conversione HP, risposta del rivelatore, scrittura di `analysis.json`;
  - `riptide::io_root`: solo il lettore `RDataFrame`, che lega le colonne per nome dallo schema
    di `common/event_schema.hpp`;
  - il buffer struct-of-arrays `EventSample` sta in `common/event_sample.hpp` e il
    campionatore gaussiano in `common/random.hpp`, così `analysis/stats` li usa senza dipendere
    da `analysis/io`;
  - il campionatore gaussiano è una Box–Muller propria su `std::mt19937_64` con uniforme a 53
    bit: lo standard fissa la sequenza del motore, non quella di `std::normal_distribution`,
    che cambia da una libreria standard all'altra;
  - il lettore rifiuta una colonna richiesta che manca nello schema del run, invece di metterla
    a 0 come `sim_reader`: nell'analisi uno 0 silenzioso sarebbe un bias;
  - il lettore riordina per `event_id` e rifiuta id ripetuti (come D23), così lo stesso run a 1 e
    a 2 thread dà lo stesso campione bit per bit.
- **Alternative:** riusare `sim_reader` (porta Geant4 nell'analisi); spostare il lettore in
  `common` (porta ROOT a valle di `analysis/stats`).

## D30 — Selezione dei candidati a singolo n–p (2026-11-22)

- **Decisione:**
  - la selezione è una catena di tagli con nome, e ognuno lascia il suo conteggio nel cut-flow:
    - `single_np`: `truth_event_class == SingleNp`;
    - `contained`: `truth_p_exits == 0`;
    - `evis_min`: `det_Evis` sopra una soglia data come parametro, senza default;
  - il default per la fase 5 è `single_np` + `contained`: il port del Caso A si confronta con il
    toy sulla sola cinematica.
- **Motivo:** gli eventi a più urti, i protoni che escono (fino al 17% sopra 15.6 MeV nella
  fase 3) e la soglia sono effetti del rivelatore da misurare dopo come sistematiche, non da
  mescolare al primo confronto con il toy.
- **Conseguenze:** il pull aperto della fase 3 sulla frazione n–C non elastica non pesa qui,
  perché la selezione usa la classe vera e non le frazioni. Resta da indagare prima della
  fase 6.

## D31 — Risposta del rivelatore in analisi (2026-11-22)

- **Decisione:**
  - prima la conversione HP, poi il rumore: `Ep_phys = physicalProtonEnergyFromHp(truth_Ep,
    truth_En, M_T)`, `reco_Ep = Ep_phys + N(0, σ_Ep)`;
  - Caso A: `reco_theta_p = truth_theta_p + N(0, σ_θ)`, lo stesso modello 1D della sua
    verosimiglianza;
  - Casi B e C: `reco_p_dir` è la traccia vera spostata nel piano tangente con due componenti
    N(0, σ_θ), come nel toy. Le estrazioni sono indipendenti da quelle di θ_p;
  - il motore è per evento, inizializzato da (seme dell'analisi, `event_id`) attraverso
    SplitMix64: cambiare un taglio non cambia il rumore degli altri eventi, e l'ordine di
    lettura non conta. Le funzioni di smearing ricevono il motore come argomento;
  - `analysis.json` registra seme, σ_Ep, σ_θ, tagli, fonte e valore di M_T, schema letto e
    `git_describe`.
- **Alternative:** un solo modello 3D con θ̂_p ricavato da `reco_p_dir`. È più realistico, ma
  vicino a θ_p ≈ 0 θ̂_p non è più gaussiano e la verosimiglianza del Caso A andrebbe rivista già
  nella fase 5. Resta il passo successivo, da misurare come sistematica.
- **Aggiunta (2026-11-22, dopo i run della fase 4):**
  - un evento `single_np` può non avere il protone (`truth_valid_np` = 0). In quel caso i file
    scrivono `truth_p_dir` = (0, 0, 0). Nella prima versione lo smearing normalizzava quel
    vettore nullo e dava `reco_p_dir` = NaN. Il NaN rendeva `null` il pull della direzione, e
    `std::max` lo nascondeva nei massimi;
  - ora un vettore nullo significa "nessuna traccia": `reco_p_dir` = (0, 0, 0), senza estrazioni
    di direzione. Le estrazioni di Ep e θ_p restano invariate;
  - una direzione vera né unitaria né nulla (entro 10⁻⁹) è un errore. Un residuo non finito nei
    riassunti è un errore;
  - i riassunti escludono questi eventi dal solo blocco `p_dir` e li contano in
    `events_without_track`, per cui `analysis_record_version` passa a 2;
  - nel campione restano eventi con `truth_Ep_phys` = 0. Il loro peso nella verosimiglianza del
    Caso A va deciso nella fase 5.

## D32 — Cinematica e griglie del Caso A (2026-11-22)

- **Decisione:**
  - la verosimiglianza usa `riptide::kinematics::protonEnergy`, la cinematica esatta con le masse
    di n e p diverse, la stessa che dà `truth_Ep_phys` nella fase 4;
  - l'inversa si aggiunge in `common/np_kinematics.hpp`: con X = p²cos²θ_p ed S = E_tot + m_t,
    dalla formula di `recoilEnergy` si ha X = Ep·S²/(2m_t + Ep), quindi cos²θ_p = X/p². Per masse
    uguali si riduce a quella del toy, Ep(2m + En)/(En(2m + Ep));
  - griglia di E_n: [E_min, E_max], parametro, default [10, 20] MeV con 500 punti. Il dominio del
    prior su E_n lasciato aperto da D25 è quindi deciso: default come il toy, con la sorgente
    della fase 6 (E_n ≤ 19.5 MeV) tutta all'interno;
  - griglia di θ_p: [0, π/2] con 500 punti;
  - griglie contigue, in doppia precisione, immutabili dopo la costruzione. Niente float32.
- **Motivo:** la cinematica del toy usa masse uguali, con uno scarto relativo ~5·10⁻⁷, ma la
  verosimiglianza deve essere coerente con la verità della fase 4. Il float32 del toy era una
  scelta di velocità di numpy; in C++ la doppia precisione non costa abbastanza da giustificarlo.
- **Alternative:** masse uguali come il toy (scarto trascurabile, ma una seconda cinematica da
  mantenere); dominio fisso [10, 20] (un dominio largo serve alla prova di robustezza al prior).
- **Conseguenze:** il confronto con il toy è entro la tolleranza statistica, non bit per bit.

## D33 — Modello e API del Caso A (2026-11-22)

- **Decisione:**
  - `CaseAModel` riceve griglie, σ_Ep, σ_θ e il log-prior su θ_p, e precalcola una volta la
    tabella Ep(E_n,i, θ_j);
  - `logLikelihood(Êp, θ̂)` restituisce un vettore sulla griglia di E_n:
    `log Σ_j exp[−½((Êp − Ep(E_n, θ_j))/σ_Ep)² − ½((θ̂ − θ_j)/σ_θ)² + log π(θ_j)]`, con un
    logsumexp stabile. Le costanti comuni a tutte le E_n sono omesse, come nel toy;
  - `singleEventLogPosterior` = loglik + log prior su E_n;
  - `combineLogLikelihoods(loglik_k, log_prior)` = log_prior + Σ_k loglik_k: il prior si conta una
    volta sola;
  - i prior sono funzioni con nome passate come argomento: `flatEnergyPrior` su E_n,
    `physicalRecoilAnglePrior` (massa esatta della cella, ½(cos 2lo − cos 2hi), da sin 2θ_p) e
    `flatRecoilAnglePrior` su θ_p. Nessun ramo sul tipo di prior;
  - ingresso: due span (`reco_Ep`, `reco_theta_p`), presi da un `EventSample`. Un valore non
    finito o una σ ≤ 0 sono errori;
  - nessun ROOT e nessun Geant4 in `analysis/stats`.
- **Motivo:** la tabella non dipende dagli eventi, quindi il costo per evento è un solo
  passaggio sulla griglia 500×500. Il prior come argomento evita gli `if` sul caso e permette la
  prova di robustezza senza toccare il modello.
- **Alternative:** ricalcolare Ep per ogni evento come il toy (più lento, stessi numeri).

## D34 — Taglio `has_proton` (2026-11-22, aggiunta a D30)

- **Decisione:**
  - un nuovo taglio con nome, `has_proton`: `truth_valid_np == 1`. Sta sulla verità, come
    `single_np`, e lascia il suo conteggio nel cut-flow;
  - il default del Caso A diventa `single_np,has_proton,contained`;
  - lo schema non cambia: `truth_valid_np` è già una colonna letta.
- **Motivo:** un evento `single_np` senza protone tracciato ha `truth_Ep_phys` = 0 e un `reco_Ep`
  fatto di solo rumore (fase 4, §2–§3: 89 / 17 / 19 eventi). Tenerlo nella verosimiglianza del
  Caso A darebbe un evento a Ep ≈ 0, cioè una misura che il modello non sa descrivere.
- **Alternative:** tenerli come parte della risposta del rivelatore (servirebbe un modello della
  loro distribuzione); escluderli in silenzio nella statistica (il conteggio sparirebbe dal
  cut-flow).
- **Conseguenze:** il motore per evento (D31) lascia identiche bit per bit le estrazioni degli
  altri eventi. La frazione persa (~10⁻⁴) si annota come sistematica per la fase 9.

## D35 — Pseudo-esperienze e prestazioni del Caso A (2026-11-22)

- **Decisione:**
  - `analysis/stats` contiene un generatore del modello del toy, che serve solo a validare il
    port: E_n uniforme sul dominio; cos θ_CM uniforme in [−1, 1], poi `protonEnergyFromCosCm`, poi
    θ_p dall'inversa di D32; infine Êp e θ̂ con il rumore gaussiano 1D di D31. Il motore di ogni
    evento viene da `perEventEngine(seme, k)`;
  - target: ≤ 0.10 s per la verosimiglianza di 1000 eventi, a thread singolo, nella build
    `release`. Prima un'implementazione semplice; poi profilo (`perf`) e ottimizzazioni esatte
    (layout, precalcolo, riuso dei buffer), rilanciando gli oracoli dopo ogni passo;
  - i tempi vanno in `docs/validation/5_results.md`, non in un test con soglia;
  - il parallelismo sugli eventi entra solo se il profilo lo giustifica, con una decisione a
    parte. Una finestra su θ_p attorno a θ̂ non è ammessa senza un limite d'errore dimostrato.
- **Motivo:** la calibrazione con verità dal prior è la prova più diretta che il port è corretto:
  i suoi numeri (pull, coverage, contrazione) si confrontano con quelli del toy.
- **Alternative:** generare con Geant4 (è la fase 6; qui mescolerebbe port e rivelatore).

### Addendum a D35 — Finestra su θ_p con limite d'errore (2026-11-22)

- **Finestra:** per ogni riga E_n,i e ogni evento la somma su θ_j si limita alle celle che possono
  contare.
  - **Limite inferiore:** sia L il valore più grande fra i termini delle celle più vicine a Êp
    lungo la riga e a θ̂. Allora L ≤ max_j t_ij.
  - **Celle tenute:** sia a_max = max_j log π(θ_j), e sia C = 50. Un termine
    t_ij = −½z_E² − ½z_θ² + log π(θ_j) può essere ≥ L − C solo se |z_E| ≤ r e |z_θ| ≤ r, con
    r = √(2(a_max + C − L)). Le celle fuori da queste due bande hanno t_ij < L − C ≤ max − C.
  - **Ricerca:** nella tabella Ep(E_n,i, θ_j) ogni riga è non crescente in θ_j. Il costruttore
    lo verifica e, se non vale, la somma si fa sull'intera riga. Le due bande sono quindi
    intervalli contigui, trovati per bisezione.
  - **Errore:** la massa tralasciata, relativa a quella tenuta, è ≤ n_θ·e^(−C) ≈ 10⁻¹⁹ con
    n_θ = 500. È anche il limite dell'errore assoluto sulla log-verosimiglianza, ed è sotto
    l'arrotondamento di un double.
  - **Ripiego:** se L non è finito (prior −inf su quelle celle) si somma l'intera riga.
- **Altre ottimizzazioni esatte:** precalcolo di 1/σ; logsumexp senza ricontrollo dentro il
  ciclo, perché gli ingressi sono validati prima; buffer riusati fra gli eventi.
- **Profilo:** `perf` non è disponibile su questa macchina, perché `perf_event_paranoid` vale 4.
  Le ottimizzazioni vengono dalla struttura del ciclo, misurando ogni passo con
  `riptide_stats_bench`. Il costo è dominato da `exp`, circa 3.7 ns a chiamata.
- **Esito:** 0.17 s per 1000 eventi, contro 1.18 s dell'implementazione semplice, quindi il target
  di 0.10 s non è raggiunto. I numeri sono in `docs/validation/5_results.md`. Con C = 50 la
  finestra tiene in media circa 95 celle per riga, cioè ±10σ_Ep lungo la cinematica. Scendere
  sotto 0.10 s richiederebbe un `exp` vettoriale, che non è esatto, oppure il parallelismo sugli
  eventi. Il parallelismo resta, come sopra, una decisione a parte.
- **Verifica:** un test Catch2 confronta la finestra con la somma su tutte le celle, con due prior,
  tre risoluzioni ed eventi ai bordi (Êp < 0, Êp sopra ogni Ep_max, θ̂ fuori da [0, π/2]). Lo
  scarto relativo è < 10⁻¹³. Sulla matrice del benchmark (1000 × 500) la finestra coincide bit per
  bit con la somma completa.

## D36 — Campioni Geant4 della fase 6 (2026-11-23)

- **Decisione:**
  - i nuovi config v3 sono copie di `config/3_ff_mono_14MeV_x.json`: stesso blocco (D27), stessa
    fisica, stesso kB (D18) e stessa funzione di passo. Cambia solo lo spettro;
  - `config/6_ff_flat_10_19p5MeV_x.json`: spettro piatto 10–19.5 MeV (D25), direzione +X, 10⁷
    eventi. Serve alla checklist per evento (bias, risoluzione, pull, coverage) con la verità
    distribuita sul dominio;
  - `config/6_ff_mono_{11,17,19}MeV_x.json`: 10⁶ eventi ciascuno. Con
    `3_ff_mono_14MeV_x_1e7`, che esiste già e ha le verifiche `c1` della fase 3, danno quattro
    energie per la contrazione 1/√N e per la robustezza al prior;
  - seme `20260907` per tutti i nuovi run, `release`, 8 thread, `scripts/run_capped.sh -m 5G`.
    Su ogni nuovo run si eseguono le verifiche `riptide_sim_check` `c1` prima dell'analisi.
- **Motivo:**
  - la direzione fissa +X è la configurazione del Caso A: Ω_n è nota e uguale in tutto il blocco;
  - 11 e 19 MeV stanno vicino ai bordi del prior, dove il toy della fase 5 vede il bias più
    grande; 14 e 17 MeV stanno nel mezzo;
  - lo stesso seme del progetto basta perché i run hanno config diversi. Il 14 MeV +X ha già un
    gemello (seme 20260908), che resta disponibile come replica.
- **Alternative:**
  - il piatto isotropo della fase 3: Ω_n cambia per evento, che il Caso A non modella;
  - mono a 10⁷ per tutte le energie: la contrazione a N = 1000 vuole al più qualche centinaio di
    esperienze, che 10⁶ eventi dovrebbero dare. Il numero di eventi selezionati si stima nelle
    previsioni dal pilota;
  - spingere la sorgente a 20 MeV: escluso da D25.
- **Conseguenze:** tre file in `config/` e le loro righe in `config/README.md`. Nessuna modifica
  a `sim/`.

## D37 — Pseudo-esperienze da campioni Geant4 (2026-11-23)

- **Decisione — dal run al campione:**
  - lettura, tagli, conversione HP e risoluzione sono quelli di `riptide_analysis_prepare`:
    catena `single_np,has_proton,contained` (D30, D34), `Ep_phys` da D24, smearing D31 con
    σ_Ep = 0.10 MeV, σ_θ = 0.08 rad e motore per evento;
  - le colonne dell'`EventSample` selezionato riempiono uno `stats::PseudoSample`:
    `truth_En`, `truth_Ep_phys` (come `truth_Ep`), `truth_theta_p`, `reco_Ep` e `reco_theta_p`.
    Il motore statistico non sa se gli eventi vengono dal toy o da Geant4.
- **Decisione — esperienze:**
  - una funzione pura in `analysis/stats` divide gli indici 0…n−1 in M blocchi disgiunti di N
    indici, dopo una permutazione con l'RNG passato come argomento e il seme registrato. Gli
    indici avanzati non si usano;
  - disgiunti vuol dire esperienze indipendenti, quindi coverage con errore binomiale onesto.
    Nessun evento entra in due esperienze della stessa riga;
  - righe come in `riptide_stats_calibrate`, per un confronto voce per voce:
    - checklist per evento sul piatto, con un sottoinsieme casuale di 20 000 eventi selezionati,
      come la fase 5, e la curva del bias in 5 bin;
    - contrazione con N = 10, 100, 1000 sui quattro mono;
    - prior su E_n ([10, 20] contro [10, 30]) con N = 1, 3, 10, 30, 100;
    - prior su θ_p (sin 2θ_p contro piatto) a 14 MeV con N = 1000;
  - M per riga è quello della fase 5, ridotto se gli eventi selezionati non bastano. La
    riduzione si scrive nelle previsioni, mai dopo il run;
  - si riusano senza modifiche `makeSummaries`, `appendSingleEventSummaries`,
    `appendCombinedSummary`, `summariseEnsemble`, `priorShift` (`calibration.hpp`) e `biasCurve`
    (`validate.hpp`). Prior e griglie come D32 e D33: E_n piatto su [10, 20] con 500 punti,
    sin 2θ_p con 500 punti.
- **Decisione — app:** `riptide_case_a_geant4`, in `analysis/apps`, una chiamata per run.
  - Scrive `case_a.json` in `statsDir(<radice>, <id>)`. Il JSON ha `run_id`, seme
    dell'analisi, seme dello split, σ, catena dei tagli con il cut-flow, schema letto e
    `git_describe`;
  - le sezioni hanno la stessa struttura di `calibration.json`;
  - `--pilot` divide i conteggi per 10 e lascia N invariato, come la fase 5;
  - come per la calibrazione, il file non contiene tempi: due esecuzioni uguali danno gli
    stessi byte.
- **Motivo:**
  - la fase 5 ha validato il motore con quelle funzioni e quel protocollo. Cambiare solo la
    sorgente degli eventi isola l'effetto del rivelatore;
  - un indice permutato e poi tagliato in blocchi è la via più semplice a esperienze disgiunte,
    ed è verificabile con un test (unione = sottoinsieme, intersezioni vuote, riproducibilità).
- **Alternative:**
  - esperienze con rimpiazzo (bootstrap): eventi condivisi, quindi esperienze correlate e
    coverage con errore sottostimato;
  - scrivere su disco il campione selezionato e smussato: un formato in più da mantenere, mentre
    la selezione di 10⁷ eventi sta in memoria (fase 4).
- **Conseguenze:**
  - la funzione di split entra in `analysis/stats` con il suo test Catch2 e compila anche in
    `ci`;
  - l'app dipende da `io_root` e compila solo con ROOT. Lo smoke test usa un run piccolo ed è
    saltato senza SSD, come gli altri test di integrazione.

## D38 — Confronto toy vs Geant4 (2026-11-23)

- **Decisione:**
  - `riptide_case_a_geant4` ha una sorgente alternativa, `--toy`. Con la stessa catena di
    esperienze, gli eventi vengono dal generatore D35 con la stessa verità del run:
    - per il piatto, E_n uniforme su [10, 19.5];
    - per un mono, la sua E_n;
    - stessi n, N, M, σ e livelli (0.68, 0.90, 0.95), con seme registrato;
  - per ogni voce della checklist, la differenza Geant4 − toy, con l'errore combinato, è il costo
    del realismo del rivelatore su quella voce. Il toy è "circolare" su una voce se lì passa e
    Geant4 fallisce: la sua validazione dipendeva dall'aver generato con il modello
    dell'inferenza;
  - le soglie si fissano nelle previsioni, con Bonferroni sul numero di voci, come nella fase 5.
- **Motivo:**
  - la fase 5 confronta il port col toy su E_n uniforme in [10, 20] e in [11, 19]. Qui servono
    la stessa E_n dei mono e il limite di 19.5 MeV del piatto, quindi `calibration.json` non è
    un riferimento diretto;
  - generare il toy nella stessa app usa lo stesso codice di esperienze e di riassunto: la
    differenza viene solo dagli eventi;
  - è anche il controllo di coerenza dell'app: con `--toy` sul dominio [10, 20] deve riprodurre
    `validation/5_results.md` entro l'errore statistico;
  - la scelta della sorgente degli eventi non è un ramo sul posterior: il posterior del Caso A
    resta uno solo.
- **Alternative:**
  - estendere `riptide_stats_calibrate` con i mono: due app con lo stesso protocollo da tenere
    allineate;
  - un confronto a mano dei due JSON fuori dal repo: non riproducibile.
- **Conseguenze:** ogni run ha due `case_a.json`, uno Geant4 e uno toy. Il confronto è uno
  script di valutazione delle soglie, come nella fase 5, i cui numeri vanno in
  `docs/validation/6_results.md`.

## D39 — Effetti attesi da Geant4 che il toy non ha (2026-11-23)

Sono ipotesi (d). Le previsioni della fase 6 le traducono in verso e grandezza attesi, prima dei
run.

- **Taglio `contained`.**
  - Fatto: la perdita di protoni dai bordi dipende da E_p e dalla direzione, fino al 17% sopra
    15.6 MeV (fase 3, §5).
  - Ipotesi: i protoni ad alto E_p e in avanti (θ_p piccolo) vengono tolti di più. La
    distribuzione selezionata in (E_p, θ_p) non segue più cos θ_CM uniforme, mentre il prior
    sin 2θ_p e la verosimiglianza la assumono.
  - Atteso: un bias di E_n verso il basso, che cresce con E_n. Pull e coverage peggiorano
    soprattutto a N grande, dove il bias pesa più della larghezza.
- **σ_np(E) non costante.**
  - Ipotesi: nel piatto, le E_n selezionate non sono uniformi, perché la probabilità di
    interazione e di contenimento dipende da E_n. La "verità dal prior" vale solo in parte.
  - Conseguenza: la checklist sul piatto si riporta anche per bin di E_n (la curva del bias), e
    la coverage globale si legge con questo limite.
  - Il modello dello spettro selezionato è della fase 9.
- **Anisotropia in CM.**
  - Fatto: i dati G4NDL lungo +X hanno a₁ = −1.425·10⁻² ± 5.7·10⁻⁴ (fase 3, (c)), mentre il
    modello ha cos θ_CM uniforme.
  - Ipotesi: è un effetto piccolo rispetto al taglio `contained`. Si stima nelle previsioni dal
    verso di a₁.
- **Dominio del prior oltre la sorgente.**
  - Fatto: il prior arriva a 20 MeV, la sorgente si ferma a 19.5 (D25, D32).
  - Ipotesi: è lo stesso effetto per toy e Geant4, perché il toy con `--toy` usa la stessa
    verità. Quindi non entra nella differenza (D38), ma sposta entrambi i bias al bordo alto.
- **Fuori dal Caso A.** L'eccesso residuo di n–C non elastico lungo +X (d, fase 3) non conta:
  `single_np` toglie gli eventi n–C. La frazione di eventi senza protone (~10⁻⁴, D34) è tolta
  da `has_proton`.

## D40 — Accettanza in θ_p nel modello del Caso A (2026-11-23)

- **Contesto.**
  - Fase 6: con Geant4 il posterior è ~10% più largo del toy, e a 14 MeV, N = 1000, c'è un bias di
    −0.023 MeV.
  - La diagnosi a 14 MeV mostra che la causa è la distribuzione in θ_p degli eventi selezionati.
    Rispetto al modello isotropo, il bin 0–5° vale 0.58, i θ intermedi 1.06, gli 85–90° 0.91. Il
    toy ripesato a quella distribuzione riproduce Geant4. La cinematica è coerente: lo scarto in
    E_p è ~10⁻⁵ MeV.
  - Il cut-flow per θ_p (a) attribuisce il deficit:
    - in avanti, a `single_np`: n–p multiplo e misto toccano il 55% a 0–5°;
    - in misura minore, a `contained`: al più il 5.6%, a 20–25°;
    - all'anisotropia del primo urto: ±6%.
  - La forma dipende da E_n in modo non liscio, χ² = 420/306 per un modello lineare in E_n.
- **Decisione.**
  - **Verosimiglianza per evento:**
    `L_i(E_n) = Σ_j q(θ_j | E_n, sel) · G(Ep_obs; E_p(E_n, θ_j), σ_Ep) · G(θ_obs; θ_j, σ_θ)`,
    dove q è la probabilità della cella θ_j fra gli eventi selezionati a quella E_n,
    normalizzata per riga. Prende il posto di `log_prior_theta`: una tabella n_En × n_θ invece di
    un vettore. q assorbe insieme anisotropia, selezione e cinematica relativistica.
  - **Niente A(E_n).** La normalizzazione globale della selezione resta fuori: la verosimiglianza
    è condizionata alla selezione, e il Caso A non usa il tasso. Lo spettro selezionato resta
    della fase 9 (D39).
  - **Stima di q.**
    - **Sorgente.** Il run piatto indipendente `6_ff_flat_10_19p5MeV_x_1e7_s20260908`, separato
      dai campioni su cui si valuta, con la catena `single_np,has_proton,contained` sulla verità
      (`truth_En`, `truth_theta_p`).
    - **Bin.** In E_n, bin da 1 MeV sul dominio della sorgente; l'ultimo è 19–19.5. In θ_p, 18 bin
      da 5° su [0, π/2]. A 0–5° sono ~140 eventi per bin, cioè l'8.5%, e l'impatto si misura
      (Conseguenze).
    - **Rapporto.** Per bin si stima r_kt = (n_kt / N_k) / m_rel(t; E_k): la frazione del bin di
      θ diviso la massa isotropa in CM con cinematica relativistica (E_p uniforme), al centro E_k
      del bin di E_n. È più liscio della densità, e vale 1 per un rivelatore ideale isotropo.
    - **Interpolazione.** r si interpola bilinearmente fra i centri dei bin (E_k, θ_t) e si tiene
      costante oltre il primo e l'ultimo centro, in E_n come in θ. Fra 19.25 e 20 MeV, e sul
      dominio largo fino a 30 MeV, vale quindi l'ultima riga: la sorgente si ferma a 19.5 (D25).
    - **Celle del modello.** Massa della cella j alla riga E_i:
      `∫_cella r(θ, E_i) dP_rel(θ | E_i)`, con il punto medio su sotto-intervalli uguali in
      u = E_p / E_p,max, poi normalizzata per riga. Con r ≡ 1 dà la massa relativistica esatta.
  - **Toy gemello (D38).** Con l'accettanza, il toy estrae l'evento come oggi e lo accetta con
    probabilità r(θ_p, E_n) / max r. Gli eventi provati usano motori per evento, con il seme
    registrato. Il toy resta il gemello «modello perfetto» del modello con q.
  - **Default invariato.** Senza `--acceptance <stats_id>` il modello resta quello attuale, con
    il prior non relativistico sin 2θ_p, e i record della fase 6 restano riproducibili bit per
    bit. Con `--acceptance`:
    - il modello «fisico» e quello del dominio largo usano q;
    - il confronto col prior piatto in θ resta com'è.
- **Motivo.**
  - Il bias e la larghezza in eccesso sono un difetto del modello, non un costo inevitabile del
    rivelatore. L'accettanza è misurabile con la simulazione e indipendente da E_n vera.
  - Condizionare sulla selezione evita di modellare il tasso. Un run indipendente evita di
    stimare q sugli stessi eventi su cui si valuta.
  - Il rapporto rispetto alla massa relativistica separa la fisica nota (cinematica) da quella
    misurata (selezione, anisotropia).
- **Alternative scartate.**
  - Documentare il bias come costo misurato: lascia nel Caso A un bias noto e correggibile.
  - q dagli stessi campioni della valutazione: è circolare.
  - Bin da 2 MeV o un modello lineare in E_n: il primo media via la struttura, il secondo la
    misura male (χ² = 420/306).
  - Un run piatto più lungo: rimandato, finché la prova di robustezza non lo chiede.
  - Correggere sin 2θ_p in relativistico sul percorso di default: cambia i record della fase 6.
    Lo scarto è ≤ 1% per cella; con `--acceptance` è già assorbito.
- **Conseguenze.**
  - Nuova app `riptide_acceptance`, che scrive `stats/<id>/acceptance.json`. Nuova opzione
    `--acceptance` di `riptide_case_a_geant4`, registrata nel record.
  - `CaseAModel` accetta un prior in θ per riga di E_n. Con righe uguali dà lo stesso risultato
    bit per bit del costruttore attuale.
  - Prova di robustezza: q stimato da metà del run (parità di `event_id`) sposta il bias del
    14 MeV a N = 1000 meno di 0.3·σ_c.
  - Lo schema delle colonne non cambia.

## D41 — Modello del Caso B (2026-11-25)

- **Contesto.** Nel Caso B la direzione Ω_n è incognita e condivisa da tutti gli eventi, mentre
  E_n è un nuisance per evento con un prior largo (assunzioni 1 + 2, campo lontano). Ogni evento
  dà (Êp, traccia 3D). Nel toy la verosimiglianza dipende da Ω_n solo tramite l'angolo θ_obs fra
  traccia e candidato, quindi si calcola una tabella per evento su θ_obs e la si interpola per
  ogni candidato.
- **Decisione.**
  - **Classe separata.** `CaseBModel` in `analysis/stats` (`case_b.hpp`, `case_b.cpp`). Non
    eredita da `CaseAModel` e non ha rami sul caso. Riusa `Grid`, `EnergyDomain`, logsumexp e la
    cinematica di D32.
  - **Griglie.** Le stesse del toy, per confrontarne i numeri alla stessa discretizzazione:
    - E_n su 500 punti in [10, 20] MeV;
    - θ_p vero su 400 punti in [0, π/2];
    - θ_obs su 801 punti in [0, π].
  - **Termine in energia, per evento:**
    g_k(θ_j) = log Σ_i exp[−½((Êp − Ep(E_i, θ_j))/σ_Ep)² + log ρ(θ_j, E_i) + log π(E_i)].
    - Il fattore relativistico della densità della traccia è
      ρ = 2m(2m + E_n) / (2m + E_n sin²θ_p)², rispetto al limite non relativistico cos θ_p/π.
      La formula è quella del toy (`log_recoil_density_ratio`) e vale |log ρ| ≲ 0.011 a 20 MeV.
    - Sta nel termine in energia perché dipende da E_n.
    - La finestra in Ep si limita alle celle che possono contare, con lo stesso limite d'errore
      dell'addendum a D35, scritto per la somma su E_n a θ_j fisso: le colonne Ep(·, θ_j) sono
      crescenti in E_n, il costruttore lo verifica, e se non vale si somma l'intera colonna.
  - **Kernel della traccia.** K(θ_obs | θ_j) è la convoluzione esatta sulla sfera dell'errore
    von Mises–Fisher di concentrazione κ = 1/σ_θ², integrata sull'azimut, come nel toy:
    K = κ/(1 − e^{−2κ}) · sin θ · exp(−2κ sin²((θ_obs − θ)/2)) · I0e(κ sin θ_obs sin θ).
    - Il kernel porta la densità cos θ/π e i pesi trapezoidali in θ.
    - Non dipende dall'evento, quindi si calcola una volta per modello.
    - I0e(x) = e^{−x} I0(x). Per I0 si usa `std::cyl_bessel_i` (C++17, `<cmath>`) finché e^{−x}
      I0(x) non perde precisione. Oltre, si usa l'espansione asintotica per x grande
      (Abramowitz–Stegun 9.7.1). La soglia si fissa nel test, su valori di riferimento di
      `scipy.special.i0e` scritti nel test con la fonte.
  - **Tabella per evento.** h_k(θ_obs) = log Σ_j exp[log K(θ_obs | θ_j) + g_k(θ_j)].
    `logLikelihood(Êp, traccia, candidati)` la interpola linearmente all'angolo
    traccia–candidato. Il costo della tabella non dipende dal numero di candidati.
  - **Accettanza.** Opzione del costruttore, spenta per default (r ≡ 1). Se è accesa, r(E_n, θ_p)
    di D40 moltiplica il termine in energia e la riga si normalizza come in D40 (D45).
- **Generatore e verosimiglianza.**
  - Il generatore sposta la traccia con una gaussiana 2D nel piano tangente, σ_θ per componente
    (`smearDirection`, come `kinematics.smear_direction` del toy). La verosimiglianza usa la vMF.
  - Le due coincidono a meno di O(σ_θ⁴), come dichiara il toy. Lo scarto è ereditato e non si
    corregge.
- **Motivo.**
  - Il toy è calibrato a ogni N con questo modello (`toy_model_findings.md` §4), quindi il port
    ha un termine di confronto numerico diretto.
  - La tabella su θ_obs rende economici sia il raffinamento sia la ricerca a molti candidati.
- **Alternative scartate.**
  - Kernel 1D gaussiano con termine specchiato: trascura la curvatura della sfera. Nel toy
    produceva un bias di μ_E di circa −σ_θ² μ_E nel Caso C.
  - Valutazione diretta per candidato senza tabella: costa n_candidati volte di più.

## D42 — Griglie di direzione, prior e stima del Caso B (2026-11-25)

- **Decisione.**
  - **Sfera.** `sphereGrid(n)`, reticolo di Fibonacci con θ_i = arccos(1 − 2(i + ½)/n) e
    φ_i = i·π(3 − √5) mod 2π. I pixel hanno area uguale, quindi il prior uniforme vale 1/n per
    pixel. Default n = 3000.
  - **Calotta.**
    - `directionCapGrid(centro, raggio, n)` è un reticolo di Fibonacci ristretto alla calotta,
      con cos α uniforme in [cos raggio, 1]. Default n = 8000.
    - `directionCapRadius` è l'angolo massimo dal MAP dei pixel grossolani con log-posterior
      > max − δ_log (δ_log = 10), più due passi della griglia grossolana, limitato a π. Il MAP si
      esclude per indice: arccos(1) arrotondato non dà 0.
  - **Prior.**
    - Uniforme sulla sfera per default, cioè uniforme per pixel anche nella calotta.
    - Il prior vMF serve solo per il test di robustezza. Si valuta su entrambe le griglie,
      perché il prior uniforme per pixel nella calotta vale solo se il prior sulla sfera è
      uniforme.
  - **Combinazione.** log p(Ω | dati) = log π(Ω) + Σ_k log L_k(Ω), con il prior contato una
    volta. Si richiede un risultato finito su tutta la sfera, altrimenti eccezione.
  - **Stima e riassunti per esperienza.**
    - Stima: MAP sulla calotta.
    - Residuo angolare: angolo fra il MAP e la verità.
    - Risoluzione dichiarata: sqrt(E_post[d²]), con d l'angolo dal MAP sulla calotta.
    - Pull: residuo / risoluzione. Per una gaussiana 2D isotropa l'rms atteso è 1 e pull·√2
      segue una Rayleigh.
    - Copertura delle regioni HPD 68/90/95%: la verità sta nella regione se il pixel più vicino
      alla verità ha densità sopra la soglia HPD.
    - `truth_in_cap`: la verità cade entro la calotta.
  - **Direzioni in `common`.** Il generatore del Caso B (in `analysis/stats`) e la risposta
    (in `analysis/io`) usano le stesse funzioni di direzione, ma stats non può includere io.
    - `Direction`, `perpendicularBasis`, `directionAroundAxis` e `smearDirection` passano senza
      modifiche in `common/include/riptide/common/direction.hpp`, e `response.hpp` lo include.
    - È un refactor a parte. I test di io esistenti restano verdi e i numeri della risposta
      non cambiano.
- **Motivo.** Sono le scelte del toy, calibrate (pixel di ~3.7° contro errori di 1–3°, quindi
  serve il raffinamento). Spostare le direzioni in `common` evita una copia e rispetta le
  dipendenze verso il basso.
- **Alternative.** HEALPix: richiede una dipendenza esterna e non serve, perché Fibonacci ha
  già area uguale. Copia delle funzioni di direzione in stats: due implementazioni dello stesso
  rumore.

## D43 — Pseudo-esperienze del Caso B (2026-11-25)

- **Decisione.** Un generatore del modello del toy in `analysis/stats`, che serve solo a
  validare il port, come D35 per il Caso A.
  - **Verità.** Per ogni esperienza Ω_n è isotropa sulla sfera. La estrae il motore dell'id
    riservato alla verità (`experimentFirstEventId`, ultimo id del blocco).
  - **Eventi.** Per ogni evento, in quest'ordine:
    1. E_n ~ U(10, 20) MeV;
    2. cos θ_CM uniforme in [−1, 1], quindi Ep e θ_p esatti (D32);
    3. azimut uniforme attorno a Ω_n;
    4. Êp = Ep + N(0, σ_Ep);
    5. traccia con `smearDirection` (σ_θ per componente).
  - **Motori.** `perEventEngine(seme, id)` per ogni evento, con lo stesso schema di id di D35.
    Seme del progetto `20260907`. Risoluzioni del toy: σ_Ep = 0.10 MeV, σ_θ = 0.08 rad.
  - **Checklist.** Le stesse coppie (N, M) del toy: (10, 400), (30, 300), (100, 200),
    (300, 100), (1000, 40) e la ripetizione a N = 1000 con M = 240 su id disgiunti.
- **Motivo.** La calibrazione con la verità dal prior è la prova più diretta del port. I suoi
  numeri si confrontano con la tabella del toy (`toy_model_findings.md` §4).

## D44 — Prestazioni del Caso B e parallelismo fra esperienze (2026-11-25)

- **Chiusura del target di D35.** Il target di 0.10 s per 1000 eventi del Caso A non è stato
  raggiunto (0.17 s) e non si insegue oltre: la checklist del Caso A gira in pochi minuti. Il
  valore resta registrato in `5_results.md`.
- **Target del Caso B.** Thread singolo, build `release`: un'esperienza a N = 1000 (tabelle,
  griglia di 3000, calotta di 8000) in ≤ 2 s.
  - La stima del costo dominante è la convoluzione: 1000 × 801 × 400 termini.
  - Prima l'implementazione semplice, poi le ottimizzazioni esatte, misurando ogni passo con un
    nuovo caso `case_b` di `riptide_stats_bench` e rilanciando gli oracoli dopo ciascuno.
  - Una banda nel kernel è ammessa solo con un limite d'errore dimostrato, come la finestra
    di D35.
- **Parallelismo.**
  - Ammesso solo fra esperienze, nell'app di calibrazione: le esperienze sono indipendenti e
    ognuna ha i suoi motori, quindi i risultati non dipendono dal numero di thread. Un test lo
    verifica bit per bit con 1 e 4 thread.
  - Dentro un'esperienza resta sequenziale. Solo `std::thread`, nessuna dipendenza nuova.
- I tempi vanno in `7_results.md`, non in un test con soglia.
- **Motivo.** Il Caso B moltiplica il costo per evento della tabella e aggiunge 11 000
  candidati per esperienza. Le esperienze della checklist sono ~10⁵ eventi in tutto, quindi il
  parallelismo fra esperienze basta e non tocca i numeri.

## D45 — Caso B su campioni Geant4 (2026-11-25)

- **Decisione.**
  - **Run.** Tre run a campo lontano piatti in 10–19.5 MeV (D25, D27), ciascuno con Ω_n fissa,
    10⁷ eventi e seme registrato:
    - +X, lato lungo del blocco;
    - +Z, lato corto;
    - (3, 2, 1)/√14, obliqua senza simmetrie del blocco.
    I config sono `config/7_ff_flat_10_19p5MeV_{x,z,obl}.json`, e ogni run passa la verifica
    `c1` (`riptide_sim_check`).
  - **Simmetria.** Le riflessioni del blocco che lasciano fissa Ω_n lungo un asse lasciano
    invariata anche la distribuzione degli eventi selezionati. Per +X e +Z il bias vettoriale
    di Ω̂ è quindi nullo per simmetria (a). L'obliqua è il test vero dell'effetto della
    selezione sulla direzione.
  - **Campioni.**
    - Catena `single_np,has_proton,contained` (D30, D34), conversione HP (D24), smearing con
      `reco_p_dir` (D31).
    - Esperienze disgiunte per N = 10 … 1000, come D37.
    - Toy gemello con la stessa verità e lo stesso numero di eventi, come D38.
  - **Modello.** Ideale (r ≡ 1) per default. Con `--acceptance <stats_id>` usa la r(E_n, θ_p)
    di D40, che è misurata lungo +X: per le altre direzioni è un'approssimazione dichiarata.
  - **Parte della selezione che dipende dalla direzione.** Il termine ε(traccia | Ω)/A(Ω) resta
    fuori. Lo si modella solo se la diagnosi lo richiede, con una decisione a parte; i run
    isotropi della fase 3 ne sono una sorgente possibile.
  - **App.** `riptide_case_b_geant4` in `analysis/apps`, un record `stats/<id>/case_b.json` per
    chiamata. Lo schema delle colonne non cambia.
- **Motivo.** Stesso schema delle fasi 5 → 6: la differenza Geant4 − toy, voce per voce, misura
  il costo del realismo del rivelatore. Le tre direzioni separano ciò che la simmetria impone
  da ciò che la selezione introduce.
- **Alternative.**
  - Una sola direzione: non distingue un bias di selezione da uno zero per simmetria.
  - Accettanza in direzione modellata subito: senza un bias misurato è complessità non
    giustificata.
- **Nota (2026-11-25, prima dei run della 7b).**
  - **Run +X.** Si riusa `6_ff_flat_10_19p5MeV_x_1e7` della fase 6: stesso config, stesso seme
    e stessi eventi di un nuovo run +X, con c1 già passata. Non esiste quindi un config
    `7_ff_flat_10_19p5MeV_x.json`; gli altri due sono `7_ff_flat_10_19p5MeV_{z,obl}.json`.
  - **Accettanza nel Caso B.** Come in D40, entra solo la forma in θ_p a E_n fissa. La densità
    del rinculo attorno a Ω si moltiplica per r(E_n,i, θ_j) e si divide per
    c_i = Σ_j m_iso(cella j; E_n,i) · r(E_n,i, θ_j), con le celle di D40. Ogni riga in E_n
    resta quindi normalizzata e A(E_n) si cancella. Con r ≡ 1 il termine è nullo e il modello
    coincide con quello ideale. Il toy gemello accetta l'evento con probabilità r / max r, come
    in D40.

## D46 — Modello del Caso C a due stadi (2026-11-26)

- **Contesto.** Nel Caso C la direzione Ω_n è condivisa, come nel Caso B, ma le energie non sono
  più indipendenti: E_n^(k) ~ N(μ_E, σ_E) troncata al dominio, con (μ_E, σ_E) incogniti. Il
  posterior congiunto ha quattro dimensioni (due di Ω_n, μ_E, σ_E). Il toy lo fattorizza in due
  stadi e la fattorizzazione costa ≤ 1.5% sulle incertezze in media, ≤ 8% nel caso peggiore (c,
  `toy_model_findings.md` §4).
- **Decisione.**
  - **Classe separata.** `CaseCModel` in `analysis/stats` (`case_c.hpp/.cpp`, e
    `case_c_posterior.hpp/.cpp` per gli stadi). Nessuna ereditarietà dal Caso B e nessun ramo
    sul caso. `CaseBModel` si usa come componente e non cambia.
  - **Stadio 1, direzione preliminare ω₀.** `CaseBModel` con il prior piatto su E_n nel dominio,
    poi `refineSharedDirection` (griglia di Fibonacci, calotta di 8000 pixel, δlog = 10, D42).
  - **Base per evento.** Con Ω_n = ω₀ e θ_obs,k l'angolo fra la traccia k e ω₀:
    base[k, i] = log Σ_j exp[kernel(θ_obs,k | θ_j) + log ρ(θ_j | E_n,i) + log N(Ê_p,k; E_p(E_n,i, θ_j), σ_Ep)],
    più l'accettanza di D40 se presente. È `CaseBModel::energyLogPosterior(Ê_p, θ_obs)` con log
    prior ≡ 0: kernel valutato in θ_obs esatto, come la base gerarchica del toy. Si calcola una
    volta per esperienza e serve alle due griglie dello stadio 2.
  - **Prior di E_n dato (μ_E, σ_E).** Gaussiana normalizzata sulla griglia di E_n con i pesi
    trapezoidali di `Grid`, quindi troncata al dominio per costruzione. Tutto in log
    (log-pesi più logsumexp): a σ_E piccolo, o con μ_E lontano dalla griglia, niente underflow.
  - **Stadio 2, iperparametri.**
    - log L(μ, σ) = Σ_k log Σ_i exp(base[k, i] + log π(E_n,i | μ, σ)).
    - log post = log prior(μ, σ) + log L(μ, σ): il prior si conta una volta.
    - Griglia globale 60 × 60: μ_E lineare sul dominio di E_n, σ_E geometrica su
      [0.01, 50] MeV, layout "ij" (μ esterno).
    - Finestra: i punti con log post > max − 10, allargati di un passo per lato (passo
      logaritmico in σ_E), tagliati al dominio della griglia globale. Griglia fine 40 × 40 nella
      finestra, su cui si calcolano tutte le uscite.
  - **Prior sugli iperparametri, propri sulla griglia.** Baseline uniforme in (μ_E, log σ_E),
    cioè peso costante per punto della griglia geometrica. Variante di robustezza: uniforme in
    σ_E, peso ∝ σ_E.
  - **Stadio 1′, direzione finale Ω̂.** Momenti della predittiva di E_n dallo stadio 2:
    m = E[μ_E], s² = E[σ_E²] + Var[μ_E]. Nuovo `CaseBModel` con prior N(m, s) troncato al
    dominio, poi di nuovo `refineSharedDirection`. Lo stadio 2 non si ripete con Ω̂: il toy fa
    lo stesso giro singolo.
  - **Mai una griglia 4D.** La marginalizzazione dello stadio 2 su Ω_n (100 pixel) e il Laplace
    congiunto del toy non entrano nella baseline. Sono diagnostiche, da portare con una
    decisione a parte solo se la checklist fallisce sulla coverage.
  - **Uscite per esperienza.**
    - μ_E: media e deviazione standard a posteriori.
    - log σ_E: media, deviazione standard e mediana.
    - Intervalli credibili centrali 0.68 / 0.90 / 0.95 di μ_E e log σ_E sulle marginali della
      griglia fine.
    - Ω̂ (MAP sulla calotta dello stadio 1′), risoluzione dichiarata σ_Ω, copertura HPD e
      `truth_in_cap` come nel Caso B.
    - Le stime puntuali sono quelle del toy: μ_E = media, σ_E = exp(mediana di log σ_E),
      Ω_n = MAP.
  - **Test di limite (Catch2).**
    - σ_E → 0: il marginale per evento coincide con la verosimiglianza del Caso B a
      E_n = μ_E, entro l'errore di quadratura dichiarato.
    - σ_E → ∞: il prior diventa piatto sul dominio e la base sommata ridà il Caso B con prior
      piatto.
    - Prior troncato: integrale trapezoidale 1 entro 1·10⁻¹², anche con μ_E fuori dal dominio.
    - R1: un log post −∞ su tutta la griglia lancia un'eccezione, mai un NaN silenzioso.
    - Valori del toy riprodotti su un'esperienza fissata, con la fonte scritta nel test.
- **Motivo.** Il toy ha già validato questa fattorizzazione (c). Riusare `CaseBModel` senza
  modificarlo lascia intatti gli oracoli della fase 7. La base per evento è il solo pezzo
  costoso e si calcola una volta.
- **Alternative scartate.**
  - Griglia 4D su (Ω_n, μ_E, σ_E): vietata, costo proibitivo per un guadagno ≤ 1.5%.
  - Stadi iterati fino a convergenza: il toy non lo fa e non c'è un difetto misurato che lo
    richieda.
  - σ_E lineare: a σ_E piccolo la risoluzione della griglia non basterebbe.

## D47 — Pseudo-esperienze e checklist del Caso C (2026-11-26)

- **Decisione.**
  - **Verità per esperienza.** μ_E ~ U(13, 17) MeV, σ_E ~ log-U(0.8, 2.4) MeV, Ω_n isotropa.
    E_n per evento dalla gaussiana troncata al dominio, la stessa densità del prior di D46.
  - **Righe.** (N, M) = (50, 200), (150, 200), (300, 100), (1000, 40) come il toy, più la
    ripetizione (1000, 240) su semi disgiunti, come nella fase 7. A N = 1000 l'esito si legge
    sulla ripetizione; la riga M = 40 serve al confronto con il toy.
  - **Robustezza al prior.** A N = 150, sulle stesse esperienze, prior uniforme in σ_E contro
    la baseline.
  - **Troncatura al dominio (roadmap).** Due righe in più, a N = 300 (M = 100) e
    N = 1000 (M = 40):
    - *bordo*: μ_E ~ U(17.5, 19) MeV, stessa σ_E, verità troncata. Il modello è esatto, quindi
      ci si aspetta calibrazione anche con molta massa troncata;
    - *sborda*: come la checklist, ma con E_n dalla gaussiana **non** troncata. Gli eventi fuori
      dal dominio restano nel campione (con μ_E = 13 e σ_E = 2.4 sono circa il 10.6%), come nello
      stress test del toy. È un test di sensibilità atteso fallire: nel toy il danno cresce
      con N.
  - **Semi.** Da (20260907, N, indice dell'esperienza), indipendenti dal numero di thread.
    Thread fra esperienze come D44, con il test bit per bit a 1 e 4 thread.
  - **App.** `riptide_case_c_calibrate`, un record `stats/<id>/case_c.json` per chiamata.
  - **Ordine della checklist.** Bias → risoluzione → pull → coverage → contrazione 1/√N →
    robustezza al prior, su μ_E, log σ_E e Ω_n.
  - **Prestazioni.** Thread singolo, build `release`: un'esperienza a N = 1000 in ≤ 5 s (b:
    circa due volte il Caso B, D44, più le griglie 2D). Prima l'implementazione semplice, poi
    si profila; dopo ogni `perf` si rilanciano gli oracoli.
- **Motivo.** Le stesse verità e righe del toy rendono il confronto voce per voce. Le due righe
  di troncatura separano l'effetto del bordo (che il modello contiene) dall'errore di modello
  (spettro fuori dal dominio, che il modello non contiene).
- **Alternative scartate.**
  - Solo la riga *sborda* del toy: non dice se il bordo in sé è trattato bene.
  - Verità fissa per tutte le esperienze: copertura misurata in un solo punto dello spazio
    degli iperparametri.

## D48 — Caso C su campioni Geant4 (2026-11-26)

- **Contesto.** I run della 7b hanno E_n incidente piatta in 10–19.5 MeV. Lo spettro degli eventi
  selezionati non è piatto: segue σ_np(E_n) · A(E_n), con A l'accettanza della catena di
  selezione. Il modello di D46 descrive lo spettro **rivelato** come una gaussiana troncata,
  senza σ_np.
- **Decisione.**
  - **Nessuna simulazione nuova.** Si ricampionano i pool selezionati dei tre run della 7b
    (+X, +Z, obliqua), catena e smearing di D45.
  - **Bersaglio del ricampionamento: lo spettro rivelato.**
    - ĝ(E) = densità empirica di `truth_En` nel pool selezionato, istogramma a passo 0.1 MeV
      sul dominio.
    - Per un'esperienza con verità (μ_E, σ_E): f(E) = gaussiana troncata al dominio, e
      l'evento si accetta con probabilità f(E)/(c ĝ(E)), c = max f/ĝ sul dominio.
    - Gli eventi accettati hanno spettro rivelato f, cioè quello che il modello assume.
    - Deconvolvere σ_np (spettro incidente ≠ rivelato) resta nella fase 9, come da roadmap.
  - **Esperienze disgiunte.** Il pool si mescola con il seme di split, registrato. Le esperienze
    di una riga percorrono la lista in ordine, e ogni evento visitato, accettato o scartato,
    non si riusa in quella riga. Le esperienze sono quindi disgiunte e ognuna è un campione
    corretto di f. Se il pool finisce, la riga tiene le esperienze complete e il record lo
    scrive (`experiments` < `experiments_planned`), come nella 7b. Righe diverse riusano il
    pool.
  - **Dominio [10, 19.5] MeV** (`EnergyDomain`, D32) nella troncatura, nella griglia di E_n e
    nella griglia di μ_E. La sorgente arriva a 19.5 MeV (D25). Con [10, 20] e μ_E = 17,
    σ_E = 2.4, circa il 4% della massa del prior cadrebbe in 19.5–20 MeV, dove non ci sono
    eventi: un bias costruito. Nella 8a (contro il toy) il dominio resta [10, 20], come il toy.
  - **Modello.** Ideale (r ≡ 1) e con `--acceptance` (r di D40), come la 7b.
  - **Toy gemello.** Stesse verità e stessi N, M per esperienza, con il dominio [10, 19.5]. Il
    toy genera direttamente da f, senza σ_np: il confronto misura il costo del rivelatore, non
    quello del ricampionamento.
  - **App.** `riptide_case_c_geant4`, un record `stats/<id>/case_c.json` per chiamata.
  - **Sottofasi.** 8a: port e checklist contro il toy; 8b: Geant4 sulle tre direzioni; 8c:
    diagnosi, solo se serve.
- **Motivo.** Il ricampionamento riusa 3 × 10⁷ eventi già validati (c1). Il bersaglio rivelato
  rende il test coerente con il modello di D46. Un bersaglio incidente misurerebbe solo l'assenza
  di σ_np, che è nota e programmata per la fase 9.
- **Alternative scartate.**
  - Nuovi run con spettro gaussiano alla sorgente: una simulazione per ogni (μ_E, σ_E), e lo
    spettro rivelato non sarebbe comunque gaussiano.
  - Ricampionare lo spettro incidente (accettanza f/max f): spettro rivelato f · σ_np · A, non
    gaussiano. Resta come diagnostica possibile della fase 9.
  - ĝ da un fit liscio: introduce una forma scelta a mano; l'istogramma con ~2000 eventi per
    bin (pool più piccolo, +Z) ha un errore ~2% sui pesi, trascurabile su gaussiane larghe
    ≥ 0.8 MeV (b).
