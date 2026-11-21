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
