# Previsioni della fase 2 — scritte prima dei run

Versione 2 (2026-11-20), **corretta dopo il pilota** e prima dei run di validazione. La
versione 1 (commit `83b959d`) era scritta dal codice del commit `7de4fec`, prima del pilota. Il
pilota ha fatto fallire un oracolo della §1: un protone ha perso 2.027 MeV in un rinculo nucleare
da scattering Coulombiano (D22). Il codice è stato corretto (commit `60c37bc`, schema 3,
`checks_version` 8), e questa versione cambia solo ciò che la sezione «Correzioni della
versione 2» elenca. Il pilota non è un run di validazione, quindi nessuna soglia è stata
adattata a un risultato di validazione. Questo file va committato prima dei run di produzione; i
risultati (`2_results.md`) citeranno l'hash di quel commit. Dopo i run le soglie non si toccano.

Le fonti dei numeri sono tre:
- le attese della fase 1b (`predictions/*_v5.json`);
- la tabella PSTAR (D20);
- i due run scratch da 10⁴ eventi dello scan del passo 4 (D21) con la config di produzione,
  cioè mono 14 MeV e piatto 10–20 MeV, seme `20260907`, 1 thread. Quei run sono cancellati e i
  loro numeri sono conservati fuori dal repo. Sono una misura (c) di questo codice, ma su un
  campione indipendente da quelli qui sotto, che hanno più eventi e più thread.

## Correzioni della versione 2

- §1, bilancio dei protoni contenuti: l'oracolo `contained_balance_MeV.min` ≥ −10⁻⁹ MeV vale
  ora sui protoni senza interazioni adroniche **e senza rinculi nucleari**
  (`truth_p_n_nuclear_recoil` = 0). Nella versione 1 l'esclusione mancava, e nel pilota a
  8 thread il minimo era −2.027 MeV (evento 246185).
- §1: il conteggio `deposits.protons_nuclear_recoil` ha una previsione propria; fra i run di kB
  deve essere identico.
- Run: schema 3 al posto di 2; i due piloti si ripetono con lo schema 3, con i `run_id`
  `pilot2_v3_*`. I piloti v2 restano sull'SSD e non si verificano più: `checks_version` 8 salta
  le verifiche dei depositi sui run con schema 2.
- Nient'altro cambia: le soglie statistiche, le bande PSTAR e le attese delle §2–§7 sono quelle
  della versione 1. Il campione PSTAR perde circa un protone su 10⁵, che non sposta le medie in
  modo visibile (b).

Etichette: **(a)** fatto stabilito (sorgenti, dati, costruzione del codice); **(b)** stima con
ipotesi esplicite; **(c)** verificato con Monte Carlo; **(d)** ipotesi da mettere alla prova.

## Run previsti

Build `release`, Geant4 11.4.1, `QGSP_BIC_HP`, G4NDL 4.7.1, cubo da 100 mm, fascio lungo z
(L = 100 mm), seme del progetto `20260907`. Ogni run gira sotto `scripts/run_capped.sh -m 5G`,
uno alla volta, dopo `scripts/data_root.sh check`. La verifica è `riptide_sim_check --check-id
c1 --bins-u 50 --bins-z 50`.

| Run (`run_id`) | Config | Eventi | Thread | Ruolo |
|---|---|---|---|---|
| `pilot2_v3_mono_14MeV_1e5_1t` | `mono_14MeV.json` | 10⁵ | 1 | prestazioni, byte/evento |
| `pilot2_v3_mono_14MeV_4e5_8t` | `mono_14MeV.json` | 4·10⁵ | 8 | prestazioni, byte/evento |
| `2_mono_14MeV_1e7` | `mono_14MeV.json` | 10⁷ | 8 | riferimento, tutto HP |
| `2_flat_10_20MeV_1e7` | `flat_10_20MeV.json` | 10⁷ | 8 | sorgente di produzione |
| `2_flat_19p5_20MeV_1e7` | `flat_19p5_20MeV.json` | 10⁷ | 8 | rampa HP/CHIPS (D16) |
| `2_mono_14MeV_kB0p145_1e6` | mono, kB = 0.145 | 10⁶ | 8 | sensibilità a kB (D18) |
| `2_mono_14MeV_kB0p155_1e6` | mono, kB = 0.155 | 10⁶ | 8 | sensibilità a kB (D18) |
| `2_mono_14MeV_kB0p165_1e6` | mono, kB = 0.165 | 10⁶ | 8 | sensibilità a kB (D18) |

Tutte le config hanno `config_version` 2, schema 3, funzione di passo del protone (0.2, 0.1 mm)
(D21) e kB = 0.155 mm/MeV (D18). Fanno eccezione i due run di sensibilità a 0.145 e 0.165: le
loro config si generano fuori dal repo da `mono_14MeV.json`, cambiando solo
`scintillator.birks_kB_mm_per_MeV`, e la copia completa resta registrata in `run.json`. Il
pilota non è un run di validazione: i suoi test non entrano nei risultati.

```bash
scripts/run_capped.sh -m 5G -- build/release/sim/riptide_sim --config config/<c>.json \
  --events 10000000 --threads 8 --run-id 2_<c>_1e7
scripts/run_capped.sh -m 5G -- build/release/sim/riptide_sim_check --run-id 2_<c>_1e7 \
  --check-id c1 --bins-u 50 --bins-z 50
```

## 1. Oracoli deterministici (soglie esatte)

Valgono per tutti i run di validazione (i tre da 10⁷ e i tre di sensibilità).

| Verifica | Soglia | Etichetta |
|---|---|---|
| `results.counters.must_be_zero.total` | = 0 | (d), solo un bug li rende ≠ 0 |
| `results.deposits.counters.must_be_zero.total` | = 0 | (d), come sopra |
| `deposits.skipped`, `deposits.pstar.skipped` | `null` | (a), schema 3 |
| Σ `deposits.event_classes` | = `deposits.events` = `results.events` | (a) |
| `deposits.contained_balance_MeV.max` | ≤ 0.01 MeV (slack di D19) | (d); scan: 1.9·10⁻⁴ MeV |
| `deposits.contained_balance_MeV.min` | ≥ −10⁻⁹ MeV, protoni senza rinculi nucleari (v2) | (d); scan: −5.3·10⁻¹⁵ MeV |
| `deposits.protons_nuclear_recoil`, run da 10⁷ | fra 1 e 150 | (d), sotto (v2) |
| `deposits.pstar.samples` + `below_table` | = `contained_balance_MeV.samples` | (a), stessi protoni |
| `deposits.pstar.sha256` | `eec195fe…081b` di D20 | (a) |
| `deposits.setup.birks_kB_mm_per_MeV` | kB della config del run | (a) |
| oracoli 1b: `Ep_residuals.max_abs`, `model_mix.Ep_residuals.max_abs`, offset HP − fisico, `awr`, `model_mix.observed_non_hp` = 0 a 14 MeV | come §1 di `1b_predictions.md`, stessi intervalli | come lì |

`Edep_overflow` (eventi con `det_Edep` sopra E_n,max) non è un oracolo, perché `checks.json` ne
riporta solo il conteggio e non i valori. Le catture con il loro Q e l'energia creata dall'offset
HP (Σε ≲ 6·10⁻⁴·E_n, D14) lo rendono positivo per costruzione. L'attesa è nella §5.

La banda di `protons_nuclear_recoil` viene da un solo caso, 1 su 137 806 protoni contenuti nei
due piloti v2 a 14 MeV (D22). Il limite di Poisson al 90% su quel tasso, [0.05, 4.7] casi per
1.38·10⁵ protoni, e i circa 2.7·10⁶ protoni contenuti di un run da 10⁷ (scan, §5) danno da 1 a
circa 95 casi. Il limite alto è portato a 150 perché la sorgente piatta arriva a 20 MeV, dove il
cammino del protone è più lungo (b). Non è un oracolo esatto: fuori banda è una previsione
fallita, non un bug.

### Run di sensibilità a kB: stessi eventi

Ogni evento è seminato solo da (event_id, seme del run), qualunque sia il thread
(`event_seeds.hpp`), e kB entra solo nella luce visibile: con kB = 0 il test d'integrazione
verifica `det_Edep` identico bit per bit (a). Quindi, fra i tre run `2_mono_14MeV_kB*_1e6`:

| Verifica | Soglia | Etichetta |
|---|---|---|
| colonna `Edep` di `deposit_spectra.csv`, tutto `proton_exits.csv` | identici byte per byte | (d), atteso per costruzione |
| `event_classes`, `contained_balance_MeV`, `protons_nuclear_recoil` (v2), `pstar.range_mm`, `Edep_overflow` | identici | (d), come sopra |
| `pstar.light_MeV.relative_residual` a 0.145 e a 0.165, meno quello a 0.155 | \|Δ\| < 0.01% (assoluto) | (b), sotto |

L'ultima soglia vale perché il residuo contro PSTAR viene dalla differenza fra lo S di Geant4 e
quello di PSTAR (D21). Il suo peso nella luce scala come kB·S/(1 + kB·S), che cambia di pochi
punti percentuali fra 0.145 e 0.165. Il residuo di −0.06% dovrebbe quindi muoversi di qualche
10⁻³ %, non di più (b).

## 2. Test statistici della fase 1b sui nuovi eventi

La prima interazione del neutrone non dipende da ciò che cambia nella fase 2: kB e la luce dei
protoni (D21) non toccano il trasporto (a). Le attese restano quelle della fase 1b: rigenerate
oggi con `riptide_sim_predict` dal commit `7de4fec` per le tre config a 10⁷ eventi, coincidono
con `predictions/*_1e7_v5.json` in tutto tranne `checks_version` (a). Valgono quindi, invariati,
il piano dei test, le soglie e le attese della §2 di `1b_predictions.md`:
- test contro `jeff33`;
- Bonferroni con α di famiglia 0.01 per run (3.261 / 3.317 sui pull);
- confronto con `endfb80` della §3 come attese e non come soglie.

I run da 10⁷ sono eventi nuovi rispetto a quelli della fase 1b solo nella parte che non riguarda
il neutrone. Con lo stesso seme, la prima interazione di ogni evento dovrebbe coincidere con
quella del run 1b corrispondente. Di conseguenza `fractions.csv`, `spectrum_u.csv`,
`depth_z.csv`, `legendre.csv` e `model_mix.csv` dovrebbero essere identici a quelli di `c3` dei
run 1b (d). È una previsione e non una soglia: lo stato del generatore dopo la prima interazione
dipende dalla fisica dei secondari, e se le colonne cambiano la causa va spiegata, ma il test
statistico resta quello sopra.

## 3. Confronto con PSTAR (soglie a banda)

A 10⁷ eventi i protoni confrontati sono circa 2.8·10⁶ (scan: 2768 su 10⁴ a 14 MeV) e l'errore
statistico del residuo relativo di luce scende a circa 3·10⁻⁵ %. Il residuo è sistematico
(Geant4 contro PSTAR, D21), non rumore. Per questo le soglie sono bande, fissate sullo scan.

**Luce, `pstar.light_MeV.relative_residual`:**

| Run | Soglia | Scan (c) |
|---|---|---|
| mono 14 MeV | [−0.070%, −0.055%] | −0.0616% ± 0.0010% |
| piatto 10–20 MeV | [−0.070%, −0.052%] | −0.0600% ± 0.0010% |
| bordo 19.5–20 MeV | [−0.070%, −0.050%] | non misurato |

Per il bordo si prevede un valore un po' meno negativo del piatto, perché ha più protoni sopra
14 MeV, dove il residuo per bin sale verso −0.05% (b).

**Luce per bin di E_p.** Il residuo per bin è una funzione liscia di E_p, uguale nello scan mono
e in quello piatto dove i bin si sovrappongono (per esempio −0.0735% a 7.6–8.0 MeV in entrambi)
(c):

| E_p [MeV] | 0–0.3 | 1 | 2 | 3.3 | 5 | 7.7 | 10 | 14 | 19.5 |
|---|---|---|---|---|---|---|---|---|---|
| residuo (c) | +1.1% / +0.8% | +0.16% | +0.08% | ≈ 0 | −0.055% | −0.0735% (minimo) | −0.070% | −0.060% | −0.048% |

- **Sotto 3 MeV, positivo e crescente verso E_p → 0:** è lo stopping nucleare, incluso nello S
  totale di PSTAR ed escluso dalla tabella di luce di Geant4 (D21). Da PSTAR da sola,
  L_elettronico/L_totale − 1 vale +0.17% a 1 MeV (a). Nei primi bin pesa anche la forma dello
  spettro dentro il bin.
- **Sopra 3 MeV, negativo:** lo S di Bethe di Geant4 è un po' più alto di quello di PSTAR (d,
  D21).

Soglia per i bin con `Ep_low_MeV` ≥ 2 MeV e almeno 100 protoni: il residuo osservato coincide con
il valore dello scan per lo stesso bin entro **0.005% assoluto**. Il riferimento è lo scan mono
per il run mono e lo scan piatto per il piatto e il bordo, che hanno gli stessi bordi di bin fino
a 19.6 MeV. Lo scan aveva errore ≤ 0.0007% per bin; il resto della banda copre le differenze di
spettro dentro il bin fra i run. Per i bin del bordo sopra 19.6 MeV vale la continuità: entro
0.005% dal bin precedente. Sotto 2 MeV niente soglia; si verifica solo il segno (positivo) e la
decrescenza con E_p (d).

**Range, `pstar.range_mm.relative_residual` (`truth_p_length` contro R_CSDA):**

| Run | Soglia | Scan (c) |
|---|---|---|
| tutti e tre | [−0.20%, 0.00%] | mono −0.10% ± 0.03%, piatto −0.07% ± 0.03% |

Per bin (attesa, non soglia, (c) dallo scan): il residuo è nettamente positivo sotto 3 MeV (da +8%
nel primo bin a +0.4% a 2 MeV, nello scan entro 0.8σ–1σ di bin; la causa non è stata studiata,
(d)), e fra −0.4% e +0.4% sopra, compatibile con circa −0.1%. A 10⁷ eventi l'errore per bin scende
a circa 0.005%: la forma diventa una misura, da descrivere nei risultati.

## 4. Sensibilità a kB (attese, non soglie)

L(E_p; kB) = ∫₀^E_p dE / (1 + kB·S(E)) con lo S totale di PSTAR, interpolato log-log e integrato
con trapezi fini (calcolo fuori dal repo). Variazione relativa rispetto a kB = 0.155 (b):

| E_p [MeV] | 1 | 2 | 5 | 10 | 14 | 20 |
|---|---|---|---|---|---|---|
| L(E_p; 0.155) / E_p | 0.126 | 0.182 | 0.294 | 0.405 | 0.464 | 0.528 |
| kB = 0.145 | +5.88% | +5.40% | +4.51% | +3.68% | +3.26% | +2.82% |
| kB = 0.165 | −5.26% | −4.87% | −4.13% | −3.41% | −3.05% | −2.66% |

- Previsione (b): `pstar.light_MeV.observed_mean` cambia fra i run di sensibilità come la media
  di queste variazioni sui protoni contenuti, dato che il residuo contro PSTAR è quasi costante
  (§1, ultima soglia).
- La media di `det_Evis` (dalla colonna `Evis` di `deposit_spectra.csv`) si sposta di meno: anche
  elettroni e ioni entrano in `det_Evis`, e passano da `G4EmSaturation` con lo stesso kB. Il
  valore si riporta come sistematica di kB, senza soglia.

## 5. Attese che non sono soglie

**Classi dell'evento (D19)**, dallo scan a 10⁴ eventi (c; errore binomiale circa
√(f(1 − f)/10⁴), cioè 0.3–0.5 punti percentuali):

| Classe | mono 14 MeV | piatto 10–20 MeV |
|---|---|---|
| `no_interaction` | 37.3% | 35.7% |
| `nc_inelastic_no_np` | 17.6% | 18.6% |
| `nc_elastic_only` | 17.0% | 18.5% |
| `mixed` | 14.1% | 14.2% |
| `single_np` | 9.8% | 9.1% |
| `multiple_np` | 3.9% | 3.8% |
| `other` | 0.22% | 0.17% |

Il bordo non è stato misurato. Per la prima interazione la §2 di `1b_predictions.md` dà meno
n–p (0.165 contro 0.230) e più n–C, quindi ci si aspettano meno `single_np` (b).

Altre attese dallo scan (c, 10⁴ eventi):

| Grandezza | mono 14 MeV | piatto 10–20 MeV |
|---|---|---|
| `np_events` | 2808 | 2723 |
| protoni contenuti confrontati con PSTAR | 2768 | 2685 |
| `protons_exiting` | 16 | 17 |
| `protons_hadronic` | 24 | 21 |
| `events_secondary_n_hadronic` | 999 | non registrato |
| `Edep_overflow` | 267 (2.7%) | 2 |
| `Evis_overflow` | 0 | 0 |
| `expected_rare.p_Edep_above_Ep` | 4 | non registrato |

Con una sorgente monocromatica `det_Edep` supera E_n quando l'energia è contenuta e l'offset HP
o la cattura aggiungono energia (D14). Con la sorgente piatta succede solo per E_n vicino a
20 MeV, quindi di rado (b).

## 6. Residui noti non modellati

| Effetto | Grandezza | Dove pesa | Etichetta |
|---|---|---|---|
| Offset HP della massa di H-1 (D14) | E_p +5.4·10⁻⁴ relativo | E_p, `det_Edep` sopra E_n | (a) |
| Eccesso di `ionIoni` a fine range per i rinculi C12 / Be9 (D14) | fino al 15% di K, 0.74 MeV su 492 MeV di rinculi C12 nel run d'integrazione | `det_Edep`, poco la luce | (c), causa (d) |
| Stopping nucleare escluso dalla luce (D21) | +0.17% a 1 MeV, +0.03% a 13 MeV contro PSTAR totale | residuo di luce sotto 3 MeV | (a) |
| S di Bethe di Geant4 contro PSTAR (D21) | circa −0.07% di luce sopra 3 MeV | residuo di luce | (d) |
| kB oltre 500 MeV/cm, fuori dal campo delle misure di Pöschl (D18) | non stimato | luce a fine range | (d) |
| C-13 contato due volte nel modello HP (D15) | circa +0.5% sul non elastico del carbonio | frazioni n–C | (a) |
| n–d nella verità come n–p (D17) | n_H2/n_H1 = 1.15·10⁻⁴ | conteggio n–p | (a) |

## 7. Prestazioni attese (dallo scan, 2026-11-20)

| Grandezza | Valore | Etichetta |
|---|---|---|
| tempo per 10⁴ eventi a 1 thread, mono 14 MeV | 2.0–2.3 s (inizializzazione compresa) | (c) |
| RSS massimo, sim a 1 thread | < 92 MiB | (c) |
| `riptide_sim_check` su 10⁴ eventi | 0.22 s, 65 MiB | (c) |
| byte per evento su disco (10⁴ eventi, 1 thread) | circa 118 | (c) |

Stime per un run da 10⁷ eventi a 8 thread (b): tempo di 4–8 minuti, con la stessa scalatura
del pilota 1b (5.5·10³ eventi/s per thread); RSS sotto 0.3 GB; circa 1.2 GB su disco; la
verifica in circa un minuto. In tutto, i tre run e i tre di sensibilità occupano circa 4 GB
(1.7 TB liberi sull'SSD). Il pilota sostituisce queste stime:
- se il byte per evento misurato supera `kBytesPerEvent` (256, provvisorio), il valore si
  aggiorna prima dei run di produzione, con un commit `fix(sim)` separato;
- altrimenti resta, con il valore misurato registrato in `2_results.md`.

## 8. Se una soglia fallisce

Vale la procedura della §6 di `1b_predictions.md`:
- le soglie non si cambiano e non si ripete il run con un altro seme;
- la causa si cerca sul run esistente;
- un oracolo esatto che fallisce è un bug finché non si dimostra il contrario;
- una convenzione di Geant4 non modellata porta a una decisione, a `checks_version` + 1 e a
  un'altra verifica sugli stessi eventi;
- nei risultati si riportano entrambi gli esiti.

In più, per le bande PSTAR della §3, un fallimento si studia bin per bin di E_p, separando
l'effetto dello stopping nucleare (sotto 3 MeV) da quello di Bethe (sopra). Le bande non si
allargano: se lo scan da 10⁴ eventi non era rappresentativo, nei risultati lo si scrive come
previsione fallita.
