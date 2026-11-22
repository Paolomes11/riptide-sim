# Risultati della fase 3 — geometria del blocco e sorgente a campo lontano

Data: 2026-11-22. Previsioni e soglie: `docs/validation/3_predictions.md` versione 2, committata
in `ff6f837` **prima** dei run. La versione 1 (`8429cc3`) differiva solo per il riferimento della
luce per bin, che la v2 ha portato nel repo. Fra le previsioni e i run è entrato solo un commit di
documentazione (`ded6c88`, stato del progetto); il codice è quello di `a15d899`.

**Esito (`c1`):** nei tre run da 10⁷ eventi passano tutte le soglie:
- oracoli esatti della §2;
- test χ² e pull contro `jeff33` (Bonferroni, α di famiglia 0.01 per run);
- bande di luce e range contro PSTAR e luce per bin.

Le attese che non sono soglie sono quasi tutte confermate. Le eccezioni sono quattro, discusse
sotto e riportate nei punti aperti:
- frazione n–C non elastica lungo +X a +3.73σ da `jeff33` (§2), non una soglia;
- previsione (d) sul range non confermata (§4);
- `Edep_overflow` isotropo mono sopra la stima del pilota (§5);
- rapporto `single_np` / (`mixed` + `multiple_np`) previsto in modo sbagliato (§5).

## Run

Build `release` al commit `ded6c88` (`git describe` v0.4.0-15-gded6c88), Geant4 11.4.1,
`QGSP_BIC_HP`, G4NDL 4.7.1, 8 thread, seme 20260907, `scripts/run_capped.sh -m 5G`, un run alla
volta, dopo `scripts/data_root.sh check`. Config `config/3_ff_*.json` (`config_version` 3), schema
4, `checks_version` 10. Dati sull'SSD in `sim/<run>` e `validation/<run>/c1`. Comandi come nella
sezione «Run previsti» delle previsioni.

| Run | Eventi | Tempo | RSS massimo | Su disco | Byte/evento |
|---|---|---|---|---|---|
| `3_ff_mono_14MeV_x_1e7` | 10⁷ | 165.4 s | 169 MB | 1.10 GB | 110.1 |
| `3_ff_mono_14MeV_isotropic_1e7` | 10⁷ | 141.8 s | 177 MB | 1.11 GB | 111.0 |
| `3_ff_flat_10_19p5MeV_isotropic_1e7` | 10⁷ | 145.5 s | 180 MB | 1.19 GB | 118.7 |

Tempo da `run.json` (inizializzazione compresa), RSS da `/usr/bin/time -v`, byte dagli
`events_t*.root` (c).

**Prestazioni contro le stime** (§8 delle previsioni):
- **velocità:** 6.0–7.1·10⁴ eventi/s, cioè 7.6–8.8·10³ per thread. È più della stima, che
  scalava il pilota a 1 thread con il rapporto 8 thread / 1 thread della fase 2; il tempo
  (2.4–2.8 minuti) cade comunque nella banda stimata di 2–6 minuti;
- **memoria:** RSS sotto 0.2 GB, come stimato;
- **byte per evento:** 110.1 / 111.0 / 118.7 contro 110 / 112 / 119 dei piloti. `kBytesPerEvent`
  = 256 B resta, con un margine di più di 2 volte;
- **disco:** 3.40 GB in tutto, contro i circa 3.5 GB stimati.

`riptide_sim_check --check-id c1 --bins-u 50 --bins-z 50`: 49 s / 38 s / 40 s, RSS al più 117 MB
(c).

## 1. Oracoli deterministici — (a)

| Verifica | Soglia | +X | isotropa mono | isotropa piatta |
|---|---|---|---|---|
| `must_be_zero.total` (prima interazione) | = 0 | 0 | 0 | 0 |
| `deposits.must_be_zero.total` | = 0 | 0 | 0 | 0 |
| `deposits.skipped`, `deposits.pstar.skipped` | `null` | `null` | `null` | `null` |
| Σ classi = eventi | 10⁷ | sì | sì | sì |
| `geometric_acceptance.rays` | = eventi | 10⁷ | 10⁷ | 10⁷ |
| `no_interaction` ≥ `rays` − `hits` | vero | vero | vero | vero |
| corda +X: min, max, sd | 60 mm entro 10⁻⁹, sd < 10⁻⁹ | 60 − 7·10⁻¹⁵ mm, sd 0 | — | — |
| corda isotropa, max | ≤ 66.33 mm | — | 65.61 mm | 65.60 mm |
| `contained_balance_MeV.max` | ≤ 0.01 MeV | 8.5·10⁻⁴ | 7.6·10⁻⁴ | 8.1·10⁻⁴ |
| `contained_balance_MeV.min` | ≥ −10⁻⁹ MeV | −5.3·10⁻¹⁵ | −5.3·10⁻¹⁵ | −7.1·10⁻¹⁵ |
| `pstar.samples` + `below_table` = protoni del bilancio | uguali | 1 069 971 + 71 | 223 900 + 21 | 218 002 + 24 |
| `pstar.sha256` | `eec195fe…081b` | sì | sì | sì |
| `birks_kB_mm_per_MeV` | 0.155 | 0.155 | 0.155 | 0.155 |
| `Ep_residuals.max_abs` | < 10⁻⁹·E_n,max | 1.8·10⁻¹⁴ MeV | 2.0·10⁻¹⁴ MeV | 2.8·10⁻¹⁴ MeV |
| `model_mix.observed_non_hp` | = 0 | 0 | 0 | 0 |
| offset HP relativo, min / max (·10⁻⁴) | nell'intervallo della §2 | 5.36581 / 5.48614 | 5.36581 / 5.48614 | 5.33517 / 5.50151 |
| `hp_elastic_H1.awr` | 0.999167 | 0.999167 | 0.999167 | 0.999167 |

Intervalli dell'offset: [5.36581, 5.48614]·10⁻⁴ per i due mono, [5.33500, 5.50179]·10⁻⁴ per il
piatto. Nei due run mono minimo e massimo coincidono con gli estremi dell'intervallo alle cifre riportate,
come nei piloti.

**Tutti PASS.**

## 2. Test statistici contro `jeff33` — (c)

Soglie: p > 1.000·10⁻³ e |pull| < 3.291 lungo +X; p > 9.09·10⁻⁴ e |pull| < 3.317 nei due run
isotropi.

| Test | +X | isotropa mono | isotropa piatta |
|---|---|---|---|
| spettro u, χ²/49 (p) | 49.57 (0.450) | 37.84 (0.876) | 45.13 (0.631) |
| profondità s = z/L_i, χ²/49 (p) | 43.99 (0.676) | 36.66 (0.903) | 38.38 (0.863) |
| frazioni, χ²/4 (p) | 16.22 (**0.00274**) | 1.43 (0.838) | 1.57 (0.815) |
| pull di Legendre a₁…a₆ | 0.96, −1.01, −1.95, −0.36, 0.70, 2.00 | 0.11, −0.12, 0.00, 0.95, −0.81, 1.59 | −0.53, 0.39, −0.30, 0.62, 1.33, −0.32 |
| pull dell'accettanza | 1.435 | 0.642 | 1.205 |
| pull della corda media | — (oracolo, §1) | 0.690 | 0.534 |

**Tutti PASS.** Il più vicino alla soglia è il χ² delle frazioni lungo +X, con p = 0.0027 contro
α = 0.001.

Il χ² delle frazioni ha 4 gradi di libertà, non i 3 scritti nelle previsioni. A 10⁷ eventi la
categoria `altro` ha abbastanza conteggi attesi (circa 170 lungo +X, 37 e 54 negli isotropi) ed
entra nel test. La soglia è sulla p-value, quindi non cambia (a).

Valori misurati di a₁: −1.425·10⁻² ± 5.7·10⁻⁴ (+X), −1.466·10⁻² ± 1.2·10⁻³ (isotropa mono),
−1.512·10⁻² ± 1.2·10⁻³ (isotropa piatta). Il segnale è a 25σ lungo +X e a 12σ negli isotropi,
come previsto. Accettanza: 6 368 380 raggi sul blocco contro 6 366 198 attesi lungo +X;
4 052 213 e 4 053 088 contro 4 051 217 negli isotropi. Corda media: 17.1465 e 17.1456 mm contro
17.142857 mm (Cauchy), sd 10.490 e 10.491 mm contro 10.49 della previsione.

**Osservazione: n–C non elastico lungo +X.** Pull per categoria contro `jeff33`, nell'ordine
[sconosciuta, n–p, n–C el., n–C non el., altro]: [−1.59, −0.19, −0.48, **+3.73**, −1.51].
- osservati 666 827 eventi n–C non elastici contro 663 889 attesi, cioè +0.44%;
- i pull per categoria non sono test del piano: il test è il χ², che passa;
- gli altri run da 10⁷ con lo stesso modello di n–C (c):

| Run | pull n–C non el. contro `jeff33` |
|---|---|
| `3_ff_mono_14MeV_x_1e7` (blocco, +X) | +3.73 |
| `3_ff_mono_14MeV_isotropic_1e7` | −1.09 |
| `3_ff_flat_10_19p5MeV_isotropic_1e7` | +0.54 |
| fase 1b/2, mono 14 MeV, fascio sottile, cubo da 100 mm | −2.25 (1 462 545 contro 1 465 065, −0.17%) |

Il segno cambia da un run all'altro, quindi non è un bias fisso della sezione d'urto del carbonio
nel modello. Fra +X e il fascio sottile dello stesso mono 14 MeV, però, la differenza è di circa
4σ, e oggi la differenza di geometria non la spiega. Le cause da guardare sono due, sullo stesso
run e senza nuovo seme:
- la distribuzione della profondità dei vertici n–C non elastici;
- la loro vicinanza alle facce del blocco (rientro di secondari, conteggio della prima
  interazione).

Va fatto prima di usare le frazioni nell'analisi. Etichetta (d): ipotesi da verificare.

## 3. Confronto con `endfb80` — attese, non soglie (c)

| Verifica | +X | isotropa mono | isotropa piatta |
|---|---|---|---|
| frazioni, χ²/4 (p) | 284.5 (2.4·10⁻⁶⁰) | 29.3 (6.9·10⁻⁶) | 33.4 (1.0·10⁻⁶) |
| non centralità prevista | circa 170 | circa 35 | circa 35 |
| spettro u, p | 0.319 | 0.86 | 0.58 |

Le frazioni falliscono nei tre run, come previsto: viene dalle valutazioni del carbonio.

Pull per categoria lungo +X, nello stesso ordine della §2: [−2.31, −1.03, −7.26, +14.88, −5.33].
Rispetto ai pull contro `jeff33`, gli spostamenti sono quelli previsti in modulo:
- n–C non elastico 11.15σ, contro 11.05σ;
- n–C elastico 6.78σ, contro 6.8σ;
- n–p 0.84σ, contro 0.85σ.

La non centralità osservata lungo +X (circa 280, χ² meno i gradi di libertà) è sopra quella prevista (circa 170). È coerente
con il +3.73σ della §2, che si somma allo spostamento fra le librerie sulla stessa categoria.

Pull di Legendre lungo +X contro `endfb80`: [0.34, −1.14, −2.86, 0.11, 0.66, 2.00]. a₃ passa da
−1.95 a −2.86, cioè 0.91σ, come previsto. Su n–p e sugli a_l le due revisioni di H-1 restano
entro 1σ: a 10⁷ eventi la simulazione non le distingue, come nella fase 1b (D13).

## 4. Confronto con PSTAR — (c)

| Verifica | Soglia | +X | isotropa mono | isotropa piatta |
|---|---|---|---|---|
| luce, residuo relativo medio | [−0.070%, −0.055%] / [−0.070%, −0.052%] piatto | −0.0619% ± 0.0001% | −0.0618% | −0.0602% |
| luce per bin, max \|run − riferimento\| | < 0.005% | 0.00013% (42 bin) | 0.00012% (42 bin) | 0.00070% (43 bin) |
| range, residuo relativo medio | +X [−0.125%, −0.035%]; iso mono [−0.235%, −0.041%]; piatto [−0.125%, +0.091%] | −0.1083% ± 0.0015% | −0.1060% ± 0.0032% | −0.1026% ± 0.0035% |

Luce per bin: bin con `Ep_low` ≥ 2 MeV e almeno 100 protoni, accoppiati al riferimento per
`Ep_low` entro 10⁻⁶ MeV. Lungo +X lo scarto più grande è nel bin con `Ep_low` = 2.52 MeV.

**Tutti PASS.** La luce per protone non dipende dalla geometria: i residui medi coincidono con la
fase 2 (−0.0620% mono, −0.0600% piatto) entro 0.0002%.

**Previsione (d) sul range non confermata.** La previsione era un residuo meno negativo della
fase 2 (−0.106% mono, −0.103% piatto), perché pesa di più la parte sotto 3 MeV. Osservato:
- +X: −0.1083%, più negativo;
- isotropa mono: −0.1060%, uguale;
- isotropa piatta: −0.1026%, uguale entro l'errore.

L'esclusione dei protoni che escono non sposta il residuo medio del range nel verso previsto. La
causa non è stata studiata. Non è una soglia, e le bande passano.

## 5. Attese che non sono soglie — (c)

**Protoni che escono** (`protons_exiting` / `np_events`, protone della prima n–p):

| | Banda (±3σ del pilota) | Osservato |
|---|---|---|
| +X | [2.8%, 3.9%] | 3.171% ± 0.017% |
| isotropa mono | [2.8%, 5.2%] | 4.365% ± 0.042% |
| isotropa piatta | [3.5%, 6.2%] | 4.758% ± 0.044% |

Le tre frazioni sono nelle bande, e il piatto supera il mono isotropo, come previsto. Fra
isotropa mono e +X non era previsto il verso: l'isotropa è più alta di 26σ. È l'effetto dei bordi
della fase 3: nella fase 2 (cubo da 100 mm) la frazione era 0.47%.

Distribuzione per E_p, da `proton_exits.csv`, in cinque fette uguali di [0, E_n,max]:

| Fetta di E_p (mono) | +X | isotropa mono | Fetta di E_p (piatto) | isotropa piatta |
|---|---|---|---|---|
| 0–2.8 MeV | 0.29% | 0.25% | 0–3.9 MeV | 0.50% |
| 2.8–5.6 MeV | 1.53% | 1.53% | 3.9–7.8 MeV | 2.70% |
| 5.6–8.4 MeV | 3.28% | 3.58% | 7.8–11.7 MeV | 6.42% |
| 8.4–11.2 MeV | 5.10% | 6.36% | 11.7–15.6 MeV | 11.1% |
| 11.2–14 MeV | 5.47% | 9.74% | 15.6–19.5 MeV | 16.9% |

La frazione cresce con E_p, cioè col range. Lungo +X si appiattisce sopra 8 MeV, mentre nei run
isotropi continua a salire. Interpretazione (b): i protoni più energetici escono in avanti, a
θ_p piccolo. Lungo +X la direzione in avanti è il lato da 60 mm, nei run isotropi è casuale
rispetto al blocco e cade spesso su un lato da 20 mm. Per la selezione dell'analisi conta che la
perdita dipende da E_p e dalla direzione: sopra 11.7 MeV, nella sorgente di produzione, esce più
di un protone su dieci.

**Classi dell'evento (D19)**, contro i piloti da 10⁵ eventi:

| Classe | +X (pilota) | isotropa mono (pilota) | isotropa piatta (pilota) |
|---|---|---|---|
| `no_interaction` | 71.686% (71.72 ± 0.14) | 93.917% (93.72 ± 0.08) | 93.713% (93.67 ± 0.08) |
| `nc_elastic_only` | 10.125% (10.10 ± 0.10) | 2.241% (2.37 ± 0.05) | 2.439% (2.43 ± 0.05) |
| `single_np` | 8.523% (8.52 ± 0.09) | 1.861% (1.84 ± 0.04) | 1.812% (1.85 ± 0.04) |
| `nc_inelastic_no_np` | 7.078% (7.08 ± 0.08) | 1.488% (1.55 ± 0.04) | 1.545% (1.58 ± 0.04) |
| `mixed` | 1.690% (1.69 ± 0.04) | 0.306% (0.33 ± 0.02) | 0.308% (0.29 ± 0.02) |
| `multiple_np` | 0.897% (0.89 ± 0.03) | 0.187% (0.19 ± 0.01) | 0.181% (0.19 ± 0.01) |
| `other` | 0.002% (0.003 ± 0.002) | 0.000% (0) | 0.001% (0) |

I run sono entro l'errore dei piloti tranne due classi dell'isotropa mono, che distano circa
2.5σ del pilota: `no_interaction` (+0.20 punti) e `nc_elastic_only` (−0.13 punti). Il pilota è
contenuto nei primi 10⁵ eventi del run, quindi il confronto non è fra campioni indipendenti.

**Previsione (b) sbagliata nel testo:** le previsioni dicevano che lungo +X `single_np` vale
circa 4.6 volte (`mixed` + `multiple_np`). Il run dà 3.29, e la tabella dei piloti nello stesso
file dava già 3.30 (8.52 / 2.58). Il 4.6 è un errore di conto della previsione, non un effetto
fisico. Resta vera la conclusione qualitativa: il blocco riduce le interazioni multiple, perché
nella fase 2 il rapporto era 0.54.

**Altre attese**, scalate dai piloti a 10⁷ eventi:

| Grandezza | +X: attesa / osservato | isotropa mono | isotropa piatta |
|---|---|---|---|
| `np_events` | 1.11·10⁶ / 1 110 802 | 2.36·10⁵ / 235 355 | 2.33·10⁵ / 230 180 |
| eventi n–p HP (attesa `jeff33`) | 1.0437·10⁶ / 1 043 882 | 2.245·10⁵ / 224 776 | 2.191·10⁵ / 219 093 |
| protoni contenuti confrontati con PSTAR | 1.07·10⁶ / 1 070 042 | 2.26·10⁵ / 223 921 | 2.20·10⁵ / 218 026 |
| `protons_hadronic` | 5 100 ± 700 / 5 624 | 900 ± 300 / 1 186 | 1 500 ± 400 / 1 229 |
| `events_secondary_n_hadronic` | 1.38·10⁵ / 136 775 | 3.0·10⁴ / 27 571 | 2.8·10⁴ / 27 023 |
| `Edep_overflow` | 7 300 ± 850 / 6 622 | 900 ± 300 / **1 550** | < 300 / 1 |
| `Evis_overflow` | 0 / 0 | 0 / 0 | 0 / 0 |
| `np_first_without_proton` (n–d, D17) | circa 120 / 108 | circa 26 / 21 | circa 25 / 23 |
| `protons_nuclear_recoil` (D22) | circa 16 / 14 | circa 3 / 2 | circa 4 / 5 |

Tutto entro circa 1σ della stima, tranne `Edep_overflow` nell'isotropa mono: 1 550 contro
900 ± 300, cioè +2.2σ. L'errore della stima viene dai 9 eventi del pilota, quindi lo scarto è
compatibile con una fluttuazione del pilota (c). `Edep_overflow` conta gli eventi con deposito
oltre il fondo scala dell'istogramma e non tocca nessuna soglia.

## Ordine dei fatti

1. **D27**, config v3, schema 4 (`b9bbb5d`), `checks_version` 10 (`a15d899`).
2. **Piloti** da 10⁵ eventi a 1 thread (`pilot3_*_1e5_1t`), build `release` di `a15d899`. Non
   sono run di validazione e i loro test non entrano qui.
3. **Previsioni v1** (`8429cc3`) e **v2** (`ff6f837`), prima dei run. La v2 aggiunge al repo il
   riferimento della luce per bin, e le soglie non cambiano.
4. **Produzione `c1`** (2026-11-22): i tre run sopra, build di `ded6c88`.

Nessuna soglia ha fallito, quindi la procedura della §9 delle previsioni non è servita. Non ci
sono run cancellati.

## Punti aperti

- **n–C non elastico lungo +X a +3.73σ** (§2): da studiare sullo stesso run, con profondità e
  vicinanza alle facce dei vertici, prima che l'analisi usi le frazioni.
- **Range contro PSTAR** (§4): la previsione (d) sul verso dello spostamento è fallita. Resta
  aperto anche il range sotto 3 MeV della fase 2.
- **Perdita di protoni dai bordi dipendente da E_p e direzione** (§5): fino al 17% sopra
  15.6 MeV nella sorgente di produzione. La selezione dei candidati a singolo n–p nella fase 4
  deve tenerne conto.
- **Massa HP di H-1 per l'analisi** (D14, D24): oggi l'AWR è solo in `checks.json`
  (`hp_elastic_H1.awr`), non in `run.json`. Va deciso come arriva all'analisi nella fase 4.
- **Offset HP di H-1** (D14), **bordo a 20 MeV** (D16), **eccesso di `ionIoni` sui rinculi C12**
  (D14) e **C-13** (D15): invariati, da decidere prima dell'analisi.
