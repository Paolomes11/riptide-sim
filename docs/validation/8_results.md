# Risultati della fase 8 — Caso C

Data: 2026-11-26. Previsioni e soglie: `docs/validation/8_predictions.md`, versione 1 in
`58bb185`, committata **prima** del codice del Caso C e dei run; versione 2 (§2.1, diagnosi
8c) in `daa2491`, committata prima del run della 8c. Dopo i run le soglie non sono
state toccate. Versione 3 (§4.1–4.3, 8b) in `3818815`, committata prima dei run della 8b.
Questo documento copre la 8a (port contro toy), la diagnosi 8c e la 8b (Geant4, §7, del
2026-11-27).

**Esito 8b (§7).** Le 12 chiamate passano il criterio su +X (3 fallimenti contro un quantile
di 6) e sull'obliqua (7 contro 8). Su +Z **scatta** (10 contro 9): quattro voci sono a bassa
potenza, e sei sono rapporti rms con un gemello ideale di rms bassa, che col gemello con
accettanza passano. Si apre la 8c su +Z con una decisione a parte. In Geant4 il bias di μ_E a
N = 1000 è negativo su +Z e sull'obliqua (−0.05 / −0.04 MeV, −0.03 / −0.02 con l'accettanza),
entro soglia.

**Esito 8a.**
- Il primo run (`8a`) ha trovato un difetto di quadratura nello stadio 2: la riga **bordo**
  falliva 20 voci su 24. È corretto in D49 (`813d154`) e la checklist è stata rilanciata per
  intero (`8a_w`). Il difetto è anche nel toy (§4).
- Con il modello corretto:
  - checklist **64/68** voci contate, nessuna obbligatoria fallita;
  - robustezza al prior **5/5**;
  - bordo **23/24**;
  - sborda: calibrazione rotta come previsto, coverage 68% di μ_E 0.50 a N = 1000 (§3).
- I 4 fallimenti della checklist sono al bordo superiore dei 2–4 falsi allarmi attesi, ma
  riguardano tutti μ_E. Tre su quattro sono la stessa riga, N = 150: la larghezza del pull e la
  coverage a 0.68 e 0.90, misurate sulle stesse esperienze. Per la regola della §2 delle
  previsioni ("fallimenti concentrati su un parametro o un N") **si apre la diagnosi 8c** (§5).

**Esito 8c: non conclusiva** per il criterio della §2.1 delle previsioni. Su 200 esperienze
nuove a N = 150 fallisce 1 voce su 17, la larghezza del pull di μ_E, ma dall'altra parte: 0.913
contro 1.070 della 8a. Le coverage di μ_E tornano nominali (0.700 / 0.935 / 0.965). Il
difetto previsto dal criterio (pull di nuovo sopra 1.07 e coverage a 0.68 di nuovo bassa) non
c'è, e le 400 esperienze insieme sono calibrate (riferimento, §5). Non si apre una diagnosi del
modello.

Etichette: **(a)** fatto stabilito; **(b)** stima con ipotesi esplicite; **(c)** verificato con
Monte Carlo; **(d)** ipotesi da mettere alla prova.

## 1. Run

Build `release` a `813d154` (`git describe` v0.9.0-12-g813d154, `project_version` 0.9.0).

```
scripts/run_capped.sh -m 5G -- /usr/bin/time -v \
    build/release/analysis/apps/riptide_case_c_calibrate --calibration-id 8a_w --threads 4
```

- Impostazioni del record:
  - seme `20260907`, σ_Ep = 0.10 MeV, σ_θ = 0.08 rad;
  - dominio di E_n [10, 20] MeV su 500 punti;
  - 3000 direzioni sulla sfera, calotta di 8000 punti (Δlog = 10);
  - tabella θ / θ_obs di 400 × 801 punti;
  - griglia degli iperparametri globale 60 × 60 (σ_E in 0.01–50 MeV, geometrica), finestra a
    Δlog = 10 e griglia fine 40 × 40.
- Output sull'SSD: `stats/8a_w/case_c.json`. Il run prima della correzione è in
  `stats/8a/case_c.json` (build `9e2227f`), il pilot in `stats/8a_pilot/`.
- Le soglie sono applicate dall'app, che scrive l'esito accanto a ogni numero.
- Tempo: 23 min 15 s a parete con 4 thread, RSS massimo 100 MB (a).

## 2. Port contro toy (§2 delle previsioni)

Stime: μ_E = media a posteriori, log σ_E = mediana, Ω_n = MAP. Il toy è la tabella di
riferimento delle previsioni. A N = 1000 l'esito si legge sulla ripetizione (M = 240).

| N | M | bias μ_E (MeV) | bias log σ_E | rms μ_E (std media) | rms log σ_E (std media) | rms Ω_n (°) |
|---|---|---|---|---|---|---|
| 50 | 200 | +0.068 ± 0.033 | −0.043 ± 0.058 | 0.477 (0.782) | 0.826 (0.725) | 1.68 |
| 150 | 200 | −0.003 ± 0.024 | −0.005 ± 0.017 | 0.344 (0.306) | 0.237 (0.215) | 0.95 |
| 300 | 100 | +0.010 ± 0.017 | −0.017 ± 0.013 | 0.173 (0.180) | 0.134 (0.130) | 0.64 |
| 1000 | 40 | +0.003 ± 0.013 | −0.003 ± 0.011 | 0.082 (0.087) | 0.067 (0.064) | 0.31 |
| 1000 | 240 | +0.005 ± 0.006 | −0.006 ± 0.004 | 0.090 (0.088) | 0.067 (0.064) | 0.36 |

| N | M | pull μ_E media / larghezza | pull log σ_E media / larghezza | pull rms Ω_n | coverage μ_E | coverage log σ_E | coverage Ω_n |
|---|---|---|---|---|---|---|---|
| 50 | 200 | +0.09 / **0.89** | +0.06 / 1.00 | 0.98 | 0.705 / 0.880 / 0.955 | 0.535 / 0.825 / 0.895 | 0.680 / 0.925 / 0.970 |
| 150 | 200 | +0.06 / **1.070** | +0.02 / 1.00 | 1.01 | **0.590** / **0.855** / 0.935 | 0.645 / 0.890 / 0.930 | 0.680 / 0.885 / 0.935 |
| 300 | 100 | +0.07 / 1.04 | −0.11 / 0.97 | 0.96 | 0.720 / 0.870 / 0.920 | 0.700 / 0.860 / 0.940 | 0.690 / 0.900 / 0.950 |
| 1000 | 40 | +0.07 / 0.94 | −0.03 / 1.07 | 0.86 | 0.675 / 0.925 / 0.975 | 0.675 / 0.900 / 0.950 | 0.775 / 0.975 / 1.000 |
| 1000 | 240 | +0.05 / 0.98 | −0.10 / 1.05 | 1.01 | 0.671 / 0.917 / 0.950 | 0.637 / 0.887 / 0.950 | 0.675 / 0.887 / 0.942 |

In grassetto le voci fallite.

- **Bias:** tutti entro 0.25 × std media (c).
- **rms C++ / toy:** tutti entro 1 ± 2√(1/M). μ_E da 0.83 a 0.98, log σ_E da 0.90 a 1.19,
  Ω_n da 0.91 a 1.07 (c).
- **`truth_in_cap`:** 1 in tutte le esperienze (a).
- **Pendenze log-log su 150–1000, con la ripetizione:** μ_E −0.69 (toy −0.71), log σ_E −0.65
  (toy −0.73), Ω_n −0.498 (soglia −0.5 ± 0.07). Passano tutte (c).
- **Pendenza di Ω_n sulle sole righe del toy (con M = 40):** −0.58, fuori soglia di 0.011. È solo di
  riferimento: a N = 1000 conta la ripetizione, e la riga M = 40 ha pull rms 0.86, cioè la
  stessa fluttuazione della riga da 40 esperienze.
- **Eccezione dichiarata, coverage di log σ_E a N = 50:** 0.535 / 0.825 / 0.895, contro il toy
  0.61 / 0.81 / 0.87. Fuori dai 2σ del toy solo a 0.68 (−2.3σ). Si riporta senza esito.
- **Pull medio comune di μ_E:** +0.068 ± 0.037 su tutte le righe, contro +0.053 ± 0.024 nel
  toy (d): coincide entro gli errori.

### 2.1 I quattro fallimenti

| Voce | Valore | Soglia | Toy |
|---|---|---|---|
| N = 50, larghezza del pull di μ_E | 0.893 | 1 ± 0.07 | 0.88 |
| N = 150, larghezza del pull di μ_E | 1.0702 | 1 ± 0.07 | 1.00 |
| N = 150, coverage μ_E 0.68 | 0.590 | 0.614–0.746 | 0.68 |
| N = 150, coverage μ_E 0.90 | 0.855 | 0.858–0.942 | 0.90 |

- **N = 50:** il toy ha 0.88, che fallirebbe la stessa soglia: il C++ riproduce il toy (c).
  A N = 50 la std a posteriori di μ_E (0.78 MeV) è più larga dell'rms (0.48): domina il prior
  uniforme su [10, 20] MeV, e il pull è stretto.
- **N = 150:** tre voci sulle stesse 200 esperienze, quindi un solo fatto.
  - La coverage a 0.68 è −2.7σ binomiali; quella a 0.90 è −2.1σ.
  - La larghezza del pull supera la soglia di 0.0002. Con M = 200 la soglia ± 0.07 vale circa
    1.4σ statistici.
  - L'rms di μ_E coincide con il toy (0.344 contro 0.350), ma la std a posteriori media è
    0.306 contro 0.328: posterior circa il 7% più stretti, in accordo con la larghezza del pull.
  - Sulle stesse esperienze, la riga di robustezza con il prior uniforme in σ_E dà coverage di
    μ_E 0.61 / 0.86 / 0.94: anche lì basso a 0.68.
  - Ipotesi (d): fluttuazione delle 200 esperienze, come la 7b → 7c. In alternativa, un difetto
    a N intermedio che non si vede a N = 300 e 1000. Non è verificata: è la 8c (§5).
- La correzione D49 non c'entra. Le stesse voci a N = 150 erano già fuori prima (0.590 / 0.850
  a `9e2227f`). La larghezza del pull (1.069 prima, 1.070 dopo) supera la soglia solo dopo, ed
  è l'unica differenza nel conteggio fra 3 e 4 (a).

## 3. Robustezza al prior e troncatura (§3 delle previsioni)

**Prior uniforme in σ_E, N = 150, M = 200** (stesse esperienze della baseline): 5/5.

| Voce | Valore | Soglia |
|---|---|---|
| Spostamento massimo di μ_E | 0.40 σ | < 0.5 σ |
| Spostamento medio di log σ_E | +0.20 σ | +0.21 ± 0.10 σ |
| Coverage log σ_E | 0.625 / 0.860 / 0.920 | 2σ binomiali dal nominale |

Prima della correzione la coverage di log σ_E con il prior alternativo falliva a 0.90 e 0.95.
La differenza viene dai pesi, che cambiano la regione a σ_E grande, dove i nodi estremi pesano
(a).

**Bordo** (μ_E ~ U(17.5, 19), verità troncata): 23/24.

| N | M | pull μ_E media / larghezza | pull log σ_E larghezza | pull rms Ω_n | coverage μ_E | coverage log σ_E | coverage Ω_n |
|---|---|---|---|---|---|---|---|
| 300 | 100 | +0.42 / **0.81** | 0.97 | 0.97 | 0.69 / 0.93 / 0.97 | 0.71 / 0.87 / 0.96 | 0.69 / 0.92 / 0.95 |
| 1000 | 40 | +0.06 / 0.94 | 0.94 | 1.04 | 0.75 / 0.93 / 0.95 | 0.78 / 0.95 / 0.95 | 0.73 / 0.88 / 0.95 |

- A N = 300 la stima di μ_E è spostata verso il bordo: bias +0.20 ± 0.03 MeV, pull medio
  +0.42. Il pull è stretto (0.81, soglia 0.86–1.14), ma la coverage passa.
- Ipotesi (d): la gaussiana troncata vicino al bordo lascia μ_E e σ_E degeneri. Il posterior è
  asimmetrico e chiuso a μ_E = 20 dal prior. La media si sposta verso il bordo, e la std
  supera la dispersione delle stime. A N = 1000 l'effetto scompare (bias +0.04 ± 0.03, larghezza
  0.94).
- Non c'è un valore del toy per il bordo: la soglia è (b).

**Sborda** (verità non troncata, modello troncato), test di sensibilità senza esito.

| N | M | bias log σ_E | pull μ_E / log σ_E larghezza | coverage μ_E | coverage log σ_E |
|---|---|---|---|---|---|
| 300 | 100 | +0.066 ± 0.018 | 1.10 / 1.16 | 0.56 / 0.85 / 0.91 | 0.57 / 0.84 / 0.88 |
| 1000 | 40 | +0.049 ± 0.015 | 1.31 / 1.23 | 0.50 / 0.825 / 0.875 | 0.525 / 0.675 / 0.80 |

- Calibrazione rotta e peggiore a N = 1000, come previsto. La coverage 68% di μ_E a N = 1000
  vale 0.50, nell'intervallo 0.38–0.55 del toy (c).
- σ_E è sovrastimata: gli eventi oltre il dominio mancano, e il modello li attribuisce a una
  gaussiana più larga, troncata (d).

## 4. Il difetto del bordo e D49

Il prior π(E_n | μ_E, σ_E) era normalizzato con i trapezi (estremi a peso h/2), mentre la somma
su E_n dello stadio 2 era senza pesi. Con σ_E sotto il passo di E_n (0.01 < 0.02 MeV) tutta la
massa di π cade su un nodo. Su un estremo la densità è doppia, e ogni evento guadagna log 2.

**Prima della correzione** (`stats/8a`, build `9e2227f`), riga bordo, 20/24 voci fallite (a):

| N | larghezza pull μ_E | larghezza pull log σ_E | coverage μ_E | coverage log σ_E |
|---|---|---|---|---|
| 300 | 1.2·10⁹ | 277 | 0.36 / 0.46 / 0.49 | 0.32 / 0.43 / 0.49 |
| 1000 | 1.2·10²⁵ | 959 | 0.35 / 0.53 / 0.55 | 0.40 / 0.50 / 0.55 |

Le larghezze enormi e le coverage vicine a 0.5 indicano che una parte grande delle esperienze
collassava su (μ_E = 20, σ_E = 0.01), con std a posteriori quasi nulla. Il record non salva le
stime per esperienza, quindi la frazione esatta non è misurata.

- **Controllo sul toy (c):** le funzioni del toy, con N = 300, danno gli argmax di D49:
  - (20.00, 0.010) per le verità (18.5, 1.5) e (18.0, 1.0);
  - (19.15, 1.81) e (18.31, 1.56) con i pesi aggiunti;
  - nessuna differenza per una verità interna, (15.0, 1.5).

  Il toy ha lo stesso difetto; la sua checklist non ha una riga bordo.
- **Effetto sulle righe interne (a):** minimo. I numeri della §2 cambiano alla terza cifra; il
  conteggio passa da 3 a 4 fallimenti, per la larghezza del pull a N = 150 (§2.1).
- **Riproduzione del toy (b):** il test su un'esperienza fissata (50 eventi) usa le funzioni del
  toy con gli stessi pesi aggiunti al prior. Coincide con il C++ entro 1e-3 relativo e con lo
  stesso pixel per ω₀ e ω̂.

## 5. Diagnosi 8c: μ_E a N = 150 con altri semi (§2.1 delle previsioni)

**Run.** Build `release` a `035bb44` (`git describe` v0.9.0-15-g035bb44), che aggiunge
all'app l'opzione `--row-events` per lanciare una sola riga della checklist.

```
scripts/run_capped.sh -m 5G -- /usr/bin/time -f "RSS %M kB, %e s" \
  build/release/analysis/apps/riptide_case_c_calibrate \
  --calibration-id 8c_n150 --row-events 150 --seed 20260908
```

- La riga è N = 150, M = 200, con il prior baseline. Sezione, verità, risoluzioni e griglie sono
  quelle della 8a; cambia solo il seme.
- Output: `stats/8c_n150/case_c.json`, e il pilot (M = 20) in `stats/8c_pilot/`.
- Tempo 79 s con 4 thread, RSS massimo 46 MB (a).

**Voci della §2 a N = 150 (c):** fallisce 1 di 17, nessuna obbligatoria.

| Voce | 8a (`8a_w`) | 8c | Soglia (8c) | Toy |
|---|---|---|---|---|
| bias μ_E (MeV) | −0.003 ± 0.024 | +0.032 ± 0.022 | \|·\| < 0.080 | +0.029 ± 0.025 |
| bias log σ_E | −0.005 ± 0.017 | −0.024 ± 0.019 | \|·\| < 0.057 | +0.004 ± 0.019 |
| rms μ_E / toy | 0.98 | 0.87 | 0.86–1.14 | 1 |
| rms log σ_E / toy | 0.90 | 1.01 | 0.86–1.14 | 1 |
| rms Ω_n / toy | 1.07 | 1.07 | 0.86–1.14 | 1 |
| larghezza del pull di μ_E | **1.070** | **0.913** | 1 ± 0.07 | 1.00 |
| larghezza del pull di log σ_E | 1.00 | 0.97 | 1 ± 0.07 | 0.95 |
| pull rms Ω_n | 1.01 | 0.99 | 1 ± 0.07 | 0.93 |
| coverage μ_E 68/90/95 | **0.590** / **0.855** / 0.935 | 0.700 / 0.935 / 0.965 | ± 2σ binomiali | 0.68 / 0.90 / 0.95 |
| coverage log σ_E 68/90/95 | 0.645 / 0.890 / 0.930 | 0.650 / 0.860 / 0.920 | ± 2σ binomiali | 0.67 / 0.87 / 0.94 |
| coverage Ω_n 68/90/95 | 0.680 / 0.885 / 0.935 | 0.670 / 0.900 / 0.960 | ± 2σ binomiali | 0.73 / 0.95 / 0.98 |

Le soglie del bias dipendono dalla std media di ciascun run: nella 8a erano 0.077 e 0.054. La
coverage di log σ_E a 0.90 (0.860) e a 0.95 (0.920) della 8c
passa per meno di 0.002.

**Esito per il criterio della §2.1: non conclusiva.**
- *Fluttuazione* chiede che passino le tre voci di μ_E fallite nella 8a. Passano le due
  coverage, ma la larghezza del pull fallisce di nuovo, a 0.913 (−1.7σ statistici, con
  σ ≈ 1/√(2M) = 0.05).
- *Difetto possibile* chiede la larghezza sopra 1.07 **e** la coverage a 0.68 sotto banda. Non
  succede: la larghezza è sotto 1 e la coverage è 0.700.
- Quindi, come scritto nelle previsioni, la si riporta così e non si cambia seme.

**Di riferimento, senza esito (b).**
- Std a posteriori media di μ_E: 0.321 MeV, contro 0.306 nella 8a e 0.328 nel toy. La std
  più stretta della 8a non si ripete: è variata fra i due campioni di 200 esperienze.
- Le 400 esperienze della 8a e della 8c insieme, per μ_E:
  - larghezza del pull 0.995;
  - coverage 0.645 / 0.895 / 0.950, entro 2σ binomiali a M = 400 (0.633–0.727, 0.870–0.930,
    0.928–0.972);
  - rms 0.325 e std media 0.313 MeV.
- La larghezza del pull si è spostata da +1.4σ a −1.7σ attorno a 1. Un difetto del modello
  sposterebbe i due campioni dalla stessa parte; qui si spostano da parti opposte. Questo
  sostiene la fluttuazione (d), ma la 8c, letta con il criterio scritto, non la conferma.
- La soglia ± 0.07 della larghezza a M = 200 vale circa 1.4σ, quindi una riga corretta la
  supera con probabilità ≈ 16% (`8_predictions.md` §2.1). È la stessa soglia stretta che fa
  fallire la riga N = 50 nel toy (0.88).

**Che cosa resta.** I fallimenti della 8a restano registrati come sono. Nessuna diagnosi del
modello si apre: il criterio del difetto non scatta. Si passa alla 8b con le previsioni
versione 3.

## 6. Prestazioni (§5 delle previsioni)

`riptide_stats_bench --case c`, un'esperienza a N = 1000, thread singolo, build `release` a
`813d154`, migliore di 3 ripetizioni (a):

| Stadio | Tempo (s) |
|---|---|
| base per evento | 3.21 |
| stadio 1 (Caso B, prior piatto) | 1.85 |
| stadio 2 (globale + finestra + fine) | 1.43 |
| stadio 1′ (Caso B, prior predittivo) | 1.88 |
| totale | 8.37 |

- Sopra l'obiettivo di 5 s.
- Lo stadio 2 è sceso da 10.8 a 1.4 s con il prodotto scalare in spazio lineare (`9e2227f`).
  Il resto è codice del Caso B (`energyLogPosterior`), che la 8a non tocca. I pesi di D49 non
  cambiano i tempi.
- La checklist completa costa 23 min con 4 thread, contro ~8 min del toy con 3 processi.

## 7. 8b, Geant4 (§4 delle previsioni, versione 3)

**Run.** Build `release` a `3818815` (`git describe` v0.9.0-21-g3818815), le 12 chiamate della
§4.2 delle previsioni in sequenza, 4 thread, semi di default, ognuna sotto
`scripts/run_capped.sh -m 5G -- /usr/bin/time -v`, il 2026-11-27 dalle 10:07 alle 12:23 (a).
- Output: `stats/8b_{x,z,obl}{,_acc}{,_toy}/case_c.json`.
- Tempi a parete per chiamata: 8 min 23 s – 8 min 57 s su +Z e sull'obliqua, 16 min 29 s –
  16 min 46 s su +X, per gemello e Geant4 (stima: ~8 e ~15 min). RSS massimo 1.41 GB per
  Geant4, 0.13 GB per il toy.
- Le M eseguite sono quelle del dry run: 200 / 200 / 100 / 40 sulle righe, ripetizione 240 su
  +X, 62 su +Z, 73 sull'obliqua (a).

**Conteggio per il criterio della §4.3 (a).** Voci contate fallite, senza la larghezza del pull
di μ_E a N = 50 (fallimento atteso dal toy, fuori conteggio):

| | +X | +Z | obliqua |
|---|---|---|---|
| ideale (`8b_d`) | 2 / 61 | 8 / 61 | 4 / 60 |
| accettanza (`8b_d_acc`) | 1 / 61 | 2 / 61 | 3 / 60 |
| totale | **3** | **10** | **7** |
| quantile 99% di Poisson | 6 | 9 | 8 |
| attesi (b) | 2.2 | 4.0 | 3.5 |
| esito | passa | **scatta** | passa |

- **Larghezza del pull di μ_E a N = 50 (c):** fallisce in tutte e sei le chiamate Geant4, come
  previsto: 0.880 / 0.873 (+X), 0.870 / 0.871 (+Z), 0.894 / 0.892 (obliqua). I gemelli danno
  0.81–0.94.
- **Per la regola scritta, su +Z si apre la 8c**, con una decisione a parte. +X e l'obliqua si
  leggono voce per voce.

### 7.1 Le voci fallite

| Direzione | Chiamata | Voce | Valore | Soglia | Bassa potenza |
|---|---|---|---|---|---|
| +X | ideale | N = 150, rapporto rms log σ_E | 1.203 | ≤ 1.15 | |
| +X | ideale | N = 300, larghezza pull μ_E | 0.834 | 1 ± 0.1 | sì |
| +X | accettanza | N = 300, larghezza pull μ_E | 0.834 | 1 ± 0.1 | sì |
| +Z | ideale | N = 50, rapporto rms log σ_E | 1.246 | ≤ 1.15 | |
| +Z | ideale | N = 150, rapporto rms log σ_E | 1.180 | ≤ 1.15 | |
| +Z | ideale | N = 300, rapporto rms μ_E | 1.490 | ≤ 1.15 | |
| +Z | ideale | N = 300, rapporto rms log σ_E | 1.309 | ≤ 1.15 | |
| +Z | ideale | ripetizione, larghezza pull μ_E | 0.882 | 1 ± 0.1 | sì |
| +Z | ideale | ripetizione, rapporto rms μ_E | 1.281 | ≤ 1.15 | sì |
| +Z | ideale | ripetizione, rapporto rms log σ_E | 1.265 | ≤ 1.15 | sì |
| +Z | ideale | ripetizione, bias Ω_n componente 2 | 1.25·10⁻³ rad | ≤ 2 errori (1.02·10⁻³) | |
| +Z | accettanza | ripetizione, larghezza pull μ_E | 0.868 | 1 ± 0.1 | sì |
| +Z | accettanza | ripetizione, bias Ω_n componente 2 | 1.22·10⁻³ rad | ≤ 2 errori (1.03·10⁻³) | |
| obliqua | ideale | N = 150, larghezza pull log σ_E | 0.857 | 1 ± 0.1 | |
| obliqua | ideale | N = 150, rapporto rms log σ_E | 1.274 | ≤ 1.15 | |
| obliqua | ideale | ripetizione, larghezza pull μ_E | 0.866 | 1 ± 0.1 | sì |
| obliqua | ideale | ripetizione, pull rms Ω_n | 0.890 | 1 ± 0.1 | |
| obliqua | accettanza | N = 150, larghezza pull log σ_E | 0.857 | 1 ± 0.1 | |
| obliqua | accettanza | ripetizione, larghezza pull μ_E | 0.857 | 1 ± 0.1 | sì |
| obliqua | accettanza | ripetizione, pull rms Ω_n | 0.872 | 1 ± 0.1 | |

La componente 2 del bias di Ω_n è il secondo vettore tangente (indice 1, «component 1» nel
`tally` del record; errore 5.1·10⁻⁴ rad). Le voci obbligatorie non ci sono: la 8b ha solo voci
contate e di riferimento.

### 7.2 Lettura nell'ordine della checklist

**Bias.** Tutte le voci passano: μ_E e log σ_E entro 0.5 std, Ω_n sull'obliqua entro
0.5 σ_c, tranne la componente di Ω_n su +Z (sopra).
- **μ_E, negativo in Geant4 (a):** alla ripetizione a N = 1000

  | | +X | +Z | obliqua |
  |---|---|---|---|
  | ideale | −0.013 ± 0.007 | −0.050 ± 0.014 | −0.042 ± 0.011 |
  | accettanza | +0.007 ± 0.008 | −0.032 ± 0.014 | −0.023 ± 0.011 |
  | gemello ideale | −0.001 ± 0.006 | −0.003 ± 0.012 | +0.002 ± 0.012 |

  Il segno è negativo anche a N = 150 e 300 su +Z e sull'obliqua (−0.03…−0.08 MeV), mentre i
  gemelli sono positivi o nulli. Su +Z ideale il bias è −0.0504 contro una soglia di 0.0505:
  passa per 10⁻⁴ MeV.
- **Lettura (d):** con l'accettanza il bias scende di circa 0.02 MeV su tutte e tre le
  direzioni, come nella fase 6 senza q (−0.023 MeV a 14 MeV). Resta un residuo di −0.02…−0.03
  MeV su +Z e sull'obliqua, non su +X. La candidata è la stessa parte della selezione che
  dipende dalla direzione (D45, punto aperto della 7c); non è verificata.
- **Ω_n su +Z:** la componente 2 è a +2.4σ in tutte e due le chiamate, sulle stesse
  esperienze: un solo fatto. Norma ≈ 1.2·10⁻³ rad contro σ_c ≈ 6·10⁻³ rad (b), circa 0.2 σ_c.

**Risoluzione: i rapporti rms.**
- Su +Z il rapporto è rms Geant4 / rms del gemello. L'rms Geant4 è lo stesso nelle due
  chiamate (stesse esperienze; l'accettanza nel modello lo cambia alla terza cifra), ma il
  gemello ideale e quello con accettanza hanno eventi diversi (a):

  | +Z, rms μ_E (MeV) | N = 150 | N = 300 | N = 1000, M = 40 | ripetizione |
  |---|---|---|---|---|
  | Geant4 ideale | 0.380 | 0.292 | 0.145 | 0.117 |
  | Geant4 accettanza | 0.385 | 0.293 | 0.146 | 0.118 |
  | gemello ideale | 0.356 | **0.196** | **0.078** | **0.091** |
  | gemello accettanza | 0.353 | 0.306 | 0.103 | 0.103 |

  I sei rapporti rms falliti su +Z vengono da un gemello ideale con rms bassa: con il gemello
  con accettanza gli stessi rapporti passano. μ_E: 0.957 a N = 300, 1.138 alla ripetizione.
  log σ_E: 0.997 a N = 50, 1.118 a N = 150, 1.023 a N = 300, 0.957 alla ripetizione. Per
  log σ_E l'rms del gemello ideale è 0.810 / 0.264 / 0.130 / 0.058, contro 1.016 / 0.281 /
  0.168 / 0.079 per il gemello con accettanza.
- **Lettura (d):** i due gemelli hanno le stesse verità, ma con l'accettanza il toy ha un altro
  campione e un'altra verosimiglianza, quindi rms diverse non sono un difetto di per sé. A
  N = 300 la differenza fra i due gemelli è grande rispetto alla dispersione attesa
  (sd(ln r) ≈ 0.1 a M = 100); la causa non è verificata.
- Fuori dalla ripetizione e dalle righe di riferimento, su +X e sull'obliqua i rapporti
  rms di μ_E stanno fra 0.80 e 1.10.
- Il rapporto delle std a posteriori di Ω_n (G4 / gemello) è 0.975–1.004 ovunque (a).

**Pull.**
- La larghezza del pull di μ_E alla ripetizione è bassa su +Z e sull'obliqua (0.857–0.882,
  bassa potenza) e alta su +X (1.055–1.059, entro soglia).
- Sull'obliqua, N = 150, la larghezza del pull di log σ_E è 0.857 in tutte e due le chiamate.
- Il pull rms di Ω_n alla ripetizione sull'obliqua è 0.89 / 0.87, contro 1.03–1.11 dei
  gemelli; su +X e +Z è 0.96–1.03.

**Coverage.** Tutte le voci contate passano in tutte e sei le chiamate.

**Contrazione 1/√N.** Non è una voce della 8b. L'rms Geant4 di μ_E scende da 0.53 a 0.11–0.12
MeV fra N = 50 e la ripetizione, come i gemelli (a).

**Robustezza al prior.** Misurata nella 8a (§3), non ripetuta.

### 7.3 Esito

**+X e obliqua: passano** il criterio (3 contro 6, 7 contro 8), con le voci della §7.1
registrate come sono.

**+Z: il criterio scatta** (10 contro 9), quindi si apre la 8c con una decisione e previsioni
scritte prima dei run. Che cosa si sa già (a):
- quattro dei dieci fallimenti sono voci a bassa potenza;
- i due bias di Ω_n sono un solo fatto;
- sei fallimenti dell'ideale sono rapporti rms con un gemello di rms bassa, che con il gemello
  con accettanza passano. Due di questi (alla ripetizione) sono anche fra le quattro voci a
  bassa potenza.
- In sintesi: 4 rapporti rms a N ≤ 300 legati al gemello, 2 rapporti rms e 2 larghezze del
  pull alla ripetizione a bassa potenza, 2 bias di Ω_n che sono un solo fatto.

Il bias negativo di μ_E su +Z e sull'obliqua (§7.2) passa le soglie, ma va guardato nella 8c.
