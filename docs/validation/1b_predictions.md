# Previsioni della fase 1b — scritte prima dei run

Versione 1 (2026-11-19). Calcolate dal codice del commit `c95eb94` **senza simulare nulla**:
sezioni d'urto e coefficienti di Legendre delle valutazioni (D13), densità del materiale Geant4,
massa HP di H-1 da G4NDL (D14), carbonio composto come HP (D15), rampa HP/CHIPS (D16). Questo
file e i JSON in `predictions/` vanno committati **prima** di lanciare i run; i risultati
(`1b_results.md`) citeranno l'hash di quel commit. Dopo i run le soglie non si toccano.

Etichette: **(a)** fatto stabilito (sorgenti, dati, costruzione del codice); **(b)** stima con
ipotesi esplicite; **(c)** verificato con Monte Carlo; **(d)** ipotesi da mettere alla prova.
Qui non c'è nulla di (c): lo diventerà ciò che i run confermano.

## Run previsti

Tre run da 10⁷ eventi, seme del progetto `20260907`, cubo da 100 mm, fascio lungo z
(L = 100 mm), `QGSP_BIC_HP`, Geant4 11.4.1, G4NDL 4.7.1. Prima, un pilota da 10⁵ eventi
(mono 14 MeV) per le prestazioni.

| Run | Config | Ruolo |
|---|---|---|
| mono 14 MeV | `config/mono_14MeV.json` | riferimento, tutto HP |
| piatto 10–20 MeV | `config/flat_10_20MeV.json` | sorgente di produzione |
| bordo 19.5–20 MeV | `config/flat_19p5_20MeV.json` | rampa HP/CHIPS (D16) |

Riproduzione (≈ 4 s e ≈ 63 MB a run, nessuna scrittura su disco):

```bash
for c in mono_14MeV flat_10_20MeV flat_19p5_20MeV; do
  build/debug/sim/riptide_sim_predict --config config/$c.json --events 10000000 \
    --bins-u 50 --bins-z 50 > docs/validation/predictions/${c}_1e7.json
done
```

I JSON contengono tutto (conteggi, spettro u, errori, separazione fra librerie); qui sotto
solo i numeri che decidono l'esito. `riptide_sim_check` va lanciato con gli stessi
`--bins-u 50 --bins-z 50`.

## 1. Oracoli deterministici (soglie esatte)

Valgono per tutti e tre i run.

| Verifica | Soglia | Etichetta |
|---|---|---|
| `counters.must_be_zero.total` | = 0 | (d) — solo un bug li rende ≠ 0 |
| `Ep_residuals.max_abs` (eventi HP) | < 10⁻⁹ · E_n,max (come lo smoke test) | (d); smoke a 14 MeV: 1.24·10⁻¹⁴ MeV (a) |
| `model_mix.Ep_residuals.max_abs` (non HP) | < 10⁻¹¹ · E_n,max | (a) — vale per costruzione della classificazione |
| `Ep_physical_offset.relative_min/max` | dentro l'intervallo previsto, tabella sotto | (b) |
| `hp_elastic_H1.awr` | 0.999167 | (a) |
| `model_mix.observed_non_hp` nel run mono | = 0 | (a) — 14 MeV è sotto la finestra |

Intervallo previsto dell'offset relativo HP − fisico (θ_p in [0, π/2), su tutta la sorgente):

| Run | min | max |
|---|---|---|
| mono 14 MeV | 5.3658·10⁻⁴ | 5.4861·10⁻⁴ |
| piatto 10–20 MeV | 5.3322·10⁻⁴ | 5.5032·10⁻⁴ |
| bordo 19.5–20 MeV | 5.3322·10⁻⁴ | 5.5032·10⁻⁴ |

(b): entrambi gli estremi cadono a E_n,max, quindi piatto e bordo hanno lo stesso intervallo, più
largo di quello a 14 MeV. I valori osservati devono cadervi dentro, con margine di arrotondamento
relativo 10⁻⁹. Unica eccezione ammessa: un evento ambiguo (`expected_rare.np_model_ambiguous` > 0),
che abbassa `relative_min` (D16; attesi circa 0.06 nel run di bordo, 0 negli altri, (b)).

## 2. Test statistici contro `jeff33`

Il confronto con `jeff33` verifica l'implementazione: G4NDL 4.7.1 deriva da JEFF-3.3 (D13).
Ipotesi (d): la simulazione campiona esattamente le attese di `jeff33`.

Piano dei test (Bonferroni, α di famiglia 0.01 per run; i pull di Legendre sono correlati e la
correzione resta conservativa):

| Run | χ² | pull | α per test | soglia \|pull\| |
|---|---|---|---|---|
| mono 14 MeV | 3 (u, z, frazioni) | 6 (a₁…a₆) | 1.111·10⁻³ | 3.261 |
| piatto 10–20 MeV | 4 (+ `model_mix`) | 7 (+ pull `model_mix`) | 9.09·10⁻⁴ | 3.317 |
| bordo 19.5–20 MeV | 4 (+ `model_mix`) | 7 (+ pull `model_mix`) | 9.09·10⁻⁴ | 3.317 |

**Esito:** un χ² passa se `p_value` > α per test, un pull se |pull| < soglia. Un run passa se
passano tutti i suoi test e gli oracoli della §1. Con un sistema corretto la probabilità di
almeno un falso allarme sui tre run è ≈ 3% (b).

### Attese principali (jeff33, 10⁷ eventi) — (b)

Ipotesi: prima interazione da sezioni d'urto `jeff33` lungo L = 100 mm, nessun contributo del
deuterio (§4), quadratura in E_n con errore trascurabile (≤ 0.0026σ sulle frazioni, misurato
raddoppiando i nodi).

| | mono 14 | piatto 10–20 | bordo 19.5–20 |
|---|---|---|---|
| P(n–p) | 0.230406 | 0.219082 | 0.164717 |
| P(n–C el.) | 0.247378 | 0.267195 | 0.298039 |
| P(n–C non el.) | 0.146509 | 0.150003 | 0.157614 |
| P(altro) | 3.10·10⁻⁵ | 5.00·10⁻⁵ | 7.76·10⁻⁵ |
| P(nessuna) | 0.375675 | 0.363670 | 0.379552 |
| σ relativo P(n–p) | 5.8·10⁻⁴ | 6.0·10⁻⁴ | 7.1·10⁻⁴ |
| σ relativo P(n–C non el.) | 7.6·10⁻⁴ | 7.5·10⁻⁴ | 7.3·10⁻⁴ |
| eventi n–p HP | 2.304·10⁶ | 2.150·10⁶ | 8.27·10⁵ |
| eventi n–p non HP | 0 | 41 008 ± 117 | 820 158 ± 524 |
| quota non HP | 0 | 1.87% | 49.8% |
| conteggio u minimo per bin | 44 380 | 41 579 | 16 094 |

Momenti di Legendre ⟨P_l(μ)⟩ = a_l sugli eventi HP, con l'errore atteso a 10⁷ eventi:

| l | mono 14: a_l ± σ | piatto 10–20: a_l ± σ | bordo: a_l ± σ |
|---|---|---|---|
| 1 | −1.4798·10⁻² ± 3.8·10⁻⁴ | −1.4513·10⁻² ± 3.9·10⁻⁴ | −1.6091·10⁻² ± 6.4·10⁻⁴ |
| 2 | 1.893·10⁻³ ± 3.0·10⁻⁴ | 2.542·10⁻³ ± 3.1·10⁻⁴ | 7.123·10⁻³ ± 4.9·10⁻⁴ |
| 3 | −6.61·10⁻⁴ ± 2.5·10⁻⁴ | −7.40·10⁻⁴ ± 2.6·10⁻⁴ | −1.425·10⁻³ ± 4.2·10⁻⁴ |
| 4 | 2.52·10⁻⁴ ± 2.2·10⁻⁴ | 3.33·10⁻⁴ ± 2.3·10⁻⁴ | 8.68·10⁻⁴ ± 3.7·10⁻⁴ |
| 5 | −1.0·10⁻⁵ ± 2.0·10⁻⁴ | −1.6·10⁻⁵ ± 2.1·10⁻⁴ | −5.3·10⁻⁵ ± 3.3·10⁻⁴ |
| 6 | 1.4·10⁻⁶ ± 1.8·10⁻⁴ | 2.4·10⁻⁶ ± 1.9·10⁻⁴ | 8.1·10⁻⁶ ± 3.1·10⁻⁴ |

Conseguenza (b): a 10⁷ eventi solo a₁ è misurato con significatività (≈ 39σ a 14 MeV); a₂ e
a₃ sono a 6.4σ e 2.7σ, a₄–a₆ sono compatibili con zero. Il test sui pull resta comunque valido:
verifica che il campione non si allontani dalle attese.

Il pull di `model_mix` usa la varianza condizionata Σ p(1 − p) agli E_n campionati, che non
dipende dalla libreria; ±117 e ±524 sopra sono invece gli errori incondizionati, più larghi (b).

## 3. Confronto con `endfb80` — attese, non soglie

`endfb80` non è un criterio di accettazione: per H-1 è una revisione della stessa analisi
(D13), per il carbonio una valutazione diversa. Si riportano le differenze attese, in unità
dell'errore statistico del run (b):

| Run | max \|Δ frazioni\| | max \|Δa_l\| | non centralità del χ² di u |
|---|---|---|---|
| mono 14 MeV | 17.1σ (n–C non el.) | 1.35σ (a₃) | 3.2 |
| piatto 10–20 MeV | 15.8σ (n–C non el.) | 1.43σ (a₃) | 3.7 |
| bordo 19.5–20 MeV | 14.5σ (n–C non el.) | 1.63σ (a₃) | 4.8 |

Previsioni (b):

- il χ² delle frazioni contro `endfb80` **fallirà** nettamente: n–C non elastico −17σ, n–C
  elastico +11σ a 14 MeV. È la differenza fra le valutazioni del carbonio, non un errore;
- su n–p le due librerie restano entro ~2σ (frazione n–p: +1.5σ a 14 MeV, +0.9σ sul piatto,
  −2.2σ al bordo): a 10⁷ eventi la simulazione **non distingue** le due revisioni di H-1;
- lo spettro u contro `endfb80` dà un eccesso di χ² di 3–5 su ≈ 50 bin: non
  rilevabile;
- con a_l: a₁ vale −1.4447·10⁻² (endfb80) contro −1.4798·10⁻² (jeff33) a 14 MeV, cioè 0.9σ.

Questo quantifica il limite già scritto in D13: la fase 1b verifica l'implementazione, non i
dati n–p.

## 4. Residui noti non modellati

Sono effetti conosciuti, che le attese `jeff33` non includono. Entrano nel χ² come piccola
non centralità e non giustificano un fallimento.

| Effetto | Grandezza | Rispetto all'errore a 10⁷ | Etichetta |
|---|---|---|---|
| C-13: inelastico di G4NDL ≠ MT3 − MT102 di JEFF-3.3 (D15) | ≈ 2.9·10⁻⁴ relativo sul non elastico del carbonio a 14 MeV | ≈ 0.4σ (σ_rel 7.6·10⁻⁴) | (b) |
| Scattering n–d: nella verità è n–p (Z = 1), nelle attese no | n_H2/n_H1 = 1.150·10⁻⁴ → ≈ 265 / 247 / 95 eventi su scala n–p (mono / piatto / bordo) | ≈ 0.2σ sul conteggio n–p (σ ≈ 1330 eventi a 14 MeV) | (b): l'ordine di grandezza usa σ_np al posto di σ_nd, che non è valutata qui |
| Massa HP di H-1 (D14) | E_p più alto di ≈ 5.4·10⁻⁴ relativo | modellato nelle verifiche; aperto per la fase 2 | (a) |

Gli eventi n–d con deutone di rinculo non hanno protone: non entrano nello spettro u né nei
residui, e sono contati in `expected_rare` (a).

## 5. Prestazioni attese (dal pilota, 2026-11-19)

Pilota mono 14 MeV, build `release` del commit `6de0124`, sotto `scripts/run_capped.sh -m 5G`,
RSS da `/usr/bin/time -v`, tempo da `run.json` (`elapsed_s`, inizializzazione compresa).

| Grandezza | 10⁵ eventi, 1 thread | 4·10⁵ eventi, 8 thread |
|---|---|---|
| tempo | 12.8 s | 9.1 s |
| eventi/s | 7.8·10³ | 4.4·10⁴ (5.5·10³ per thread) |
| RSS massimo | 89 MB | 139 MB |
| byte/evento su disco | 54.5 | 55.1 |

`riptide_sim_check` sul pilota da 4·10⁵: 1.8 s, RSS 65 MB.

Stime per un run da 10⁷ con 8 thread (b): circa 4 minuti, RSS sotto 0.2 GB, circa 0.55 GB
su disco; la verifica, circa 45 s. Il run piatto e quello di bordo hanno E_n più alto e un po'
più di interazioni sul carbonio: il tempo può crescere, non l'ordine di grandezza. Si usano
8 thread su 20 core, con tetto di 5G a run, un run alla volta.

Sul pilota l'attesa n–p di `riptide_sim_check` (somma evento per evento) coincide con la
previsione scalata a 4·10⁵ eventi: 92 162.5678 contro 92 162.5678 (c). Il pilota non è un run
di validazione e i suoi test statistici non entrano nei risultati.

## 6. Se una soglia fallisce

1. **Le soglie non si cambiano** dopo aver visto i dati, e non si ripete il run con un altro
   seme per farlo passare.
2. Si cerca la causa sul run esistente: quale bin, categoria o grado; con quale segno; se
   dipende da E_n (piatto contro bordo contro mono). Gli oracoli della §1 sono esatti: un
   fallimento lì è un bug finché non si dimostra il contrario.
3. Se la causa è una convenzione di Geant4 non ancora modellata (come D14, D15, D16), la si
   verifica sui sorgenti, si scrive una decisione, si corregge il codice delle attese
   (`checks_version` + 1) e si rilancia `riptide_sim_check` **sugli stessi eventi**. Nei
   risultati si riportano entrambi gli esiti, con l'ordine dei fatti.
4. Un solo fallimento statistico senza causa trovata si riporta come tale, con la sua
   probabilità di falso allarme (≈ 1% a run), e non si nasconde.
