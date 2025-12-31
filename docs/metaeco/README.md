# Meta-ecology tool concordance (portfolio-safe)

Markdown conversion of the original LaTeX notes. Paths/identifiers have been generalized.

## What this pipeline does (in one paragraph)

Given a merged abundance table from two profiling tools (e.g., , ) and optionally per-sample read depths \((N)\), this pipeline (i) compares tools in a compositional space (CLR/ALR/PLR) and reports Pearson/Spearman plus the Aitchison distance, (ii) estimates unbiased factorial moments to recover mean proportions \(_i\) for each taxon, (iii) fits abundance fluctuation distributions (AFDs) using a Negative Binomial (Gamma–Poisson) model and optionally Poisson–lognormal, (iv) performs goodness-of-fit and posterior predictive checks, and (v) derives detection probabilities, \(N^{*}(q)\) sequencing-depth targets, and richness forecasts \(R(N)\) with slopes \(R'(N)\). Utilities include metadata-driven sample selection, negatives/decontam, stratified diagnostics, Taylor’s law, and MAD lognormal checks.

## Compositional background

Let a sample be a composition \(=(x_1,,x_D)\) with \(x_i>0\) and \(_i x_i=1\).

\[
(x_i) \;=\;  x_i - {D}_{j=1}^D  x_j.
\]
Distances in CLR space correspond to Aitchison geometry on the simplex. In practice, zeros require replacement.

Choose a reference part \(r\) (or an average of a set of references) and set
\[
_r(x_i) \;=\;  x_i -  x_r.
\]
Robust if the reference is stable.

For taxa \(i,j\) present in both tools:
\[
{x_j} \;=\;  x_i -  x_j.
\]
Operates on pairs; robust to sparsity; ignores taxa not present in both vectors.

 For two compositions (after a common subcomposition/zero handling),
\[
d_A(,) \;=\; \| () - () \|_2.
\]
The pipeline reports this as  ( over the per-sample common subcomposition).

## AFDs: Gamma--Poisson (Negative Binomial) and alternatives

We model per-sample counts \(X_{is}\) for taxon \(i\) in sample \(s\) with depth \(N_s\). Let \(_i\) be the mean proportion across samples (estimated from unbiased factorial moments). The mean count is
\[
_{is} \;=\; _i\,N_s.
\]

### Negative Binomial parameterization

We use the mean--size parameterization:
\[
X  (,), 
[X]=,
(X)=+{}.
\]
The pmf:
\[
(X=n ,) \;=\;
{()\,n!}\,
({+})^{}
({+})^n.
\]
Zero probability and its Poisson limit (\(\)):
\[
p_0(,) = ({+})^{}
\;\; e^{-}.
\]

### Detection probability, occupancy-as-detection

Given \(_{is}=_i N_s\):
\[
(X_{is}>0) \;=\; 1 - p_0(_{is},_i).
\]
Empirical ``occupancy'' across samples is the mean detection probability:
\[
_i \;=\; {S}_{s=1}^S (1-p_0(_i N_s,_i)).
\]
(Here this is a  notion, not a hierarchical occupancy model.)

### Depth for target detection {N*(q)}}
Solve \((X>0)  q\) for \(N\). With \(=N\):
\[
N^{*}(q) \;=\;

-{}, & = ,\\[1em]
{}, & <.

\]

### Richness and slope vs depth

Let \(p_{0,i}(N)=p_0(_i N, _i)\). Then
\[
R(N) \;=\; _i (1-p_{0,i}(N)),

R'(N) \;=\; _i { N}(1-p_{0,i}(N)).
\]
Closed forms:
\[
R'(N) =

_i _i\,e^{-_i N}, & _i=,\\[0.5em]
_i _i(1+_i N}{_i})^{-(_i+1)}, & _i<.

\]

## Unbiased factorial-moment estimators

Given counts \(n_{is}\) with depths \(N_s\), define totals
\[
D_1=_s N_s,
D_2=_s N_s(N_s-1),
D_3=_s N_s(N_s-1)(N_s-2).
\]
Then for each taxon \(i\),
\[
_1=}{D_1},
_2=(n_{is}-1)}{D_2},
_3=(n_{is}-1)(n_{is}-2)}{D_3}.
\]
Finite-depth ``naive'' moments (on proportions \(n_{is}/N_s\)) have curvature away from the factorial moments; the code predicts that curvature and reports residuals with bootstrap CIs.

## Goodness-of-fit and diagnostics

 Bin counts and compare observed vs NB-expected frequencies (merged bins ensure expected \( 5\)); adjust df for fitted parameters; BH-FDR \(q\)-values are reported.

 Simulate zero fractions under fitted parameters; two-sided tail probability for observed zero fraction.

 Regress \(\) vs \(\) across taxa; test \(H_0:\) slope \(=2\).

 Anderson–Darling test for normality of \( _1\) (Shapiro--Wilk when \(n 5000\)).

## Agreement metrics between tools

For each sample, transform both abundance vectors with chosen transform, then report:
[leftmargin=1.1em]
- Pearson correlation (linear) and Spearman correlation (rank).
- : Euclidean norm of CLR differences on that sample’s common subcomposition.

Tip: with many zeros, PLR or ALR with a top-\(k\) pivot is more stable than CLR.

## Memory \& numerics (practical notes)

[leftmargin=1.1em]
- Zero probabilities use IEEE-safe exponent bounds (around \( 745\) for double precision) to avoid overflow/underflow.
- Presence matrices are avoided for large data; per-fold occupancies are computed by grouping where possible.
- Block bootstraps () respect clustered designs (subjects/sites/time).

## Workflow recipes (do-this-to-answer-that)

Below,  file names and flags to your context. All commands assume Python entrypoint .

### Am I sequencing deep enough? What depth hits 95\% detection for key taxa?

python p4_eco.py \
  --merged merged.tsv --counts counts.tsv \
  --tool-pair metaphlan,bracken \
  --afd nb --target-detect 0.95 \
  --out-dir results/depth

Read :
[leftmargin=1.1em]
- : the \(N^{*}(q)\) per taxon.
- : at median observed \(N\), do you meet the target?

  	extbf{If many taxa are insufficient}, consider higher depth or focus on taxa with higher \(_i\) or lower overdispersion (larger \(_i\)).

### How many more species do I gain if I increase depth by ?}

python p4_eco.py \
  --merged merged.tsv --counts counts.tsv \
  --afd nb --depths 5e6,1e7,2e7,4e7 \
  --out-dir results/richness

Inspect :
[leftmargin=1.1em]
- : expected richness at each depth.
- : marginal gain (\(dR/dN\)); when this is near zero, you’re in diminishing returns.

### Does NB fit my counts? Should I prefer Poisson--lognormal for some taxa?

python p4_eco.py \
  --merged merged.tsv --counts counts.tsv \
  --afd both --gof-alpha 0.05 --ppc-draws 1000 \
  --out-dir results/modelcheck

**Look at** :
[leftmargin=1.1em]
-  (BH-FDR): small values indicate lack of fit.
-  (AIC): NB vs PLN preference per taxon.
- : PPC on zero fraction (extreme tails suggest misfit).

If many failures: enable , check denominators, or accept PLN for those taxa.

### Which correlation (Pearson vs Spearman) should I use at a planned depth?

python p4_eco.py \
  --merged merged.tsv --counts counts.tsv \
  --forecast-depths 5e6,1e7,2e7 \
  --out-dir results/forecast

See  and the summary:
when predicted zero-fraction is low, Pearson is suggested; else Spearman.

### I have lots of zeros. How do I compare tools robustly?

python p4_eco.py \
  --merged merged.tsv --tool-pair metaphlan,bracken \
  --agreement-mode alr --alr-topk 3 \
  --out-dir results/agreement_alr

Or use PLR:

python p4_eco.py \
  --merged merged.tsv --tool-pair metaphlan,bracken \
  --agreement-mode plr --plr-maxpairs 50000 \
  --out-dir results/agreement_plr

### Decontam with negatives (prevalence-based is safer first)

python p4_eco.py \
  --merged merged.tsv --counts counts.tsv \
  --metadata meta.xlsx --sample-col SAMPLE_CODE \
  --decontam prevalence --prevalence-thresh 0.5 \
  --out-dir results/decontam_prev

If you must subtract:

python p4_eco.py \
  --merged merged.tsv --counts counts.tsv \
  --metadata meta.xlsx --sample-col SAMPLE_CODE \
  --decontam subtract --out-dir results/decontam_sub

Check  and .

### Clustered study (subjects/time/sites): honest uncertainty

python p4_eco.py \
  --merged merged.tsv --counts counts.tsv \
  --metadata meta.xlsx --sample-col SAMPLE_CODE \
  --bootstrap-block-col subject_id \
  --cv-folds 5 --cv-metric occupancy_mse \
  --out-dir results/clustered

This widens CIs (good!) and makes occupancy calibration out-of-sample.

## Flag guide (what they do \& how tweaking changes behavior)

### Compositional transforms

[leftmargin=1.1em]
- }:
CLR is powerful with few zeros; ALR stabilizes via a pivot; PLR is robust for sparse data.
-  (CLR/ALR only): larger values damp zero influence; smaller values emphasize rare parts (risking noise).
- , : pick a stable biological reference or average of top-\(k\) abundant taxa.
- : cap pairs for speed; larger caps reduce variance of the correlations.

### AFD fitting \& model checking

[leftmargin=1.1em]
- }: choose NB, PLN, or compare by AIC.
- : per-taxon \(_i\); needs adequate data per taxon; improves fit heterogeneity.
- : robust grid search for shared \(\).
- , : k-fold calibration for occupancy/zero fraction; higher \(k\) = more compute, less optimistic.
- : more draws = stabler PPC, higher runtime.
- : Hermite nodes; higher = more accurate PLN likelihood, slower.

### Detection planning \& forecasting

[leftmargin=1.1em]
-  \(q\): fraction for \(N^{*}(q)\); larger \(q\) larger \(N^{*}\).
- : compute \(R(N),R'(N)\) at these \(N\).
- : predict zero fractions and suggest Pearson vs Spearman per depth.

### Selection, metadata, negatives

[leftmargin=1.1em]
- : quick subsetting by sample IDs.
- , , , , , : inclusion/exclusion and negatives tagging.
- : fallback negatives detection by ID patterns.
- }, : recommended to start with prevalence; subtraction can over-correct at low \(N\).
- : model-guided filtering threshold before re-agreement (higher threshold \(\) fewer taxa, often higher agreement).

### Uncertainty \& stats

[leftmargin=1.1em]
- : CIs for curvature residuals.
- , , : clustered resampling and early stopping.
- , : GOF threshold and conservative multiplicity option.
- , : per-batch Taylor/MAD reporting.

## Interpreting \& validating

 If  $$ , relationship is monotonic but nonlinear; zeros likely dominate.

 Tiny \(_i\) indicate heavy overdispersion; cross-check GOF and PPC. Use \(N^{*}(q)\) and \(R(N)\) to plan depth. 

[leftmargin=1.1em]
- **Mock communities**: compare predicted vs observed \(N^{*}(q)\); check if \(R(N)\) tracks known richness as depth increases.
- **Spike-ins**: verify \(_i\) and \(_i\) recovery; PPC should not systematically reject spikes.
- **Train/test splits**: use  and compare predictive zero fractions to held-out.

## Caveats \& extensions

[leftmargin=1.1em]
- ``Occupancy'' here is detection; full occupancy--detection hierarchies are future work.
- Compositional constraints persist; for across-sample comparability in agreement metrics, pre-filter to a fixed subcomposition (e.g., taxa present in \(\)X\% samples).
- Phylogeny, ecological nulls, and hierarchical batch models are not yet integrated.

*Reproducibility tips:** pin , version input TSVs and metadata, keep  with the exact command used.
