
<!-- README.md is generated from README.Rmd. Please edit that file -->

This folder contains all input and output related to the sensitivity
analysis on hyper-prior specification in the multilevel HMM.

## Simulation files

The files `simV4_3st_Core*.R` are used to generate 100 independent
intensive longitudinal datasets for each sensitivity simulation
scenario. From the baseline scenarios, we vary the number of persons,
time points, and pairwise DKL (see below for exact specifications),
while fixing the true number of states to 3 and the number of dependent
variables to 4. The sensitivity analysis is performed in three blocks:

#### 1. Block 1 - Emission means

Sensitivity check on the specification of the hyper-priors on the
emission means. Specifically, the hyper-prior value of the variable and
state specific means `Mu_0` are varied to:

- `Mu_0` = equal across all states within a variable (mean of the Depmix
  inferred state means across one variable)
- `Mu_0` = shrunken such that the state means are closer together, with
  a factor of 0.5. For example, if before it would be 3, 5, 9 it is now
  4, 5, 7.

In addition, we vary the value of K0 (only for the scenario in which we
shrink the state means to be closer together) to:

- `K0` = 1 (default setting)
- `K0` = 20% of the sample size, meaning `n` \* 0.2.

#### 2. Block 2 - Between-subject variability

The hyper-parameters reflecting the between-subject variability in the
emission means and the transitions are set to reflect:

- Very little between-person variability:
  - Between-subject variance in emission means = 0.01^2
  - Between-subject variance in logit intercepts = 0.01^2
- Very large between-person variability (i.e., predicted means and
  transitions are all over the place across individuals):
  - Between-subject variance in emission means = 1^2
  - Between-subject variance in logit intercepts = 3^2

Important to note here is that both settings are very far from the
reality in our simulated dataset, so this truly relates to misspecified
prior values.

#### 3. Block 3 - Transition probabilities

Sensitivity check on the specification of the hyper-priors on the
transition probabilities, which are on the multinomial logit scale.
Specifically, the hyper-prior value of the logit intercepts are set such
that all state transitions (including self-transitions) have an equal
probability. In addition, we vary the value of K0 to:

- `K0` = 1 (default setting)
- `K0` = 20% of the sample size, meaning `n` \* 0.2.

Varying the pairwise DKL and number of persons and time points is done
by reading line by line from the file `Core*_settings.csv`, where the
columns contain the number of persons `n`, the number of time points
`n_t`, the pairwise KLD, and the settings of the hyper-parameter values
described above (see **Result tables** below for a detailed description
of each column).

## Post-processing

The file `Post-processing sensitivity.R` contains the code to process
the raw output simulation result files from the cluster into
`extracted_results_Sensitivity.RDS` (step 1), and further convert these
into files containing the relevant performance metrics only, which are
stored in the subfolder `results_tables`.

## Result tables

Results are contained in the subfolder `result_tables`, which contains
the following .RDS files:

- `Label_switch_proxy_Sensitivity.RDS`: information on severe label
  switching occurrence.
- `emission_performance_Sensitivity.RDS`: information on parameter
  estimation performance for the emission distribution.
- `gamma_performance_Sensitivity.RDS`: information on parameter
  estimation performance for the transition probabilities.
- `state_decoding_Sensitivity.RDS`: information on state decoding
  performance.

Each of the matrices starts with the following columns:

- `n`: number of subjects in the used sensitivity scenario (possible
  values: 30, 120)
- `n_t`: number of observations per subject in the used sensitivity
  scenario (possible values: 50, 100, 200)
- `KL_div`: KL divergence used in the sensitivity scenario (possible
  values: 5, 7)
- `mu_sc`: specification of the hyper-prior means for the emission
  distributions (0 = default setting, i.e., data-driven starting values
  from a pooled HMM; 1 = means set equal across states within a
  variable; 2 = state means within a variable shrunken to be closer
  together by a factor of 0.5)
- `mu_K0_sc`: specification of K0 for the emission mean prior (0 =
  default value of 1; 1 = 20% of the number of subjects, i.e., `n` \*
  0.2). Only varied when `mu_sc` = 2.
- `trans_unif`: specification of the hyper-prior means for the
  transition probabilities on the multinomial logit scale (0 = default
  setting, informed by self-transition probabilities; 1 = set such that
  all transitions including self-transitions are equally probable)
- `alph_K0_sc`: specification of K0 for the transition probability prior
  (0 = default value of 1; 1 = 20% of the number of subjects, i.e., `n`
  \* 0.2). Only varied when `trans_unif` = 1.
- `var`: specification of the hyper-parameters for the between-subject
  variability in both emission means and transition probabilities (0 =
  default setting; 1 = very small between-subject variability, with
  variance set to 0.01^2 for both emission means and logit intercepts; 2
  = very large between-subject variability, with variance set to 1^2 for
  emission means and 3^2 for logit intercepts)

#### Emission performance

This result table contains information on bias, label switching problem,
precision (empirical SE), and coverage of the 95% credibility interval
for both the state dependent emission distribution means and standard
deviations. For the means, bias is the form of the median absolute
relative bias. Relative because the values of the true state means vary
over the simulation runs, and cannot be compared between runs. Absolute
as the true state means are both positive and negative, as such we want
to inspect the magnitude of the bias and not direction. Median because
the relative bias is sensitive to small true values (resulting in very
large relative bias), using the median instead of mean relative bias,
this influence is mitigated. The label switching problem denotes the
proportion of instances for which the RMSE over the state means of one
dependent variable \< 0.20 (cut off was empirically investigated) for
ALL dependent variables in a particular model. Note that for the state
dependent means, the precision (empirical SE) cannot be obtained as the
true values vary over the simulation runs. For the SD, bias is included
as median relative bias to facilitate same plotting as the emission
means and usage of the 10% max relative bias criteria, and as mean bias
(true value is always 1). The two values differ slightly as the first is
based on median bias, and the second on mean bias.

#### Gamma performance

This result table contains information on bias, precision (empirical
SE), and coverage of the 95% credibility interval. As the true values of
the diagonal of the transition probability matrix are quite high (0.7
and 0.8) compared to the off-diagonal values (0.1 and 0.2), all results
are presented separately for the diagonal and off diagonal values of the
transition probability matrix. Performance metrics are not presented
separately for each separate transition probability parameter, as this
does not add information. Results are presented on the probability scale
to aid interpretation, note that estimation is on the multinomial logit
(MNL) scale. In addition, inspecting the MNL scale may present
misleading results as different combinations of MNL intercepts may
result in very similar probability values. As the true values of the
transition probabilities over the simulation runs and scenarios are
identical within one level of the factor ‘true number of states’, bias
is presented on the original scale of the parameter, as mean bias.

#### State decoding

This results table contains information on how well the Viterbi
algorithm (using the estimated model parameters) infers the true state
sequence. State decoding performance is presented using three metrics:

1)  `mean_prop_correct`: the proportion of times the correct state is
    assigned to a time point over all subjects and simulation runs
    within a give scenario,
2)  `mean_kappa`: as the proportion correct does not adjust for the fact
    that some of the correct state assignment could be due to change, I
    also included Cohen’s kappa, which corrects for chance,
    incorporating the number of possible outcome values (here, true
    number of states), and
3)  `mean_prop_over_20`: the proportion of times the correct state
    receives a state probability over 0.20, often used as extra
    performance criterion on state decoding to investigate when state
    assignment was not correct, at least the correct state received some
    support in being probable.
