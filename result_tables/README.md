
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Included files

Result tables are separate for the scenarios pertaining to a true number
of states 1, 2, and 3. Results are contained in the top .RDS file
`Result_table_*`, where \* takes on `1st`, `2st`, or `3st`.
`Result_table_*` is a list, containing the following matrix elements (=
result tables):

- `model_selection_*` : information on model fit criteria over each of
  the fitted states (1 to 4), both AICc and AIC.
- `emission_performance_*`: information on parameter estimation
  performance for the emission distribution.
- `gamma_performance_*`: information on parameter estimation performance
  for the emission distribution.
- `state_decoding_*`: information on state decoding performance.

Each of the matrices starts with the following columns:

- true_m: true number of states in the used scenario,
- n: number of subjects in the used scenario (possible values: 15, 30,
  60, 120)
- n_t: number of obervations per subject in the used scenario (possible
  values: 50, 100, 200, 400, 800)
- n_dep: number of dependent variables used in the scenario (possible
  values: 2, 4, 8)
- KL_div: KL divergence in the used in the scenario (possible values: 3,
  5, 7)

Note: for the scenarios pertaining to a true number of states = 1, the
elements `gamma_performance_1st` and `state_decoding_1st` are omitted as
they are obsolete. In additon, for the emission distribution
performance, information on the label switching problem is omitted.

### Model selection

This result table contains two types of info: means and proportions.
Means pertain to the average AICc/AIC over the 100 simulation runs for
that particular number of states model. The proportions denote the
proportion that this particular state solution resulted in the lowest
AICc/AIC over the 100 simulations.

### Emission performance

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

### Gamma performance

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

### State decoding

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
