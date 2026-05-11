#import "utils/general-utils.typ": * 
#import "template/styles.typ": *

#set document(title: "Thesis Proposal") // Note: this title is visible in the PDF viewer

#show: styles

#set align(center)
#text(
  heading("Combining overlap correction and back-to-back regression methods to analyse eye-tracking-EEG co-registration data during natural reading", numbering: none, outlined: false), size: 1.15em
)
#v(14pt)

#text("<MSc> Thesis Proposal - <Yanhong Xu>
May 2026")
\


#line(length: 100%, stroke: gray)

#set align(left)
#set heading(numbering: "1.")


= Introduction

== Motivation
Natural reading is a complex cognitive process involving many causal and non-causal factors that may be partially correlated and unfold across at visual, attentional, lexical-semantic, and oculomotor levels. Co-registration of eye movements (EMs) and electroencephalography (EEG) enables the investigation of fixation-related potentials (FRPs) during natural reading. A key challenge in FRP data is that rapid fixations give rise to strong temporal overlap, while lexical factors are often highly collinear. Disentangling the respective contributions of correlated factors, such as word surprisal, word frequency, and word length, from multivariate EEG signals therefore remains a central methodological challenge.



#cite(<king_back--back_2020>, form: "prose") introduced back-to-back regression (B2B), a linear regression framework designed to disentangle covarying predictors in multidimensional neural data and to obtain interpretable predictor-wise estimates. Yet in EMs–EEG co-registration, rapid event sequences elicit temporally overlapping responses, which can bias coefficient estimates and further complicate inference. A principled pipeline that jointly addresses temporal overlap and predictor collinearity is therefore needed, and it remains unclear how robust and interpretable predictor estimates can be obtained in practice.

This project evaluates whether overlap correction using linear deconvolution within a regression event-related potential framework (rERP), implemented in Unfold @ehinger_unfold_2019, improves B2B’s ability to disentangle collinear predictors in fixation-locked EEG during natural reading, using simulations with known ground truth and a real-data application to _Zurich Cognitive Language Processing Corpus_ (ZuCo) @hollenstein_zuco_2018. To the best of our knowledge, this combination has not been systematically benchmarked in the literature. This work also aims to support more interpretable analyses of NLP-derived lexical predictors (e.g., surprisal) in FRPs, potentially informing both linguistic and cognitive perspectives on natural reading. 


== Related work
The ZuCo dataset provides co-registered eye-tracking and EEG data during natural reading, making it a suitable benchmark for fixation-locked analyses with lexical and oculomotor predictors. In this project, ZuCo is used as a real-data application. The relevant methodological background on overlap correction and B2B is summarised below.



=== FRPs and temporal overlap in natural reading 
FRPs are event-related EEG responses time-locked to the onset of eye fixations, with FRP data, it is possible to investigate neural correlates of specific aspects of oculomotor behaviour, while also providing insight into the processes underlying reading. For instance, as shown in @fig:frp, reading-related EEG responses include well-established components such as the N400, which is commonly associated with late semantic integration and is often modulated by lexical variables such as word frequency and predictability, as well as earlier visual components linked to occipital processing @dimigen_coregistration_2011 @degno_eye_2020.

However, FRP estimates in eye-tracking–EEG co-registration are complicated by two related challenges: temporally overlapping responses from successive fixations, and multiple lexical and oculomotor covariates that can modulate the EEG signal and covary with one another. These challenges motivate regression-based approaches to FRP analysis.

#figure(
  grid(
    columns: (1fr, auto),
    gutter: 6pt,
    [
      #image("images/FRP.png", width: 100%)
    ],
    [     
      #align[#v(70pt)
          #rotate(-90deg)[
            #text(size: 7pt)[Source from #cite(<dimigen_coregistration_2011>, form: "prose")]
          ]
      ]
    ],
  ),
  caption: [
    *grand average FRP in natural reading.* _λ: early visual component. N400: lexico-semantic response component_.
  ],
) <fig:frp>

=== Regression-based ERP and linear deconvolution
Regression-based ERP methods, or rERPs, estimate predictor-specific EEG effects by modelling event-related responses within a regression framework rather than relying only on condition-wise averaging, which allows predictors of interest and control covariates to be included in the same model. 

However, in eye-tracking-EEG co-registration, not only multiple factors and covariates, but also temporally overlapping responses can modulate the EEG signal and confound event-related estimates. Thus, an additional complication in co-registration is how to correct for temporal overlap in fixation-related EEG responses. 

#cite(<ehinger_unfold_2019>, form:"prose") developed the Unfold toolbox, which implements deconvolution-based rERP that extends the rERP framework by modelling continuous EEG with a time-expanded design matrix, allowing overlapping event-related responses to be estimated jointly and provides overlap-corrected estimates, thereby enabling less biased and more reliable interprtation of event-related effects @dimigen_regression-based_2021.







=== General Linear Model (GLM)-reconstructed single trials for decoding

Standard deconvolution-based rERP analyses typically return time-resolved beta coefficients, that is, estimated response waveforms for specific predictors or event types. These estimates are suitable for interpreting overlap-corrected FRP effects, but they are not event-level EEG epochs that can be directly used in a standard epoched decoding analysis.

#cite(<r_decoding_2023>, form: "prose") combined GLM-based deconvolution with multivariate pattern analysis (MVPA) by reconstructing overlap-corrected single trials. The key idea is to use the GLM estimates to reconstruct the activity overlapping with a target fixation, subtract this activity from the continuous EEG signal, and then epoch the target fixation. The resulting GLM-reconstructed single trial preserves a single-fixation EEG observation while reducing contamination from temporally overlapping responses.

This distinction is important for the present project because the two overlap-corrected B2B pipelines use deconvolution differently. In the one-step pipeline, overlap correction is incorporated directly into the B2B estimation procedure. In the two-step pipeline, GLM-based deconvolution is first used to reconstruct overlap-corrected single-fixation EEG responses, and B2B is then applied to these reconstructed event-level responses.


=== Back-to-Back regression <sec:b2b>
B2B consists of a decoding step and an encoding step (@fig:b2b), typically using regularised linear solvers (e.g., ridge regression or a linear SVM). To avoid overfitting, the samples $(X, Y)$ are repeatedly split into disjoint subsets $(X_1, Y_1)$ and $(X_2, Y_2)$.

In the first step, a backward linear regression model predicts coefficients $hat(G)$ from multivariate neural observations $Y_1$ to predictors $X_1$. The matrix $hat(G)$ defines feature-specific supervised readout dimensions in $Y$-space. The decoded estimates are then used to construct the predictions $hat(X) = Y_2 hat(G)$. In the second step, a forward linear regression model uses the original predictors $X_2$ to predict $hat(X)$, yielding $hat(H)$, so that $X_2hat(H)$ best approximates $hat(X)$. The diagonal of $hat(H)$ is the factor-wise, interpretable estimates:

$
hat(G) &= (Y_1^T Y_1 + Lambda_Y)^(-1) Y_1^T X_1 \
hat(H) &= (X_2^T X_2 + Lambda_X)^(-1) X_2^T Y_2 hat(G) \
hat(S) &= "diag"(hat(H))
$


B2B combines both the benefits of forward and backward modelling to provide interpretable and unidimensional coefficient estimates for each factor. B2B has been shown to be a promising analytical method, and it was applied to a reading task to decompose different word features that relate to different brain responses. In their results, word length and word frequency were identified as key effects, with evidence consistent with earlier effects in visual processing and later evoked N400 responses.


#figure(
  grid(
    columns: (1fr, auto),
    gutter: 6pt,
    [
      #image("images/b2b.png", width: 100%)
    ],
    [     
      #align[#v(70pt)
          #rotate(-90deg)[
            #text(size: 6pt)[CC BY-NC-ND 4.0 · King et al. (2020)]
          ]
      ]
    ],
  ),
  caption: [
    *B2B algorithm overview.* _F denotes an encoding mapping capturing the spatiotemporal EEG pattern associated with each predictor_.
  ],
) <fig:b2b>





\
= Planned Project
== Research Question
1. In simulations mimicking co-registration data with known ground truth, does combining overlap correction (deconvolution rERP) with B2B improve the disentanglement and recovery of collinear predictors when fixation-related responses overlap in time?

2. When applied to ZuCo, does the best-performing pipeline yield stable and interpretable estimates for lexical predictors (surprisal, frequency, length) while controlling for oculomotor covariates?



== Goals



=== Main Goals <mainGoals>
#v(0.3em)
#set enum(numbering: "A.")

#[
  #show figure: set align(left) 
  + #goal("A simulation framework with controllable temporal overlap and predictor collinearity is implemented using Unfold toolbox.") <goal1>  
  + #goal([Analysis pipelines spanning the overlap × collinearity conceptual method map (@fig:methodmap), including baseline decoding, deconvolution rERP decoding, B2B, and rERP + B2B variants, are systematically benchmarked against known ground truth.]) <goal2>
  + #goal("A fixation-aware surprisal predictor is implemented, distinguishing first fixations and refixations.") <goal3>
  + #goal("The selected best-performing pipeline(s) will be applied to the ZuCo dataset, including lexical predictors such as word surprisal (from a pretrained language model), log word frequency, and word length, while controlling for oculomotor covariates.") <goal4> 
]


=== Stretch Goals <stretchGoals>
#v(0.3em)
#set enum(numbering: "A.", start: 5) // continue the numbering from where the main goals left off. Adjust `start` depending on how many main goals you have.
#[
  #show figure: set align(left) 
  + #goal("Different solver variants in the B2B backward step (L2/ridge vs. linear SVM) are evaluated and compared to assess how solver choice affects decoding accuracy and result stability in the simulations.") <goal5>
]

== Approach <approach>

=== Overview

#figure(
  image("images/b2b_project_workflow.svg"),
  caption: [*Project workflow overview*: _Simulation benchmarking → select best pipelines → ZuCo application_]
)

The project is organised in two stages. First, simulations with known ground truth will be used to evaluate whether overlap correction improves B2B’s performance under controlled temporal overlap and predictor collinearity (*@goal1*,*@goal2*,*@goal5*). Second, the best-performing pipeline(s) will be applied to a real natural reading dataset (ZuCo) as a practical demonstration. (*@goal3*, *@goal4*)

=== Toy simulation sanity check

Before implementing the full simulation benchmark, a minimal toy simulation is used as a sanity check for the B2B logic.

In this toy simulation, the predictor matrix $X$ contains three features, $x_1$, $x_2$, and $x_3$. The observations $Y$ are generated only from the truly contributing predictors $x_1$ and $x_2$, whereas $x_3$ is correlated with $x_1$ but has no direct contribution. Specifically, $x_3 = 0.8x_1 + 0.4 epsilon$, so $x_3$ contains both a component shared with $x_1$ and a residual component that is specific to $x_3$.

This setup illustrates the distinction between decodability and contribution. Because $x_3$ is correlated with $x_1$, it can still be partly decoded from $Y$. However, the $x_3$-specific residual never contributes to $Y$, and should therefore not be treated as a true recovered contribution. 

The decoded representation $hat(X) = Y G$ is therefore constrained by the information present in the observations $Y$. Since the recoverable signal is generated from $x_1$ and $x_2$, $hat(X)$ should mainly reflect a low-dimensional recoverable structure rather than the full original predictor space. As shown in @fig:toy #footnote[
  In this toy figure, $X_("test")$ and $Y_("test")$ correspond to the B2B subsets $X_2$ and $Y_2$ from @sec:b2b..
], $X_("test")H$ and $Y_("test")G$ are approximately aligned on a plane with this recoverable structure. Thus, B2B aims to suppress the apparent contribution of the correlated but non-causal feature $x_3$, rather than preserving all variance in $X$.




#figure(
  image("images/b2b_continous_3effects_5panels.svg"),
  caption: [*Minimal toy simulation:* illustrating that a correlated but non-causal feature can be decodable, while B2B aims to recover only the features with true contribution.]
)<fig:toy> 



=== Simulation design (*@goal1*)
Synthetic co-registration-like data are generated to mimic fixation-locked responses with controllable overlaps and correlations. Two main factors are manipulated:

- Temporal overlap: none vs. present (implemented by varying event rate / inter-event interval (ISI) distributions).
- Predictor collinearity: none / low / high (implemented by controlling correlations among predictors, using covariate design in Unfold toolbox.)


=== Analysis pipelines (*@goal2*)

#figure(
  image("images/method_map.svg", width: 80%),
  caption: [*Conceptual method map: overlap x collinearity.*]
) <fig:methodmap>

The method map provides a conceptual overview of candidate analysis strategies along the dimensions of overlap correction and predictor disentanglement. In total, five representative pipelines are compared:

- *Standard decoding:* a baseline decoding model without explicit overlap correction or B2B.
- *Plain B2B:* B2B applied directly to fixation-locked data without explicit overlap correction.
- *Overlap-corrected decoding:* deconvolution-based rERP decoding without B2B.
- *One-step overlap-corrected + B2B:* a joint approach in which linear deconvolution-based overlap correction is incorporated directly into the B2B estimation procedure.
- *Two-step overlap-corrected → B2B:* a sequential pipeline in which a deconvolution GLM is first used to reconstruct overlap-corrected single-fixation EEG trials. B2B is then applied to these reconstructed event-level responses.

=== Evaluation criteria (*@goal2*, *@goal5*)
- *Simulation evaluation:* In simulations, ground truth is defined in terms of multi-source predictor-specific waveform shapes. Since the simulated predictors may be correlated, evaluation must assess not only whether true effects are recovered, but also whether effects are assigned to the correct predictors. The primary simulation metric will therefore be a cross-talk matrix. For each recovered predictor-specific waveform and each true waveform, the matrix entry is defined as: 

  $ M_(i j) = upright("corr")(hat(w)_i(t), w_j^(upright("true"))(t)) $

  where $hat(w)_i(t)$ denotes the waveform recovered for predictor $i$ and $w_j^(upright("true"))(t)$ denotes the true simulated waveform for predictor $j$. Ideally, the matrix should show high values on the diagonal and low values off the diagonal, indicating that each true effect is recovered for the correct predictor with little leakage to other predictors.

- *Real-data evaluation:* For the ZuCo application, ground-truth neural effects are unknown. Evaluation will therefore rely on held-out predictive performance and predictor knockout analyses. This knockout procedure follows the feature-importance evaluation proposed in the original B2B framework, illustrated in @fig:evaluation. For each lexical predictor, the contribution of that predictor will be estimated by comparing the performance of the full model with a model in which the predictor is knocked out:
  $
  Delta R_i = R_("full") - R_("knockout"(i))
  $


  Here, $R$ denotes the correlation-based held-out prediction score. A larger positive $Delta R_i$ indicates that predictor $i$ contributes to held-out prediction performance. Statistical testing will be used to assess the reliability of the knockout-based performance change.

#figure(
  grid(
    columns: (1fr, auto),
    gutter: 6pt,
    [
      #image("images/corr.png", width: 80%)
    ],
    [     
      #align[#v(70pt)
          #rotate(-90deg)[
            #text(size: 6pt)[CC BY-NC-ND 4.0 · King et al. (2020)]
          ]
      ]
    ],
  ),
  caption: [*Held-out performance algorithm from B2B paper*
  ],
) <fig:evaluation>

=== Word surprisal (*@goal3*)
Word surprisal will be computed using a pretrained transformer language model, such as GPT-2, via Hugging Face. For each word, surprisal is defined as the negative log-probability of each word given its preceding context, $-log p(w_i | w_(<i))$. The resulting word-level values will then be aligned to the corresponding fixation events in ZuCo and included in the model as a lexical predictor. A fixation-aware variant will distinguish first fixations from refixations, so that the same word-level predictor can be modelled differently depending on fixation status.


=== Real-data application: ZuCo (*@goal4*)
The best-performing configuration(s) from the simulation benchmark will be applied to the ZuCo dataset. Fixation-locked EEG responses will be modelled using the selected pipeline(s), with lexical predictors including word surprisal, log word frequency, and word length, while controlling for oculomotor variables such as saccade amplitude, refixation status, and landing position. The real-data application will assess whether the selected method can disentangle the influence of correlated lexical factors from multivariate fixation-locked EEG observations during natural reading.



=== Optional extension: solver comparison (*@goal5*) 
Solver variants in the B2B backward step may be compared to assess sensitivity to modelling choices, focusing on L2/ridge regression versus linear SVM.

\

= Plan
\
- Review the relevant literature on B2B, FRPs, and the deconvolution rERP.
- Implement a simulation benchmark with known ground truth, and evaluate alternative pipelines under controlled overlap and collinearity. 
- Apply the best-performing configuration(s) to the ZuCo dataset. 
- Result consolidation, documentation, and thesis writing. 



#timeliney.timeline(
  show-grid: true,
  {
    import timeliney: *
      
    headerline(group([*Apr*],[*May*],[*June*],[*Jul*],[*Aug*],[*Sep*]))
    
    task("Literature review, writing proposal", (0, 1), style: (stroke: 2pt + gray))
    task("Simulation + analysis (Goal A-B)", (0.5, 3), style: (stroke: 2pt + gray))
    task("ZuCo application (Goal C-D)", (2.5, 4.5), style: (stroke: 2pt + gray))
    task("Buffer / Review / writing thesis / (Goal E)", (4,6), style: (stroke: 2pt + gray))

    milestone(
      at: 4,
      style: (stroke: (dash: "dashed")),
      align(center, [
        *Main goal completion*\
        Sep 2026
      ])
    )
  }
)




#line(length: 100%, stroke: gray)

#bibliography("refs.bib", style: "american-psychological-association")
