#import "utils/general-utils.typ": * 
#import "template/styles.typ": *

#set document(title: "Thesis Proposal") // Note: this title is visible in the PDF viewer

#show: styles

#set align(center)
#text(
  heading("Combine overlap correction and back-to-back regression methods to analyse eye-tracking-EEG co-registration data during natural reading", numbering: none, outlined: false), size: 1.15em
)
#v(14pt)

#text("<MSc> Thesis Proposal - <Yanhong Xu>
March 2026")
\


#line(length: 100%, stroke: gray)

#set align(left)
#set heading(numbering: "1.")


= Introduction

== Motivation
Natural reading as a complicated cognitive process, contains numerous causal and non-causal factors, which can be partially correlated, unfolding at visual, attentional, lexical-semantic, and oculomotor levels. Eye movements (EMs) and electroencephalogram (EEG) co-registration research provides access to  fixation-related potentials (FRPs) during natural reading. A key challenge in FRP data is that rapid fixations induce strong temporal overlap, while lexical factors are often highly collinear. At the same time, a central methodological challenge is to disentangle the respective contributions of collinear factors (e.g., word surprisal, word frequency, word length) from multivariate EEG signals.




#cite(<king_back--back_2020>, form: "prose") introduced  a back-to-back regression (B2B) to disentangle covarying predictors in multidimensional neural data and to obtain interpretable factor-wise estimates. Yet in EMs–EEG co-registration, rapid event sequences elicit temporally overlapping responses, which can bias coefficient estimates and further complicate inference. A principled pipeline that jointly addresses temporal overlap and predictor collinearity is therefore needed, and it remains unclear how robust and interpretable factor estimates can be obtained in practice.

This project benchmarks whether overlap correction via linear deconvolution-based regression event-related potential (rERP) (implemented in Unfold) @ehinger_unfold_2019 improves B2B’s ability to disentangle collinear factors in fixation-locked EEG during natural reading, using simulations with known ground truth and a real-data application to ZuCo @hollenstein_zuco_2018. To the best of current knowledge, this combination has not been systematically benchmarked in the literature. This work also aims to support more interpretable analyses of NLP-derived lexical predictors (e.g., surprisal) in FRPs, potentially informing both linguistic and cognitive perspectives on natural reading. 


== Other work
_Zurich Cognitive Language Processing Corpus_ (ZuCo) dataset provides co-registered eye-tracking and EEG data during natural reading, making it a suitable benchmark for fixation-locked analyses with lexical and oculomotor predictors. In this project, ZuCo is used as a real-data application. 
Together with overlap correction and B2B, the relevant background is summarised below.


=== Temporal overlap in FRPs (Overlap)
FRPs are event-related EEG responses time-locked to the onset of eye fixations, with FRP data, it is possible to investigate neural correlates of specific and different aspects of oculomotor behaviour, while also providing novel scientific insight into the processes underlying reading. For example, reading-related EEG responses include well-established components such as the N400, which is commonly associated with late semantic integration and is often modulated by lexical variables such as word frequency and predictability, as well as earlier visual components linked to occipital processing @dimigen_coregistration_2011 @degno_eye_2020.

An important issue in co-registration is how to separate temporally overlapping neural components in the EEG. Because multiple covariates can modulate the  EEG signal and confound event-related estimates, #cite(<ehinger_unfold_2019>, form:"prose") introduced the Unfold toolbox, which intergrates a regression-ERP(rERP) approach based on linear deconvolution(deconvolution rERP), enabling less biased estimation of event-related effects and more reliable interpretation.

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

=== Back-to-Back regression (B2B)
#cite(<king_back--back_2020>, form: "prose") developed a B2B linear regression method to disentangle the respective contributions of covarying features from multidimensional observations, which combines both the benefits of forward and backward modelling to provide interpretable and unidimensional coefficient estimates for each factor. B2B has been shown to be a promising analytical method, and it was applied to reading task to decompose different word features that relate to different brain responses. In their results, word length and word frequency emerged as key effects, with evidence consistent with earlier effects in visual processing and later evoked N400 responses. 

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

Decoding in B2B refers to a backward model that predicts putative predictors from multivariate neural observations (i.e., $Y -> hat(X)$), typically using a regularised linear solver (e.g., ridge regression or a linear SVM). The decoded estimates are then used in a forward step to obtain factor-wise, interpretable importance estimates.




\
= Planned Project
== Research Question
1. In simulations mimicking co-registration data with known ground truth, does combining overlap correction (deconvolution rERP) with B2B improve the disentanglement and recovery of collinear predictors under temporal overlap?
2. How does solver choice in the B2B backward step  (L2/ridge vs. linear SVM) affect decoding accuracy and result stability in these simulations?
3. When applied to ZuCo, does the best-performing pipeline yield stable and interpretable estimates for lexical predictors (surprisal, frequency, length) while controlling for oculomotor covariates?



== Goals



=== Main Goals <mainGoals>
#v(0.3em)
#set enum(numbering: "A.")

#[
  #show figure: set align(left) 
  + #goal("A simulation framework with controllable temporal overlap and predictor collinearity is implemented using Unfold toolbox.") <goal1>  
  + #goal([Analysis pipelines spanning the overlap × collinearity method map (@fig:methodmap), including baseline decoding, deconvolution rERP, B2B, and deconvolution rERP + B2B variants, are systematically benchmarked against known ground truth.]) <goal2>
  + #goal("A fixation-aware surprisal predictor is implemented, distinguishing first fixations and refixations, and is compared against word-linked surprisal in the ZuCo analysis.") <goal3>
  + #goal("The best-performing pipeline is applied to the ZuCo dataset, including lexical predictors such as word surprisal (from a pretrained language model), log word frequency, and word length, while controlling for oculomotor covariates.") <goal4> 
]


=== Stretch Goals <stretchGoals>
#v(0.3em)
#set enum(numbering: "A.", start: 5) // continue the numbering from where the main goals left off. Adjust `start` depending on how many main goals you have.
#[
  #show figure: set align(left) 
  + #goal("Solver variants in the B2B backward step (L2/ridge vs. linear SVM) are evaluated and compared.") <goal5>
]

== Approach <approach>

=== Overview

#figure(
  image("images/b2b_project_workflow.svg"),
  caption: [*Project workflow overview*: _Simulation benchmarking → select best pipelines → ZuCo application_]
)

The project is organised in two stages. First, simulations with known ground truth are used to evaluate whether overlap correction improves B2B’s performance under controlled temporal overlap and predictor collinearity (*@goal1*,*@goal2*,*@goal3*). Second, the best-performing pipeline is applied to a real natural reading dataset (ZuCo) as a practical demonstration. (*@goal4*, *@goal5*)

=== Simulation design (*@goal1*)
Synthetic co-registration-like data are generated to mimic fixation-locked responses with controllable overlaps and correlations. Two main factors are manipulated:

- Temporal overlap: none vs. present (implemented by varying event rate / inter-event interval (ISI) distributions).
- Predictor collinearity: none / low / high (implemented by controlling correlations among predictors, using covariate design in Unfold toolbox.)


=== Analysis pipelines (*@goal2*)

#figure(
  image("images/method_map.png"),
  caption: [*Method map: overlap x collinearity.*]
) <fig:methodmap>

Five pipelines are compared:

- *Standard decoding:* a baseline decoding model without explicit overlap correction or B2B.
- *Plain B2B:* B2B applied directly to fixation-locked data without explicit overlap correction.
- *Overlap correction only:* deconvolution-based rERP in Unfold without B2B.
- *Two-step overlap-corrected B2B:* overlap correction via deconvolution-based rERP, followed by B2B on overlap-corrected single-trial responses.
- *one-step overlap-corrected B2B:* a joint approach in which overlap correction is incorporated directly into the B2B estimation procedure.


=== Evaluation criteria (*@goal2*, *@goal5*)
- In simulations, ground truth is defined in terms of predictor-specific waveform shapes; performance is therefore evaluated by the time-course similarity between $hat(beta)_i (t)$ and $beta_i^("true")(t)$, using e.g. correlation or RMSE.

- In ZuCo, where ground truth is unknown, lexical predictor effects will be assessed using knockout-based held-out performance loss:
$
Delta R_i = R_("full") - R_("knockout"(i))
$

=== Word surprisal (*@goal3*)
Word surprisal ($-log p(w_i | w_(<i)))$, derived from a pretrained language model as the negative log-probability of each word given its preceding context, will be aligned to the corresponding fixation events in ZuCo and included in the model as lexical predictors. A fixation-aware variant will additionally distinguish first fixations from refixations. 


=== Real-data application: ZuCo (*@goal4*)
The best-performing configurations are applied to the ZuCo dataset. Fixation-locked events are modelled with Unfold for overlap correction, and B2B is then used to assess factor estimates for lexical predictors (word surprisal, log word frequency and word length) while controlling for oculomotor variables (saccade amplitude, refixation, landing position).



=== Solver comparison (*@goal5*) 
To assess sensitivity to modelling choices, solver variants in the B2B backward step are compared, focusing on L2/ridge versus linear SVM.

\

= Plan
\
- Review the relevant literature on B2B, FRPs, and deconvolution rERP.
- Implement a simulation benchmark with known ground truth, and evaluate alternative pipelines and solver variants under controlled overlap and collinearity. 
- Apply the best-performing configuration to the ZuCo dataset. 
- Result consolidation, documentation, and thesis writing. 



#timeliney.timeline(
  show-grid: true,
  {
    import timeliney: *
      
    headerline(group([*Mar*],[*Apr*],[*May*],[*June*],[*Jul*],[*Aug*]))
    
    task("Literature review, writing proposal", (0, 1), style: (stroke: 2pt + gray))
    task("Simulation + analysis (Goal A-B)", (0.5, 3), style: (stroke: 2pt + gray))
    task("Zuco application (Goal C-D)", (2.5, 4), style: (stroke: 2pt + gray))
    task("Buffer / Review / Stretch goals / (Goal E)", (4,6), style: (stroke: 2pt + gray))

    milestone(
      at: 4,
      style: (stroke: (dash: "dashed")),
      align(center, [
        *Main goal completion*\
        Aug 2026
      ])
    )
  }
)




#line(length: 100%, stroke: gray)

#bibliography("refs.bib", style: "american-psychological-association")
