// A central place where libraries are imported (or macros are defined)
// which are used within all the chapters:
#import "utils/global.typ": *


// Fill me with the Abstract
#let abstract = [#lorem(150)]

// Fill me with acknowledgments
#let acknowledgements = [#lorem(50)]


// if you have appendices, add them here
#let appendix = [
  = Appendices
  //#include "./chapters/appendix.typ"
]

// Put your abbreviations/acronyms here.
// 'key' is what you will reference in the typst code
// 'short' is the abbreviation (what will be shown in the pdf on all references except the first)
// 'long' is the full acronym expansion (what will be shown in the first reference of the document)
//
// In the text, call @eeg or @uniS to reference  the shortcode
#let abbreviations = (
  (
    key: "eeg",
    short: "EEG",
    long: "Electroencephalography",
  ),
  (
    key: "uniS",
    short: "UoS",
    long: "University of Stuttgart",
  ),
)

#show: thesis.with(
  author: "<author>",
  title: "<title>",
  degree: "<degree>",
  faculty: "Faculty of Electrical Engineering and Computer Science",
  department: "Computational Cognitive Science",
  major: "<major>",
  supervisors: (
    (
      title: "Main Supervisor",
      name: "Benedikt Ehinger",
      affiliation: [Computational Cognitive Science \
        Faculty of Electrical Engineering and Computer Science, \
        Department of Computer Science
      ],
    ),
    (
      title: "Second Supervisor",
      name: "Max Mustermann",
      affiliation: [Computational Cognitive Science \
        Faculty of Electrical Engineering and Computer Science, \
        Department of Computer Science
      ],
    ),
  ),
  epigraph: none,
  abstract: abstract,
  appendix: appendix,
  acknowledgements: acknowledgements,
  preface: none,
  figure-index: false,
  table-index: false,
  listing-index: false,
  abbreviations: abbreviations,
  date: datetime(year: 2025, month: 6, day: 1),
  bibliography: bibliography("refs.bib", title: "Bibliography", style: "american-psychological-association"),
)

// Code blocks
#codly(
  languages: (
    rust: (
      name: "Rust",
      color: rgb("#CE412B"),
    ),
    // NOTE: Hacky, but 'fs' doesn't syntax highlight
    fsi: (
      name: "F#",
      color: rgb("#6a0dad"),
    ),
  ),
)

// If you wish to use lining figures rather than old-style figures, uncomment this line.
// #set text(number-type: "lining")

// import custom utilities
#import "utils/general-utils.typ": *

// Main Content starts here
= Introduction <chp:introduction>

Why does the thesis matter? Give context for the thesis, explain the research gap and provide a summary.


#underline[*Note:*] The suggested structure in this template is loosely based on the summarized figure presented in #pc[@Mensh2017-vd]. It can be adapted to suit the requirements of a Bachelor/Master thesis. 

Below is the updated version of the figure (citation included in the caption).

#figure(
  image("template/assets/pcbi.1005830.g001.png", width: 70%),
  caption: "Figure describing a suggested paper structure.
    Reproduced from " + pc[@The_PLOS_Computational_Biology_Staff2017-gr] + "."
)


== Background
Describe the background of the thesis topic, along with the research gap to be filled or the problem to be solved.

== Research questions
You may write about your research questions here.

== Summary
Sum up the approach and results.

//
// NOTE:
// It's important to have explicit pagebreaks between each chapter,
// otherwise header stylings from the template might break
#pagebreak()
= Chapter on a specific topic
Write about a particular topic if you think it needs a separate chapter, e.g. an important concept having a lot of background to cover.

#pagebreak()
= Methods
Describe the method(s) used.

== Execution details
If applicable, explain how the method was followed.

== Quick tips on using this template

By default, using the `@` syntax to reference a source from the bibliography shows the citation in brackets @Mensh2017-vd. To get a prose citation, use the utility function `#pc[@source]` like so: #pc[@Mensh2017-vd].

Mark TODOs using the `#todo("Do this!")` function. It shows up in large red letters #todo("Do this!") and breaks up the paragraph.

Inline TODOs can be added using the custom utility function`#ilt[TODO text]` which shows the #ilt[TODO text] within the same line yet clearly standing out.

#pagebreak()
= Results
Present, describe and summarize the results obtained. 


#pagebreak()
= Discussion

== Interpretation
Write about the interpretation of the results and explain how the identified research gap was filled.

== Limitations
What limitations does the interpretation have? How can this be fixed?

== Outlook / Future Scope
What difference was made by the work done on the thesis, and how can it be expanded on in the future?


#pagebreak()
= Summary
Recap what was discussed.

= Other chapters
The further chapters give detailed information on the template. They have been included from a separate typst file and can be omitted by commenting out or deleting the `#include` line following this chapter in the source code. 

#include "template/demo/main.typ"
