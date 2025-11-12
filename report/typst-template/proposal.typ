#import "utils/general-utils.typ": * 
#import "template/styles.typ": *

#set document(title: "Thesis Proposal") // Note: this title is visible in the PDF viewer

#show: styles

#set align(center)
#text(
  heading("A Bachelor/Master Thesis Proposal", numbering: none, outlined: false), size: 1.15em
)
#v(14pt)

#text("<degree> Thesis Proposal - <student name>
June 2025")
\


#line(length: 100%, stroke: gray)

#set align(left)
#set heading(numbering: "1.")


= Introduction

== Motivation
This is the motivation for your project. Why is it interesting?

== Other work
Here you add all the background findings, citing different studies and the like. You do not need to include text-book knowledge e.g. how exactly EEG works, except if you are working on the topic and require the detail. Think of writing this proposal for a fellow student, but one that was in the EEG/ET-courses already.

=== Subheading 1
Use this if you need to write about a topic in more detail. You can also refer to sources in your bibliography @QayyumY-sac2022. Prose citations, e.g. "according to #pc[@QayyumY-sac2022]" are also possible.

=== Subheading 2
Yet another subtopic.

\
= Planned Project
== Research Question
What is the main guiding question(s) you want to answer? For example,
+ Is it easy to generate research questions?
+ How many research questions does a proposal need?
+ What is the answer to the question of Life, the Universe, and Everything?

== Goals
Here you define the goals you really want to reach, and the stretch goals. Please link the respective goal to the respective paragraphs in the section `Approach` (@approach), via e.g. wording or using the syntax shown.

=== Main Goals <mainGoals>
#v(0.3em)
#set enum(numbering: "A.")

#[
  #show figure: set align(left) 
  + #goal("I want to finish this Proposal.") <goal1>  
  + #goal("I want to finish my Thesis.") <goal2>
]


=== Stretch Goals <stretchGoals>
#v(0.3em)
#set enum(numbering: "A.", start: 3) // continue the numbering from where the main goals left off. Adjust `start` depending on how many main goals you have.
#[
  #show figure: set align(left) 
  + #goal("I'll do another thesis if there is time.") <goal3>
  + #goal("Maybe another one as well.") <goal4>
]

== Approach <approach>
How are you going to answer the research questions? For example:

A thesis proposal will be written at first (*@goal1*). For *#underline[@goal2]*, further steps are required. Here in the source code you can see how to link back to the goals A,B,...,N. Click on the corresponding goal-number text (#underline[@goal3], *@goal4*) to go to the respective goal in the list, and change the formatting as per your wish - see the source code.

\

= Plan
\
Write some 3-4 sentences on the order, but mainly focus on the time-table. Note that every planning time-table is subject to change and adaptation during the project. It is a first orientation, but should be discussed and adapted regularly.



#timeliney.timeline(
  show-grid: true,
  {
    import timeliney: *
      
    headerline(group([*Jun*],[*Jul*],[*Aug*],[*Sept*],[*Oct*],[*Nov*]))
    
    task("Literature review", (0, 1), style: (stroke: 2pt + gray))
    task("Writing Proposal (Goal A)", (0.5, 2), style: (stroke: 2pt + gray))
    task("Thesis work (Goal B)", (1, 4), style: (stroke: 2pt + gray))
    task("Buffer / Review / Stretch goals (Goal C-D)", (4,6), style: (stroke: 2pt + gray))

    milestone(
      at: 4,
      style: (stroke: (dash: "dashed")),
      align(center, [
        *Main goal completion*\
        Sept 2025
      ])
    )
  }
)




#line(length: 100%, stroke: gray)

#bibliography("refs.bib", style: "american-psychological-association")