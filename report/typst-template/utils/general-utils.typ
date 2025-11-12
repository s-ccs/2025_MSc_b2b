#import "@preview/big-todo:0.2.0": *
#import "@preview/timeliney:0.2.1"

// prose citation helper function (from https://github.com/typst/typst/issues/2716#issuecomment-1817870741)
// To use: #pc[@refname]
#let pc = (citation) => {
  set cite(form: "prose")
  citation
}

// inline todo helper function
// usage: #ilt("todo text") or #ilt[todo text]
#let ilt = (text) => {
  todo(text, inline: true)
}

// make it possible to click and refer back to goals by creating a custom figure type 
#let goal = (text) => {
  figure(
    [
      #text
    ],
    supplement: [Goal],
    kind: "goal",
    outlined: false,
    numbering: "A",
    // caption: [#text]
  )
}