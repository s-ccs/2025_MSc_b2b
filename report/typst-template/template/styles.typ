#let styles = rest => {
  // show heading.where(level: 1): it => {
  //   v(10pt)
  //   text(it, size: 1em)
  //   // v(14pt)
  // }
  
  show heading.where(level: 2): it => {
      v(5pt, weak: false)
      text(it)//, size: 1.2em)
      numbering("1")
      v(0pt, weak: true)
  }
  
  set par(justify: true)

  set page(numbering: "1")
  
  rest
}