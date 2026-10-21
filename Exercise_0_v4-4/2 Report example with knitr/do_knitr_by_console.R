library(knitr)
# generate the tex file from the Rnw file:
knit(input = "knitr_example.Rnw", output = "knitr_example.tex")

# generate the PDF file from the Rnw file:
knit2pdf(input = "knitr_example.Rnw", output = "knitr_example.tex")

# extract R code:
purl("knitr_example.Rnw")
