context("collapse orthogroups")

test_that("collapse orthogroup requires expression and ortholog data", {
          expect_error(collapse_orthologs(make_ogset()))
          })
