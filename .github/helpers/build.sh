#!/bin/bash
R -e "library('devtools')" -e "devtools::build(binary = TRUE, args = c('--preclean'))"
