#!/bin/bash
Rscript -e "library('pkgdown')" -e "pkgdown::build_site()"
