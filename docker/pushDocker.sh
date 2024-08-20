#!/bin/bash
#
TAG=""
cp ../manuscript/full_analysis/analysis.py .
docker build --progress=plain -t fcunial/kanpig_experiments .
docker push fcunial/kanpig_experiments${TAG}
rm -f analysis.py
