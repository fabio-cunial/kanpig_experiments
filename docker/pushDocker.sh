#!/bin/bash
#
TAG=""
docker build --progress=plain -t fcunial/kanpig_experiments .
docker push fcunial/kanpig_experiments${TAG}
