#!/usr/bin/env bash
set -e

docker build -t sesame-workflow-shiny .
docker run --rm -p 3838:3838 sesame-workflow-shiny
