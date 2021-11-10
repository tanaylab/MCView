#! /bin/bash

set -e

docker build -f Dockerfile -t tanaylab/mcview:latest .

docker push tanaylab/mcview:latest