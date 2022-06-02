#! /bin/bash

set -e

git clone git@github.com:tanaylab/MCView

clean_up() {
    ARG=$?
    rm -rf MCView
    exit $ARG
} 
trap clean_up EXIT

pushd MCView

docker build --no-cache -f Dockerfile -t tanaylab/mcview:latest .
docker push tanaylab/mcview:latest

popd # MCView