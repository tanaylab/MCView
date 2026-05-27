#! /bin/bash

set -e

# Builds and publishes the MCView Docker image to Docker Hub, tagged with both
# the package version (from DESCRIPTION) and `latest`.
#
# Requires a reasonably current Docker (>= 20.10.10) -- see the note in Dockerfile.

git clone git@github.com:tanaylab/MCView

clean_up() {
    ARG=$?
    rm -rf MCView
    exit $ARG
}
trap clean_up EXIT

pushd MCView

VERSION=$(grep -i '^Version:' DESCRIPTION | awk '{print $2}')
GIT_SHA=$(git rev-parse HEAD)
echo "Building tanaylab/mcview:${VERSION} (commit ${GIT_SHA}) and :latest"

# Pin the image to the exact cloned commit so it matches the version tag.
docker build --no-cache -f Dockerfile \
    --build-arg "MCVIEW_REF=${GIT_SHA}" \
    -t "tanaylab/mcview:${VERSION}" \
    -t tanaylab/mcview:latest .

docker push "tanaylab/mcview:${VERSION}"
docker push tanaylab/mcview:latest

popd # MCView
