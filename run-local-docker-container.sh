#!/bin/bash

: ${IMAGE:=nicolasbock/qmd-progress:latest}

docker pull ${IMAGE}
docker run --interactive --tty --rm \
  --volume ${PWD}:/qmd-progress --workdir /qmd-progress \
  --user $(id --user):$(id --group) \
  ${IMAGE}
