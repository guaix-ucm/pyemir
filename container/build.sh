#!/bin/bash
# Build the container image of pyemir locally with podman
#
# Run it from the root of the repository. The image is tagged with the
# version of pyemir and with latest, for ghcr.io and quay.io. The
# release workflow publishes the image in both registries.

set -euo pipefail

IMAGE_NAMES=("ghcr.io/guaix-ucm/pyemir" "quay.io/guaix-ucm/pyemir")

python -m build --wheel --outdir dist
WHEEL=$(ls -t dist/pyemir-*.whl | head -n 1)
# The version is in the name of the wheel: pyemir-VERSION-py3-none-any.whl
VERSION=$(basename "$WHEEL" | cut -d- -f2)
# '+' (local versions of development) is not valid in a tag
TAG=${VERSION//+/-}

TAGS=()
for name in "${IMAGE_NAMES[@]}"; do
    TAGS+=(-t "$name:$TAG" -t "$name:latest")
done

podman build --format docker \
    -f container/Containerfile \
    --label org.opencontainers.image.version="$VERSION" \
    --label org.opencontainers.image.created="$(date -Iseconds)" \
    "${TAGS[@]}" \
    .
