#!/bin/sh
# Write SHA256SUMS for a tag's ordering files.
#
#   tools/make_sums.sh v1 data/
#
# The tag goes in a leading "# <tag>" comment, NOT a third column: a third
# column would be parsed as part of the filename and break `shasum -c`, losing
# the ability to verify a store by hand with standard tools. A comment keeps the
# file both self-describing and checkable:
#
#   $ shasum -a 256 -c SHA256SUMS
#   EPIC.ordering.tsv.gz: OK
set -eu
tag=${1:-v1}
dir=${2:-data}
cd "$dir"
{
    echo "# $tag"
    LC_ALL=C ls *.ordering.tsv.gz | LC_ALL=C sort | xargs shasum -a 256
} > SHA256SUMS.new
mv SHA256SUMS.new SHA256SUMS
echo "wrote $dir/SHA256SUMS for $tag"
