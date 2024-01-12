#!/usr/bin/env bash


base="/local/storage/liang/Downloads/Datasets/polyline_wkt"

./build/lsi "$base/dtl_cnty/dtl_cnty.wkt" "$base/USAZIPCodeArea/USAZIPCodeArea.wkt" |& tee "dtl_cnty_USAZIPCodeArea.log"
./build/lsi "$base/USACensusBlockGroupBoundaries/USACensusBlockGroupBoundaries.wkt" "$base/USADetailedWaterBodies/USADetailedWaterBodies.wkt" |& tee "USACensusBlockGroupBoundaries_USADetailedWaterBodies.log"

for cont in Africa Asia Australia Europe North_America South_America; do
  ./build/lsi "$base/lakes/$cont/lakes_${cont}.wkt" "$base/parks/$cont/parks_${cont}.wkt" |& tee "lakes_parks_${cont}.log"
done
