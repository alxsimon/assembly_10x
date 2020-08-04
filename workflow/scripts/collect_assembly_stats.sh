#!/usr/bin/env bash

echo "species,version,$(head -n1 results/supernova_assemblies/gallo_v1/outs/summary.csv)" \
> results/assemblies_stats.csv

for S in 'gallo' 'edu' 'tros'; do
	for V in 'v1' 'v2'; do
		echo "${S},${V},$(tail -n +2 results/supernova_assemblies/${S}_${V}/outs/summary.csv)" \
		>> results/assemblies_stats.csv
	done
done

