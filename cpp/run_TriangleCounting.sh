#!/bin/bash -x

if [ $# -ne 1 ]; then
    echo "USAGE: run_TriangleCounting.sh [Dataset]"
    exit 1
fi

./TriangleCounting ../data/${1}/edges.csv 10000 1-0.001 2 6-150 10-1 2
./TriangleCounting ../data/${1}/edges.csv 10000 1-0.001 2 6-150 10-1 3
./TriangleCounting ../data/${1}/edges.csv 10000 1-0.001 2 6-150 10-1 4
./TriangleCounting ../data/${1}/edges.csv 10000 1-0.001 0 0 10-1 2
./TriangleCounting ../data/${1}/edges.csv 10000 1-0.001 0 0 10-1 3
./TriangleCounting ../data/${1}/edges.csv 10000 1-0.001 0 0 10-1 4
