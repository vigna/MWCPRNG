#!/bin/bash

# Compile searchMWC with the desired BITS specification

if [[ "$2" == "" ]]; then
	echo "$(basename $0) SEED MAXDIM" 1>&2
	exit 1
fi

MAXDIM=$2
FILE=$(mktemp)

./searchMWC $1 | while read M A1; do ./printdat 1 $MAXDIM 18446744073709551616 $M $FILE; echo -ne "$M\t$A1\t"; lattestZZRR $(for((i=2; i<=$MAXDIM; i++)); do echo -n " $FILE-$i.dat"; done) | grep Figure | egrep -o "[0-9]+[.]?[0-9e-]+" | tr '\n' '\t'; echo ; done | python3 <(cat <<EOF
import sys

for line in sys.stdin:
    a = line.split('\t')
    m = int(a[0])
    a1 = int(a[1])

    sys.stdout.write("0x%x\t0x%x\t%s" % (m, a1, '\t'.join(a[2:])))
    sys.stdout.flush()
EOF
)
