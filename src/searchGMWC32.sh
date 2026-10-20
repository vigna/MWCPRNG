#!/bin/bash

# Compile searchGMWC with the desired BITS specification

if [[ "$3" == "" ]]; then
	echo "$(basename $0) SEED MASKBITS MAXDIM" 1>&2
	exit 1
fi

MAXDIM=$3
FILE=$(mktemp)

./searchGMWC $1 $2 | while read M MINUSA0 A0INV A1; do ./printdat 1 $MAXDIM 4294967296 $M $FILE; echo -ne "$M\t$MINUSA0\t$A0INV\t$A1\t"; lattestZZRR $(for((i=2; i<=$MAXDIM; i++)); do echo -n " $FILE-$i.dat"; done) | grep Figure | egrep -o "[0-9]+[.]?[0-9e-]+" | tr '\n' '\t'; echo ; done | python3 <(cat <<EOF
import sys

for line in sys.stdin:
    a = line.split('\t')
    m = int(a[0])
    minus_a0 = int(a[1])
    a0inv = int(a[2])
    a1 = int(a[3])

    sys.stdout.write("0x%x\t0x%x\t0x%x\t0x%x\t%s" % (m, minus_a0, a0inv, a1, '\t'.join(a[4:])))
    sys.stdout.flush()
EOF
)
