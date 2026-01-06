#!/bin/bash

awk '/^>/ {
    sub(/^>/, "", $0);
    split($0, a, "|");
    print ">" a[1]
    next
}1' /16S/16S.fasta > /16S/16S_fixed-header.fasta
