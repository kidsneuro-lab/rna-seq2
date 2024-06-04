#!/bin/sh

echo "Generating MultiQC report"

which multiqc 1> /dev/null 2> /dev/null

if [ $? -eq 0 ]; then
    multiqc \
      --force \
      --outdir results/multiqc_report \
      --dirs \
      results
else
    echo "MultiQC not available"
    exit 1
fi

