#!/bin/sh

# Shell script for programatic access to Revigo. Run it with:
# revigo.sh example.csv

# Submit job to Revigo and fetch results
curl "http://revigo.irb.hr/Revigo" -X POST --silent --data-urlencode "cutoff=$2" --data-urlencode "valueType=pvalue" --data-urlencode "speciesTaxon=0" --data-urlencode "measure=SIMREL" --data-urlencode "goList@$1" --data-urlencode "removeObsolete=$3" --header "Content-Type: application/x-www-form-urlencoded"
