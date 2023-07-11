#!/bin/bash

for list in *txt
do
	while read file; do
		aws s3 cp $file inputs/
	done <$list
done
