#!/bin/bash

for file in inputs/*
do
if [ $(gfastats $file | grep 'Total scaffold length: ' | tr -d 'Total scaffold length: ') == 0 ]; then
echo $file
fi
done
