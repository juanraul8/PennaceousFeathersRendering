#!/bin/bash

for file in ../bin/*.svg; do
    [ -f "$file" ] || break
    inkscape --export-type="png" "$file"
done
