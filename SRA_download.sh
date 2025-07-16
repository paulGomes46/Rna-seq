#!/bin/bash

seq=()  # Initialize the array



# Populate the array with filenames
for i in {28364867..28364870}; do
    seq+=("SRR$i")  # Add each filename to the array
done

# Iterate over the array and process each file
for file in "${seq[@]}"; do  
    echo "starting"
    fasterq-dump "$file" --progress
    tar -czvf "$file.tar.gz" "$file"
    rm "$file"
    echo "$file compressed"
done

