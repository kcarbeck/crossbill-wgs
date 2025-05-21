#!/usr/bin/env bash

SOURCE_DIR="/workdir/kcarbeck/sorted"
DEST_DIR="/lustre2/home/lc736_0001/crossbill/sortedBam"

# We'll record size and path for each file in source and destination.
# Then compare sorted lists.

# 1. Gather file sizes and paths for source
find "$SOURCE_DIR" -type f -printf "%s %p\n" \
    | sed "s|$SOURCE_DIR/||" \
    | sort > /workdir/kcarbeck/tmp/source_sizes.txt

# 2. Gather file sizes and paths for destination
find "$DEST_DIR" -type f -printf "%s %p\n" \
    | sed "s|$DEST_DIR/||" \
    | sort > /workdir/kcarbeck/tmp/dest_sizes.txt

# 3. Compare the two size lists
diff /tmp/source_sizes.txt /tmp/dest_sizes.txt
if [ $? -eq 0 ]; then
    echo "All files match in size and path."
else
    echo "File size mismatch found (or missing files)."
fi