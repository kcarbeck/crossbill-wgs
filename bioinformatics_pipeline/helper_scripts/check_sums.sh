# 1. Generate checksums for /path/to/source
find /workdir/kcarbeck/sorted -type f -print0 \
    | xargs -0 -n1 -P 4 sha256sum \
    > /tmp/source_sums.txt

# 2. Generate checksums for /path/to/destination
find /lustre2/home/lc736_0001/crossbill/sortedBam -type f -print0 \
    | xargs -0 -n1 -P 4 sha256sum \
    > /tmp/destination_sums.txt

# -type f only finds files.

# -print0 combined with -0 in xargs helps handle filenames with spaces and special characters properly.

# -n1 makes sure that xargs spawns one sha256sum process per file.

# -P 4 indicates the number of processes you want to run in parallel (tweak as needed for your CPU).


# For /path/to/source checksums, remove /path/to/source prefix
sed "s|/path/to/source/||" /tmp/source_sums.txt \
    > /tmp/source_sums_rel.txt

# For /path/to/destination checksums, remove /path/to/destination prefix
sed "s|/path/to/destination/||" /tmp/destination_sums.txt \
    > /tmp/destination_sums_rel.txt

# Sort and compare
sort /tmp/source_sums_rel.txt -o /tmp/source_sums_rel.txt
sort /tmp/destination_sums_rel.txt -o /tmp/destination_sums_rel.txt

diff /tmp/source_sums_rel.txt /tmp/destination_sums_rel.txt
