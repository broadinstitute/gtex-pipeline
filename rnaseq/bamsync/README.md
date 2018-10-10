### bamsync
This utility adds the following information that is lost when realigning BAM files:
* Copies the "read fails platform/vendor quality checks" flag (0x200) from the original BAM
* Adds a read group ID to each read based on the read ID, in the format RG:Z:xxxxx.x

Command:
```
bamsync reference.bam target.bam
```
