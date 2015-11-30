use warnings;
use strict;

foreach my $e (`tabix /mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz chr1:34047-34047`) {
	print $e;
}
