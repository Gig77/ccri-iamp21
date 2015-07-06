use warnings;
use strict;

my $jaspar_motifs = "/mnt/projects/generic/data/opossum3/jaspar_motifs.txt";
my $gmt = "/mnt/projects/generic/data/opossum3/jaspar_core.gmt";

my $cmd = "mysql --user=opossum_r --host=opossum.cmmt.ubc.ca --database=oPOSSUM31_human --column-names=FALSE --connect-timeout=10000 -B -A -e \"select distinct tf_id from conserved_tfbss where tf_id like 'MA%';\" > $jaspar_motifs";
print("$cmd\n");
system($cmd);

unlink($gmt);

my $motifs = 0;
open(M, "$jaspar_motifs") or die "$!\n";
while(<M>) {
	chomp;

	$cmd = "mysql --user=opossum_r --host=opossum.cmmt.ubc.ca --database=oPOSSUM31_human --column-names=FALSE --connect-timeout=10000 -B -A -e 'SET SESSION group_concat_max_len = 1000000; select tf_id, CONCAT(\"http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=\", tf_id, \"&rm=present&collection=CORE\"), GROUP_CONCAT(external_id) from (select distinct conserved_tfbss.tf_id, external_gene_ids.external_id from conserved_tfbss join external_gene_ids on external_gene_ids.gene_id=conserved_tfbss.gene_id join external_gene_id_types on external_gene_id_types.id_type=external_gene_ids.id_type join genes on genes.gene_id=conserved_tfbss.gene_id where external_gene_id_types.name=\"HGNC symbol\" and genes.biotype=\"protein_coding\" and conserved_tfbss.conservation >= 0.60 and conserved_tfbss.rel_score >= 0.85 and conserved_tfbss.start <= 2000 and conserved_tfbss.tf_id=\"$_\") x group by tf_id;' | sed 's/,/\\t/g' >> $gmt.part";
	print("$cmd\n");
	system($cmd);

	$motifs ++;
}
close(M);

print("$motifs motifs written to file $gmt.part\n");

# translate Jaspar motif ID to TF name
print("Translating Jaspar motif names to TF names...\n");
my %motif2tf;
open(J, "/mnt/projects/generic/data/jaspar-5.0alpha/MATRIX.txt");
while(<J>) {
	chomp;
	my ($id, $db, $motif, $version, $tf) = split("\t");
	$motif2tf{$motif} = $tf;
} 
close(J);

open(OUT, ">$gmt") or die "Could not write to file $gmt\n";
open(GS, "$gmt.part") or die "Could not read $gmt.part\n";
while(<GS>) {
	my ($motif, $rest) = /^([^\.\t]+)\.?\d*(.*)/;
	$motif = $motif2tf{$motif}."_JASPAR" if ($motif2tf{$motif});
	print OUT $motif,$rest,"\n";
}
close(GS);
close(OUT);

unlink("$gmt.part");

print("Finished with success.\n");