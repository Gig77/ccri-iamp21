use strict;
use warnings;

my $outfile = "/mnt/projects/iamp/results/geneset_annotation.tsv";
open(OUT, ">$outfile") or die "ERROR: Could not write to output file $outfile.\n";

print OUT "geneset\tcategory\tdescription\tlinkid\n";

# MSigDB 5.0

open(IN, "/mnt/projects/generic/data/msigdb5.0/h.all.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_H_HALLMARK\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c1.all.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C1_POSITIONAL\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c2.cgp.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C2_CGP\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c2.cp.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C2_PATHWAYS\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c3.mir.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C3_MOTIF_MIR\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c3.tft.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C3_MOTIF_TFT\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c4.cgn.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C4_COMPUTATIONAL_CGN\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c4.cm.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C4_COMPUTATIONAL_CM\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c5.bp.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C5_GO_BP\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c5.cc.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C5_GO_CC\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c5.mf.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C5_GO_MF\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c6.all.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C6_ONCOGENIC_SIGNATURES\t\t$gs\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c7.all.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C7_IMMUNOLOGIC_SIGNATURES\t\t$gs\n";
}
close(IN);

# DSigDB 1.0

my %dsigdb;
open(IN, "/mnt/projects/generic/data/DSigDB/DSigDB_All.txt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	if (/^Compound\(([^\)]+)\) : (.*)/) {
		my ($gs, $cat) = ($2, $1);
		$cat =~ s/ /_/g;
		
		if ($cat eq "TTD") {
			$cat = "D4_TTD";
		}
		elsif ($cat eq "CTD") {
			$cat = "D4_CTD";
		}
		elsif ($cat eq "BOSS") {
			$cat = "D4_BOSS";
		}
		elsif ($cat eq "D1") {
			$cat = "D1_FDA";
		}
		elsif ($cat eq "D3") {
			$cat = "D3_CMAP_PERTURBATION";
		}
		$dsigdb{$gs} = $cat;
	}
#	elsif (/^Compound : ([^\(]+)\(([^\)]+)\)/) {
#		my ($gs, $cat) = ($1, $2);
#		$cat =~ s/ /_/g;
#		$dsigdb{$gs} = $cat;
#	}
}
close(IN);

open(IN, "/mnt/projects/generic/data/DSigDB/DSigDB_v1.0_All.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	next if (!$gs);

	my $gs_short = $gs;
	$gs_short =~ s/_TTD_\d+//;
	$gs_short =~ s/_CTD_\d+//;
	
	if ($gs =~ /(.*)_MRC$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_MRC\t\t$gs\n";
		next;
	}
	elsif ($gs =~ /(.*)_LINCS$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_LINCS\t\t$gs\n";
		next;
	}
	elsif ($gs =~ /(.*)_Roche$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_ROCHE\t\t$gs\n";
		next;
	}
	elsif ($gs =~ /(.*)_GSK$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_GSK\t\t$gs\n";
		next;
	}
	elsif ($gs =~ /(.*)_FDA$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_FDA\t\t$gs\n";
		next;
	}
	elsif ($gs =~ /(.*)_RBC$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_RBC\t\t$gs\n";
		next;
	}
	elsif ($gs =~ /(.*)_Kinome Scan$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_KINOMESCAN\t\t$gs\n";
		next;
	}
	elsif ($gs =~ /(.*)_BOSS$/) {
		print OUT uc($gs),"\tDSIGDB_D4_BOSS\t\t$gs\n";
		next;
	}
	elsif ($dsigdb{$gs_short}) {
		print OUT uc($gs),"\tDSIGDB_",$dsigdb{$gs_short}, "\t\t$gs\n";
	}
	else {
		print STDERR "WARNING: Category for DSigDB gene set $gs unknown.\n";
	}
}
close(IN);

# GeneSigDB 4.0

open(IN, "/mnt/projects/generic/data/GeneSigDB/ALL_SIGSv4.nodup.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	next if (!$gs);

	my ($pmid, $rest) = $gs =~ /([^-_]+)[-_](.*)/;
	my $xmlfile = "/mnt/projects/generic/data/GeneSigDB/details/$pmid/$gs-index.xml";
	my $desc = "";
	open(F, $xmlfile) or die "ERROR: Could not open file $xmlfile\n";
	while (<F>) {
		if (/<sigDescription>(.*)<\/sigDescription>/) {
			$desc = $1;
			last;
		}
	}
	close(F);
	print OUT uc($gs),"\tGeneSigDB\t$desc\t$gs\n";
}
close(IN);

close(OUT);

print "Success. Output written to $outfile\n"; 

