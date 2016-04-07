#!/usr/bin/perl -w

my $usage = << "EOF";

Usage: $0
	-l promoter length
	-norg number of organisms
	-ncyc number of cycles
        -m matrix file
	-c motif cutoff
	-i list of tfs
        -o output file
EOF

use Getopt::Long;

my %par = (
                  );

&GetOptions("l=i" => \$par{length},
	    "norg=i" => \$par{norg},
	    "ncyc=i" => \$par{ncyc},
	    "m=s" => \$par{matrix},
	    "c=f" => \$par{cutoff},
	    "i=s" => \$par{infile},
            "o=s" => \$par{outfile});

($par{length} && $par{norg} && $par{ncyc} && $par{matrix} && $par{cutoff} && $par{infile} && $par{outfile}) or die $usage;

&go (\%par );

print STDERR "$0 done\n";

###############################################

sub go {
	my ( $par ) = @_;

    	my @motlens = (0,18,20,12,13,7,17,6,5,6,10,9,12,7,8,10,7,12,12,9,12,8,7,12,12);
	my @maxscores = (0,18.287,25.879,11.403,13.032,8.717,15.523,7.175,5.979,7.503,10.639,10.881,11.78,7.429,10.28,9.38,8.782,12.546,12.189,9.841,10.802,8.888,8.164,11.466,11.299);
	my $maxlen = 20;
	my ($mx, $motnames) = &init_mx($par->{matrix});

	open my $handle, '<', $par->{infile};
	chomp(my @tfs = <$handle>);
	close $handle;

	my $ntfs = $#tfs;
	my $organisms = init($par->{length},$par->{norg},$par->{infile},\@tfs,$mx,\@motlens);

    	open (OF, ">$par->{outfile}") or die "Cannot open outfile $par->{outfile} $!\n";

	my $date = localtime;

	print OF "Date: ".$date."\nMotifs:";
	for $ii (0 .. $ntfs) {
		my $ii2 = $motnames->{$tfs[$ii]};
		print OF "\t$ii2";
	}
	print OF "\nPromoter length: ".$par->{length}."\n";
        print OF "Number of organisms: ".$par->{norg}."\n";
        print OF "Number of generations: ".$par->{ncyc}."\n";
        print OF "Motif cutoff value: ".$par->{cutoff}."\n\n-------------------------\n\nGeneration";
	for my $ij (1 .. $ntfs+1) {
		print OF "\tMotif ".$ij."\tScore\tPosition";
	}
	print OF "\n";

	for my $cyc (1 .. $par->{ncyc}) {
		# mutate promoters
		&mutate($organisms, $par->{norg}, $par->{length});
		# check if improved
		my $ok = &match_score($organisms,$par->{norg},\@tfs,$mx,\@motlens,\@maxscores,$par->{cutoff});
		if ($ok == 1) {
			print OF "Number of generations needed for module formation: $cyc\n";
			exit(-1);
		}
		# write motifs to output
		&write_gen_info($organisms, \*OF, $ntfs, $cyc, \@tfs, \@motlens);
	}

	print OF "Number of generations needed for module formation: $par->{ncyc}\n";

    	close (OF);
}

sub init_mx {
        my $mxf = shift;
        my %mx = ();
        my %motnames = ();

        my $nrow = 0;
        my $n_motif = 1;
        my %bp = (1,"A",2,"T",3,"G",0,"C");

        open (MF, "<$mxf")  or die "Cannot open matrix file $mxf $!\n";
        while (<MF>) {
                chomp;
		if ($_ =~ /motif/) {
			my ($m, $mname, $mm) = split;
			$motnames{$n_motif} = $mname;
		}
                unless ($_ =~ /motif/) {
			$nrow++;
  	              my ($b, @w) = split;
        	        my $h = $nrow % 4;
        	        for my $i (0 .. $#w) {
        	                my $j = $i + 1;
        	                $mx{$n_motif}{$bp{$h}}{$j} = $w[$i];
        	        }
        	        if ($h == 0) {
        	                $n_motif++;
        	        }
		}
        }
        close (MF);

        return (\%mx,\%motnames);
}

sub init {
	my ($len, $norg, $inf, $tfs, $mx, $motlens) = @_;
	my %orgs = ();
	for my $i (1 .. $norg) {
		$orgs{$i}{mutpos} = 0;
		my $prom = &generate_promoter($len);
		$orgs{$i}{promoter} = $prom;
		for my $j (0 .. $#{$tfs}) {
			my $motid = $tfs->[$j];
			my ($max_score, $max_pos) = &init_match_score($prom, $motid, $mx, $motlens);
			$orgs{$i}{$j}{score} = $max_score;
			$orgs{$i}{$j}{pos} = $max_pos;
			$orgs{$i}{score} += $max_score;
		}
	}
	return \%orgs;
}

sub generate_promoter {
	my $len = shift;
	my $prom = "";
	my %bp = (1,"A",2,"T",3,"G",0,"C");
	for my $i (1 .. $len) {
		my $rn = int(rand(4));
		$prom .= $bp{$rn};
	}
	return $prom;
}

sub init_match_score {
	my ($prom, $motid, $mx, $motlens) = @_;
	my $maxscore = -1000;
	my $maxpos = 0;
	my $motlen = $motlens->[$motid];
	for my $bp (0 .. length($prom)-$motlen-1) {
		my $prom_slice = substr($prom,$bp,$motlen);
		my $score = 0;
		for my $i (1 .. $motlen) {
			my $b = substr($prom_slice,$i-1,1);
			my $val = $mx->{$motid}{$b}{$i};
			$score += $val;
		}
		if ($score > $maxscore) {
			$maxscore = $score;
			$maxpos = $bp;
		}
	}
	return ($maxscore, $maxpos);
}

sub match_score {
	my ($orgs, $norg, $tfs, $mx, $motlens, $maxscores, $cutoff) = @_;
	my $ok = 0;
	for my $o (1 .. $norg) {
		$orgs->{$o}{score} = 0;
		my $ok2 = 0;
		my $prom = $orgs->{$o}{promoter};
		my $prom_center = substr($prom,($orgs->{$o}{mutpos})-20,2*20+1);
		for my $t (0 .. $#{$tfs}) {
			my $maxscore = $orgs->{$o}{$t}{score};
			my $maxpos = 0;
			my $k = $t + 1;
			my $motid = $tfs->[$t];
			my $motlen = $motlens->[$motid];
			for my $bp (0 .. length($prom_center)-$motlen-1) {
				my $prom_slice = substr($prom_center,$bp,$motlen);
				my $score = 0;
				for my $i (1 .. $motlen) {
		                        my $b = substr($prom_slice,$i-1,1);
        		                my $val = $mx->{$motid}{$b}{$i};
                		        $score += $val;
				}
		                if ($score > $maxscore) {
                		        $maxscore = $score;
					$orgs->{$o}{$t}{score} = $score;
					$maxpos = $orgs->{$o}{mutpos} + $bp - 20;
                		}
			}
			my $mt = substr($orgs->{$o}{promoter},$orgs->{$o}{$t}{pos},$motlen);
			if ($maxscore/($maxscores->[$motid]) >= $cutoff) {
				$ok2++;
			}
		}
		if ($ok2 == $#{$tfs}+1) {$ok = 1;}
		for my $t2 (0 .. $#{$tfs}) {
			 $orgs->{$o}{score} += $orgs->{$o}{$t2}{score};
		}
	}
	return $ok;
}

sub mutate {
	my ($org, $norg, $len) = @_;
	my %bp = (1,"A",2,"T",3,"G",0,"C");
	for my $o (1 .. $norg) {
		my $prom = $org->{$o}{promoter};
		my $rn = int(rand($len))+1;
		my $bp = substr($prom,$rn,1);
		my $rn2 = int(rand(4));
		my $rbp = $bp{$rn2};
		while ($rbp eq $bp) {
			$rn2 = int(rand(4));
			$rbp = $bp{$rn2};
		}
		substr($prom,$rn,1) = $rbp;
		$org->{$o}{mutpos} = $rn;
		$org->{$o}{promoter} = $prom;
	}
}

sub write_gen_info {
	my ($organisms, $of, $ntfs, $cyc, $tfs, $motlens) = @_;
	print $of $cyc;
	for my $i (0 .. $ntfs) {
		my $j = $i + 1;
		my $tf_ind = $tfs->[$i];
		my $motlen = $motlens->[$tf_ind];
		my $motpos = $organisms->{1}{$i}{pos};
		my $motscore = $organisms->{1}{$i}{score};
		my $mot_occ = substr($organisms->{1}{promoter},$motpos,$motlen);
		print $of "\t$mot_occ\t$motscore\t$motpos";
	}
	print $of "\n";
}
