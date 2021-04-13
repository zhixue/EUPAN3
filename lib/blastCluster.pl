#!/usr/bin/perl 
use strict;
use warnings;

my $usage="$0 <origin.fa> <identity> <self-blast_output> <output_prefix>

$0 is used to remove redundant sequences from self-blastn output.
blastn should be run with \" -evalue 1e-5  -outfmt \"6 qseqid sseqid 
qlen slen length qstart qend sstart send pident evalue\" -max_target_seqs 1000\"

";

die $usage if @ARGV!=4;
my ($fa,$IDEN,$blast_out,$prefix)=@ARGV;

my %seq=readfa($fa);
my %seqL=readLength($fa);

my @cluster;
my %fflag;
my %pflag;

my $index=0;
my %now;
my $pre_q="";

open(IN,$blast_out) || die "Cannot open $blast_out!\n";
while(<IN>){
    chomp;
    my ($qname,$tname,$qlen,$tlen, $len,$qstart,$qend,$tstart,$tend,$iden,$evalue)=split /\t/,$_;
    next if $qname eq $tname;
    next if defined $fflag{$qname};
    $iden/=100;
    if($qname eq $pre_q){
	if ($qlen>=$tlen){
	    add_to_now($tname,$iden,$len,$tstart,$tend,$tlen);
	}
	else{
	    add_to_now($tname,$iden,$len,$qstart,$qend,$qlen);
	}
    }
    else{
	my $new_clus=init_item($qname);
	$fflag{$qname}=1;
	$pflag{$qname}=$index;
	foreach my $k (keys(%now)){
	    if($now{$k}->{iden}>$IDEN){
		if(defined $fflag{$k}){
		    next if $fflag{$k}>=$now{$k}->{iden};
		    $cluster[$pflag{$k}]->{$k}=0;
		    $new_clus->{$k}=1;
		    $fflag{$k}=$now{$k}->{iden};
		    $pflag{$k}=$index;
		}
		else{
		    $new_clus->{$k}=1;
		    $fflag{$k}=$now{$k}->{iden};
		    $pflag{$k}=$index;
		}
	    }
	}
	push @cluster, $new_clus;
	$index++;
	%now=();
	if ($qlen>=$tlen){
	    add_to_now($tname,$iden,$len,$tstart,$tend,$tlen);
	}
	else{
	    add_to_now($tname,$iden,$len,$qstart,$qend,$qlen);
	}
    }
    $pre_q=$qname;
}
close IN;

my $outcl=$prefix.".clust";
my $outfa=$prefix;
my %rep;

open(OUT,">$outcl");
my $i=0;
for($i=0;$i<@cluster;$i++){
    my $l=0;
    my $rep="";
    foreach my $k (keys(%{$cluster[$i]})){
	if($cluster[$i]->{$k}==1){
	    if($seqL{$k}>$l){
		$rep=$k;
		$l=$seqL{$k};
	    }
	}
    }
    $rep{$rep}=1;
    print OUT "Group ",$i+1,": ";
    foreach my $k (keys(%{$cluster[$i]})){
	print OUT " ",$k if $cluster[$i]->{$k}==1;
    }
    print OUT "\n";
}
foreach my $k (keys(%seq)){
    if(!defined($fflag{$k})){
	$rep{$k}=1;
	print OUT "Group $i: $k\n";
	$i++;
    }
}
close OUT;

open(OUT,">$outfa");
open(IN,$fa);
my $f=0;
while(<IN>){
    chomp;
    if(/^>(.+)$/){
	if(defined $rep{(split(" ",$1,2))[0]}){
	    $f=1;
	    print OUT $_,"\n";
	}
    }
    else{
	print OUT $_,"\n" if $f==1;
    }
}
close IN;
close OUT;

################################################
sub readfa{
    my %h;
    open(IN,$_[0]);
    while(<IN>){
	chomp;
	$h{(split(" ",$1,2))[0]}=1 if(/^>(.+)$/);
    }
    close IN;
    return %h;
}

sub readLength{
    my %h;
    open(IN,$_[0]);
    my $name="";
    my $idn="";
    while(<IN>){
	chomp;
	if(/^>(.+)$/){
	    $name=$1;
	    $idn=(split(" ",$1,2))[0];
	    $h{$idn}=0;
	}
	else{
	    $h{$idn}+=length($_);
	}
    }
    close IN;
    return %h;
}

sub init_item{
    my %h;
    $h{$_[0]}=1;
    return \%h;
}

sub add_to_now{
    my ($tname,$iden,$len,$tstart,$tend,$tlen)=@_;
    if(defined $now{$tname}){
	if(!check_overlap($now{$tname},$tstart,$tend)){
	    push @{$now{$tname}->{start}},$tstart;
	    push @{$now{$tname}->{end}},$tend;
	    $now{$tname}->{iden}=$now{$tname}->{iden}+$iden*$len/$tlen;
	}
    }
    else{
	$now{$tname}=init_h($iden*$len/$tlen,$tstart,$tend);
    }
}

sub init_h{
    my %h;
    $h{iden}=$_[0];
    $h{start}=init_array($_[1]);
    $h{end}=init_array($_[2]);
    return \%h;
}

sub init_array{
    my @a;
    $a[0]=$_[0];
    return \@a;
}
sub check_overlap{
    my ($ad,$s,$e)=@_;
    my $f=0;
    for(my $i=0;$i<@{$ad->{start}};$i++){
	if(check_overlap_helper($s,$e,$ad->{start}->[$i],$ad->{end}->$i)){
	    $f=1;
	    last;
	}
    }
    close $f;
}

sub check_overlap_helper{
    my ($a,$b,$c,$d)=@_;
    if($a<$c){
	return 0 if $b<$c;
	return 1 if $b>=$c;
    }
    elsif($a>$c){
	return 0 if $d<$a;
	return 1 if $d>=$a;
    }
    else{
	return 1;
    }
}
