package Snasm;

use 5.26.0;
use strict;
use warnings;
use File::Which; 
use File::Basename;
use Cwd qw(abs_path);
use JSON::PP;
#use Path::Tiny;
#use Hash::Merge qw(merge);
use Data::Dumper;
use Getopt::Std;
use File::Basename;
use Time::Piece;
use File::Copy;
use File::Temp;
use File::Path qw(make_path);
use List::Util qw(min max);
use Digest::MD5 qw(md5_hex);

use Exporter 'import';
our @EXPORT = qw(
 msg wrn err
 tsv csv 
 tsv2hash tsv2array
 csv2hash csv2array
 array2tsv hash2tsv
 run exe cpus today need
 fixfile fixdir fixdate fixint fixcvg
 fhandle fh2file tmpdir
 get_opts show_usage
 slurp slurp_fh spew spew_nl
 kovid_dir kovid_db 
 kovid_db_file kovid_seq_file kovid_csv_file
 kovid_ref_fasta kovid_ref_id kovid_ref_len
 kovid_enter kovid_exit
 fasta2array array2fasta 
 afa2array afa2hash
 fasta2hash hash2fasta
 subseq revcom translate
 base_tally site_tally
 hash_bin
);

our %CODON = (
  "TTT" => "F", "TTC" => "F", "TTA" => "L",
  "TTG" => "L", "CTT" => "L", "CTC" => "L",
  "CTA" => "L", "CTG" => "L", "ATT" => "I",
  "ATC" => "I", "ATA" => "I", "ATG" => "M",
  "GTT" => "V", "GTC" => "V", "GTA" => "V",
  "GTG" => "v", "TCT" => "S", "TCC" => "S",
  "TCA" => "S", "TCG" => "S", "CCT" => "P",
  "CCC" => "P", "CCA" => "P", "CCG" => "P",
  "ACT" => "T", "ACC" => "T", "ACA" => "T",
  "ACG" => "T", "GCT" => "A", "GCC" => "A",
  "GCA" => "A", "GCG" => "A", "TAT" => "Y",
  "TAC" => "Y", "TAA" => "*", "TAG" => "*",
  "CAT" => "H", "CAC" => "H", "CAA" => "Q",
  "CAG" => "Q", "AAT" => "N", "AAC" => "N",
  "AAA" => "K","AAG" => "K" ,"GAT" => "D",
  "GAC" => "D", "GAA" => "E", "GAG" => "E", 
  "TGT" => "C", "TGC" => "C", "TGA" => "*",
  "TGG" => "W", "CGT" => "R", "CGC" => "R",
  "CGA" => "R", "CGG" => "R", "AGT" => "S", 
  "AGC" => "S", "AGA" => "R", "AGG" => "R",
  "GGT" => "G", "GGC" => "G", "GGA" => "G",
  "GGG" => "G",
);

# used by kovid_enter/exit and tempdir
my $old_dir;

#....................................................................
sub msg { say STDERR "@_"; }

#....................................................................
sub err { msg("ERROR:", @_); exit(1); }

#....................................................................
sub wrn { msg("WARNING:", @_); }

#....................................................................
sub tsv { join("\t", @_)."\n" }

#....................................................................
sub kovid_enter {
  msg("kovid toolkit by \@torstenseemann");
  msg(@_) if @_;
}

#....................................................................
sub kovid_exit {
  #chdir("/"); # allow tmpdir removal
  if ($old_dir) {
    msg("Changing back to dir: $old_dir");
    chdir $old_dir;
  }
  msg(@_) if @_;
  msg("Done.");
  exit(0);
}

#....................................................................
sub tsv2hash {
  my($fname, $keycol, $sep) = @_;
  $sep ||= "\t";
  msg("Opening tabular file: $fname");
  open my $TABLE, '<', $fname;
  my @hdr;
  my %hash;
  my $count=0;
  while (<$TABLE>) {
    chomp;
    my @row = split m/$sep/;
    if (@hdr) {
      my %data;
      @data{@hdr} = @row; # make hash
      $count++;
      my $key = $data{$keycol} // "_${keycol}_$count";
      $hash{$key} = \%data;
    }
    else {
      @hdr = @row;
      msg("Found", 0+@hdr, "columns:",
          join(' | ', @hdr) );
    }
  }
  close $TABLE;
  msg("Loaded $count data rows.");
  return \%hash;
}

#....................................................................
sub csv2hash {
  my(@opt) = @_;
  my $tsv = File::Temp->new();
  run("csvtk csv2tab \Q$opt[0]\E > $tsv");
  $opt[0] = $tsv;
  return tsv2hash(@opt);
}

#....................................................................
sub csv2array {
  my(@opt) = @_;
  my $tsv = File::Temp->new();
  run("csvtk csv2tab \Q$opt[0]\E > $tsv");
  $opt[0] = $tsv;
  return tsv2array(@opt);
}

#....................................................................
sub tsv2array {
  my($fname, $sep, $skip_regex) = @_;
  $sep ||= "\t";
  msg("Opening tabular file: $fname");
  open my $TABLE, '<', $fname;
  my $mtx;
  my $ncol = 0;
  #my $count=0;
  while (<$TABLE>) {
    next if $skip_regex and m/$skip_regex/;
    chomp;
    my @row = split m/$sep/;
    push @$mtx, [ @row ];
    $ncol = max($ncol, 0+@row);
  }
  close $TABLE;
  my $count = scalar(@$mtx);
  msg("Loaded $count rows x $ncol cols");
  return $mtx;
}

#....................................................................
sub array2tsv {
  my($fname,$matrix,$header,$sep) = @_;
  $sep ||= "\t";
  my $rows = scalar @$matrix;
  msg("Writing $rows rows to $fname");
  my $out = fhandle($fname);
  unshift @$matrix, $header if $header;
  for my $row (@$matrix) {
    print $out join($sep, @$row),"\n";
  }
  return;
}

#....................................................................
sub hash2tsv {
  my($fname,$hash,$cols,$sep,$null) = @_;
  $sep ||= "\t";
  $null ||= '';
  my $rows = scalar keys %$hash;
  msg("Writing $rows rows to $fname");
  my $ncol = @$cols or err("No columns provided");
  msg("Columns:", @$cols);
  my $out = fhandle($fname);
  print $out join($sep, @$cols),"\n";
  for my $id (keys %$hash) {
    print $out join($sep,
      map { $hash->{$id}{$_} // $null } @$cols
    ),
    "\n";
  }
  return;
}

#....................................................................
sub today { 
  my($sep) = @_;
  $sep ||= '';
  my $t = localtime;
  my $d = $t->ymd;
  $d =~ s/-/$sep/g;
  return $d;
}

#....................................................................
sub fhandle {
  my($fname) = @_;
  open my $fh, '>', $fname 
    or err("Can't open '$fname' for writing");
  return $fh;
}

#....................................................................
sub tmpdir {
  my($change) = @_;
  # CLEANUP->1 should be default for this
  my $dir = File::Temp->newdir();
  $dir or err("Can't create temp folder: $!");
  msg("Created tempdir: $dir");
  if ($change) {
    msg("Changing to tempdir: $dir");
    $old_dir = getcwd();
    chdir $dir;
  }
  return $dir;
}

#....................................................................
# dump \*FH to a temp filename eg. \*ARGV
sub fh2file {
  my($fh, $fname) = @_;

#  $fname //= File::Temp->new();
#  copy($fh, $fname) or err($!);
#  return $fname;

  my $out = $fname ? fhandle($fname)
                   : File::Temp->new();
  my $lines=0;
  msg("Dumping input to $out ...");
  while (<$fh>) {
    print $out $_;
    $lines++;
  }
  msg("Copied $lines lines to $out");
  return $out;
}

#....................................................................
sub csv {
  # replace comma with semicolon
  my @row = @_;
  @row = map { $_ // '' } @row;
  @row = map { s/,/;/g; $_ } @row;
  return join(",", @row)."\n";
}

#....................................................................
sub csv_OLD {
  # replace comma with semicolon
  #my @row = map { s/,/;/g; $_ } @_; 
  #my @row = @_;
  my @row;
  for my $item (@_) {
    $item //= '';
    if ($item =~ m/,/) {
      $item =~ s/"/""/g;
      $item = "\"$item\"";
    }
    push @row, $item;
  }
  join(",", @row)."\n";
}

#....................................................................
sub run { 
  msg("Running: @_"); 
  system(@_)==0 or err("Could not run: @_"); 
}

#....................................................................
sub need {
  exe($_) for @_;
}

#....................................................................
sub exe {
  my($tool) = @_;
  my $path = which($tool);
  $path ? msg("Found '$tool' - $path") 
        : err("Can't find '$tool' - please install it");
  return $path;
}

#....................................................................
sub cpus { 
  my($c) = qx(nproc --all);
  chomp $c;
  return $c || 1;
}

#....................................................................
sub fixcvg { 
  my($n,$text) = @_;
  $text ||= "invalid coverage percent";
  $n =~ m/^[\d.]+$/ or err("$text");
  if ($n > 0 and $n < 1) {
    wrn("Converting coverage $n into percentage");
    $n *= 100;
  }
  if ($n < 0 or $n > 100) {
    err("Coverage $n not in [0,100]");
  }
  return sprintf "%.2f", $n;
}

#....................................................................
sub fixint { 
  my($n,$text) = @_;
  $text ||= "invalid";
  $n =~ m/^\d+$/ or err("$text - need an integer");
  return $n;
}

#....................................................................
sub fixfile { 
  my($f, $text) = @_;
  $text ||= "bummer";
  $f or err("$text - empty filename");
  -r $f or err("$text - can't read file '$f'");
  return abs_path($f);
}

#....................................................................
sub fixdir { 
  my($dir, $text) = @_;
  $text ||= "bummer";
  if (! $dir) {
    err("$text - null folder name");
  }
  elsif (-f $dir) {
    err("'$dir' is a regular file, not a fdlder");
  }
  elsif (-d $dir) {
    wrn("Folder '$dir' already exists.");
  }
  else {
    msg("Making folder '$dir'");
    make_path($dir);
  }
  return abs_path($dir);
}

#....................................................................
sub fixdate { 
  my($date, $allow_blank) = @_;
  $date ||= '';
  #my $dd = join('|', 1..31);
  #my $mm = join('|', 1..12);
  my $dd = '\d{1,2}';
  my $mm = '\d{1,2}';
  my $yy = '\d{2}|\d{4}';
  my($d,$m,$y);
  if ($date =~ m/\b($dd)\D($mm)\D($yy)\b/) {
    ($d,$m,$y) = ($1,$2,$3);
  }
  elsif ($date =~ m/\b($yy)\D($mm)\D($dd)\b/) {
    ($d,$m,$y) = ($3,$2,$1);
  }
  else {
    $allow_blank 
      ? return ''
      : err("Can't parse date: $date");
  }
  $y += 2000 if $y < 2000;
  return sprintf("%04d%02d%02d", $y, $m, $d);
}

#....................................................................
sub get_opts {
  my($u) = @_;
  my $exe = $u->{name} || basename($0);
  my %opt;
  my $optdef = 'hV';
  getopts($optdef, \%opt);
  return \%opt;
}

sub show_usage {
  my($u, $ec) = @_;
  $ec ||= 0;
  my $fh = $ec ? \*STDERR : \*STDOUT;
  select $fh;
  my $exe = $u->{name} || basename($0);
  if ($u->{desc}) {
    print $fh "SYNOPSIS\n  ",$u->{desc},"\n";
  }
  print "USAGE\n  $exe [options] ", ($u->{parm} || ''),"\n";
  print "OPTIONS\n";
  print "  -h       SHow this help\n";
  if ($u->{vers}) {
    print "  -V       Print version and exit\n"
  }
  for my $opt ($u->{opts}->@*) {
#    my($switch,$type,$desc) = split m'~', $opt;
    printf "  -%1s %-7s %s\n",
      split m'~', $opt;
  }
  print "END\n";
  exit($ec)
}

#....................................................................
sub slurp { 
  my($fname, $chomp) = @_;
  $fname or err("Filename is empty string");
  -r $fname or err("can't read file to slurp: $fname");
  open my $FILE, '<', $fname or err("Can't open file to slurp: $fname");
  my @lines = <$FILE>;
  close $FILE;
  chomp @lines if $chomp;
  return wantarray ? @lines : \@lines;
}

#....................................................................
sub slurp_fh { 
  my($fh, $chomp) = @_;
  my @lines = <$fh>;
  chomp @lines if $chomp;
  return wantarray ? @lines : \@lines;
}

#....................................................................
sub spew { 
  my $fname = shift;
  open my $FILE, '>', $fname or err("Can't spew to file: $fname");
  print $FILE @_;
  close $FILE;
}

#....................................................................
sub spew_nl { 
  my($fname,@lines) = @_;
  @lines = grep { defined } @lines;
  if (@lines) {
    spew($fname, map { "$_\n" } @lines);
  }
  else {
    wrn("No lines to write to '$fname'");
  }
}

#....................................................................
sub fasta2array {
  exe("seqkit");
  my $seq;
  for my $fname (@_) {
    msg("Loading FASTA: $fname");
    open my $TAB, '-|', "seqkit fx2tab --only-id '$fname'";
    while (<$TAB>) {
      chomp;
      push @$seq, [ split m/\t/ ];
      # replace empty seq with single N
      $seq->[-1] //= 'N';
    }
    close $TAB;
  }
  msg("Loaded", scalar(@$seq), "sequences.");
  return $seq;
}

#....................................................................
sub afa2array {
  my $seq = fasta2array(@_);
  my $N = @$seq;
  $N or err("No sequences in: @_");
  my $L = length($seq->[0][1]);
  my %seen;
  for my $s (@$seq) {
    my $len = length($s->[1]);
    $len==$L or err($s->[0], "has length $len not $L");
    !$seen{$s->[0]}++ or err("Duplicate ID '".$s->[0]."' in @_");
  }
  msg("Have $N aligned seqs with $L sites.");
  return wantarray ? ($seq, $N, $L) : $seq;
  #return ($seq,$N,$L);
}

#....................................................................
sub fasta2hash {
  my $array = fasta2array(@_);
  msg("Converting FASTA array to hash");
  my $hash = { map { ($_->[0] => $_->[1]) } @$array };
  msg("Unique sequence IDs:", scalar(keys %$hash) );
  return $hash;
}

#....................................................................
sub afa2hash {
  my($array,$N,$L) = afa2array(@_);
  msg("Converting AFA array to hash");
  my $hash = { map { ($_->[0] => $_->[1]) } @$array };
  my $HN = scalar keys %$hash;
  $HN==$N or wrn("Lost", $N-$HN, "sequences in hash conversion");
  return wantarray ? ($hash,$HN,$L) : $hash;
}

#....................................................................
sub array2fasta { 
  my($fh, $seqs) = @_;
  for my $s (@$seqs) {
    print $fh ">", $s->[0], "\n", 
              $s->[1], "\n";
  }
}
#....................................................................
sub hash2fasta { 
  my($fh, $seqs) = @_;
  for my $id (keys %$seqs) {
    print $fh ">", $id, "\n", 
              $seqs->{$id}, "\n";
  }
}

#....................................................................
sub kovid_dir { 
  my $dir = $ENV{'KOVID_HOME'};
  return $dir if $dir && -d $dir;
  # assume commands are in $KOVID_HOME/kovid/bin
  $dir = abs_path(dirname($0).'/../../');
  return $dir if $dir && -d $dir;
  err("Could not determine \$KOVID_HOME folder");
}
sub kovid_ref_id {
  return 'MN908947.3';
}
sub kovid_ref_len {
  return 29903;
}
sub kovid_ref_fasta {
  my $fn = kovid_dir().'/refs/'.kovid_ref_id().'.fna';
  -r $fn or err("Reference doesn't exsit: $fn");
  return $fn;
}
sub kovid_db_file {
  my $fn = kovid_dir().'/DB/DB.json';
  -r $fn or err("DB file doesn't exsit: $fn");
  return $fn;
}
sub kovid_csv_file {
  my $fn = kovid_dir().'/DB/DB.csv';
  -r $fn or err("CSV file doesn't exsit: $fn");
  return $fn;
}
sub kovid_seq_file {
  my $fn = kovid_dir().'/DB/DB.ffn';
  -r $fn or err("SEQ file doesn't exsit: $fn");
  return $fn;
}

#....................................................................
sub kovid_db { 
  my $fn = kovid_db_file();
  msg("Loading database: $fn");
  my $db = decode_json path($fn)->slurp;
  my $N = scalar keys %$db;
  $N or err("Database has no data");
  msg("Database contains $N entries");
  return $db;
}
#....................................................................
sub kovid_db_OLD { 
  my(%filter) = @_;

  my $dir = kovid_dir();

  msg("LOading & merging JSON databases....");
  my $mdu = decode_json path(fixfile("$dir/DB/WGS.json"))->slurp;
  my $qc  = decode_json path(fixfile("$dir/DB/QC.json"))->slurp;
  my $epi = decode_json path(fixfile("$dir/DB/EPI.json"))->slurp;
  my $sra = decode_json path(fixfile("$dir/DB/ACC.json"))->slurp;
  
#  my $db =  merge( merge($mdu, $qc), $sra); 
  my $db = $mdu;
  delete $db->{''}; #VIC01
  
  for my $id (keys %$db) {
    my $nid = $db->{$id}{'Notification_ID'} or next;
    if (exists $epi->{$nid}) {
      $db->{$id} = merge($db->{$id}, $epi->{$nid});
    }
    # delete boring data
    my @boring = grep m/_ONT|_upload|^Tas/, keys $db->{$id}->%*; 
    delete @{$db->{$id}}{@boring};
  }
  
  if (scalar keys %filter) {
    msg("Filtering database:", keys %filter);
    ID: for my $id (keys %$db) {
      for my $k (keys %filter) {
        if (!defined $db->{$id} or $db->{$id}{$k} eq $filter{$k}) {
          delete $db->{$id};
          next ID;
        }
      }
    }
  }
  
  return $db;
}

#....................................................................
sub revcom {
  my($dna) = @_;
  # ATUGCYRSWKMBDHVN atugcyrswkmbdhvn
  # TAACGRYSWMKVHDBN taacgryswmkvhdbn
  $dna =~ tr/ATUGCYRSWKMBDHVNatugcyrswkmbdhvn/TAACGRYSWMKVHDBNtaacgryswmkvhdbn/;
  $dna = reverse($dna);
  return $dna;
}

#....................................................................
sub translate {
  my($dna,$phase) = @_;
  $phase //= 0;
  $dna =~ s/[^ACGT]/N/g;
  $dna = uc substr $dna, $phase;
  my $aa='';
  while ($dna =~ m/(...)/g) {
    $aa .= $CODON{$1} // 'X';;
    #msg("codon=$1 aa=$aa");
  }
  return $aa;
}

#....................................................................
sub subseq {
  my($dna,$begin,$end,$rc) = @_;
  $end ||= $begin;
  $end >= $begin or err("subseq begin($begin) > end($end)");
  my $len = $end-$begin+1;
  $dna = substr($dna, $begin-1, $len);
  $dna = revcom($dna) if $rc;
  return $dna;
}

#....................................................................
sub base_tally {
  # expecting arrayref -> [ ID, SEQ ]
  my($seq) = @_;
  # check it is good
  my $N = scalar @$seq;
  $N > 0 or err("No seqs in alignment");
  my $L = length($seq->[0][1]);
  $L > 0 or err("Alignment has no columns");
  # count up
  msg("Tallying $N seqs x $L sites...");
  my $tally = [];
  for my $s (@$seq) {
    my $count = {};
    map { $count->{$_}++ } split m//, $s->[1];
    push @$tally, $count;
  }
  return $tally;
}

#....................................................................
sub site_tally {
  # expecting arrayref -> [ ID, SEQ ]
  my($seq) = @_;
  # check it is good
  my $N = scalar @$seq;
  $N > 0 or err("No seqs in alignment");
  my $L = length($seq->[0][1]);
  $L > 0 or err("Alignment has no columns");
  # count up
  msg("Tallying $L sites from $N seqs...");
  my $tally = [];
  for my $pos (0 .. $L-1) {
    my $count = {};
    for my $s (@$seq) {
      my $nt = substr($s->[1],$pos,1);
      $count->{ $nt // '-' }++;
    }
    push @$tally, $count;
  }
  return $tally;
}


#....................................................................
sub hash_bin {
  my($string, $len) = @_;
  $len ||= 2;
  my $hex = md5_hex($string); 
  return substr($hex,0,$len);
}

#....................................................................
1;

__DATA__



Base	Name	Bases Represented	Complementary Base
A	Adenine	A	T
T	Thymidine	T	A
U	Uridine(RNA only)	U	A
G	Guanidine	G	C
C	Cytidine	C	G
Y	pYrimidine	C T	R
R	puRine	A G	Y
S	Strong(3Hbonds)	G C	S*
W	Weak(2Hbonds)	A T	W*
K	Keto	T/U G	M
M	aMino	A C	K
B	not A	C G T	V
D	not C	A G T	H
H	not G	A C T	D
V	not T/U	A C G	B
N	Unknown	A C G T	N

ATUGCYRSWKMBDHVN
TAACGRYSWMKVHDBN
