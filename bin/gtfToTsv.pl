#Scripted by Martin Taylor for LCE paper 2020
# Convert GTF file to TSV for simple R import
use strict;
while(<STDIN>){
  chomp;
  my @sp = split /\t/;
  my @sq = split /;\s*/, $sp[8];
  my %keys = ('gene_id' => "NA",
    'gene_name' => 'NA',
    'transcript_id' => 'NA',
    'gene_biotype' => 'NA',
    );
  for my $r (0 .. $#sq){
    my @ss = split /\s+/, $sq[$r];
    if (defined $keys{$ss[0]}){
      $ss[1] =~ s/"//g;
      $keys{$ss[0]} = $ss[1];
    }
  }
  if($keys{'gene_name'} eq "NA"){
    $keys{'gene_name'} = $keys{'gene_id'};
  }
  print join "\t", ($sp[0],$sp[3],$sp[4],$keys{'gene_id'},1,$sp[6],$keys{'gene_name'},$keys{'transcript_id'},$keys{'gene_biotype'},$sp[2]);
  print "\n";
}
