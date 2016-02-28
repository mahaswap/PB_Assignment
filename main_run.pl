#!/usr/bin/perl 
# PB_assignment
#
# Created by Swapnil MAHAJAN on 9th Aug. 2011.
use File::Basename;

@args=split(/\,/,$ARGV[0]);
$pwd=`pwd`; chomp $pwd; # Current directory path

$rand="";
$rand=$ARGV[1];
chomp $rand;

# To store all arguments
$flag_pdb=0;$flag_out=0;
foreach $a(@args)
{
  @tmp_arg=split(/\=/,$a);
	if($tmp_arg[0] eq "pdb"){$pdb=$tmp_arg[1]; $flag_pdb=1;}     # pdb file
	elsif($tmp_arg[0] eq "out"){$out=$tmp_arg[1]; $flag_out=1;}  # output dir where all output files will be stored
}

# To check if all arguments have been entered or not.
if($flag_pdb==0)
{
	print "Please enter PDB file\n";
	print "USAGE\: perl main_run.pl pdb=/Users/swapnil/WORK/pentaDB/1DWJ.pdb,out=/Users/swapnil/WORK/pentaDB/OUT\n";
	exit;
}
elsif($flag_out==0)
{
	print "Please enter output directory path\n";
	print "USAGE\: perl main_run.pl pdb=/Users/swapnil/WORK/pentaDB/1DWJ.pdb,out=/Users/swapnil/WORK/pentaDB/OUT\n";
	exit;
}
else
{
  if(!-e $pdb){print "Can not find $pdb\n";exit;}
}

( $name, $path, $suffix) = fileparse($pdb, qr/\.[^.]*$/);
`mkdir -p $out`;    # Creates output directory to store PB sequences
open (PDB,$pdb)||die("$0\:$pdb not found\:$!\n");

open (OUT,">$out/$name\_$rand\.pdb")||die("$0\:tmp\.pdb\:$!");

foreach $line(<PDB>)
{
  if($line=~/^END /)
  {
    print OUT "TER\n";
    close OUT;
    $op=`perl read_pdb.pl $out/$name\_$rand.pdb $out`;   # Script to assign torsions and PB. This script can create *.pbseq, *.aaseq, *.penta, (*.tor commented) files.
    if(-e "$out/$name\_$rand.pbseq")
    {
    `cat $out/$name\_$rand.pbseq >>$out/$rand\_$name.pbseq`;
    `cat $out/$name\_$rand.aaseq >>$out/$rand\_$name.aaseq`;
    `cat $out/$name\_$rand.tor >>$out/$rand\_$name.tor`;
    `rm -f $out/$name\_*`;
    }
  }
  if($line=~/^TER / || $line=~/^ENDMDL/)
  {
    print OUT "TER\n";
    close OUT;
    $op=`perl read_pdb.pl $out/$name\_$rand.pdb $out`;   # Script to assign torsions and PB. This script can create *.pbseq, *.aaseq, *.penta, (*.tor commented) files.
    if(-e "$out/$name\_$rand.pbseq")
    {
    `cat $out/$name\_$rand.pbseq >>$out/$rand\_$name.pbseq`;
    `cat $out/$name\_$rand.aaseq >>$out/$rand\_$name.aaseq`;
    `cat $out/$name\_$rand.tor >>$out/$rand\_$name.tor`;
    `rm -f $out/$name\_*`;
    }
    open (OUT,">$out/$name\_$rand\.pdb")||die("$0\:tmp\.pdb\:$!");
  }
  if($line=~/^ATOM / || $line=~/^HETATM /)
  {
    print OUT $line;
  }
}
#`rm -f $name\_*`;
