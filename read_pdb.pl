#!/usr/bin/perl 
# PB_assignment
#
# Created by Swapnil MAHAJAN on 8th Aug. 2011.

use File::Basename;
use  Math::Trig;


$complete_path=$ARGV[0];
$out_path=$ARGV[1];
( $name, $path, $suffix) = fileparse($complete_path, qr/\.[^.]*$/);

print "############\n$name\n";


open(PDB, "$complete_path")||print("$0\:$complete_path not found: $!");

$count=0; # to count total number of residues in pdb file
$res_no=0;
$prev_res_no=-10000;  # residue numbers
$flag_ter=0; # if chain or model ends then $flag_ter=1
@res_num=(); @resi=();
@N=(); @CA=(); @C=(); @torsions=();
@pb=();
$pb[0]="Z";$pb[1]="Z";

%aa=("ALA"=>"A","ASX"=>"B","CYS"=>"C","ASP"=>"D","GLU"=>"E","PHE"=>"F","GLY"=>"G","HIS"=>"H","ILE"=>"I","LYS"=>"K","LEU"=>"L",
"MET"=>"M","MSE"=>"M","ASN"=>"N","PRO"=>"P","GLN"=>"Q","ARG"=>"R","SER"=>"S","THR"=>"T","VAL"=>"V","UNK"=>"U","TRP"=>"W","TYR"=>"Y","GLX"=>"Z");

foreach $line(<PDB>)
{
  if($line=~/^TER / || $line=~/^ENDMDL/){$flag_ter=1;}
  if (($line=~/^ATOM / || $line=~/^HETATM /) && $flag_ter==0)
  {
    $res_no=substr($line,22,5);
		# Updated by Swapnil MAHAJAN on 10th Aug. 2011. removed "or" condition from "if"; because some HETATM start without "N" at begining
		#if (($res_no ne $prev_res_no || substr($line,13,3)=~/N /) && substr($line,16,1)=~/[A| ]/) # New res starts with 'N' atom. Alternate_position (substr, 16,1) should have same res. no.
    if ($res_no ne $prev_res_no && substr($line,16,1)=~/[A| ]/) # New res starts with 'N' atom. Alternate_position (substr, 16,1) should have same res. no.
		{
			if(exists($aa{substr($line,17,3)})){$resi[$count]=$aa{substr($line,17,3)};}
			else{$resi[$count]="X";}  # Updated by Swapnil MAHAJAN on 10th Aug. 2011. To add "X" if PDB 3 letter code is not in %aa
			
      $count++;
      $tmp=$res_no;
      $tmp=~s/\s//g;
			push(@res_num, $tmp);
      #print "$tmp ";
    }
    if(substr($line,13,3)=~/N / && substr($line,16,1)=~/[A| ]/)
    {
			$N[$count][0]= substr($line,30,8);
			$N[$count][1]= substr($line,38,8);
			$N[$count][2]= substr($line,46,8);
			
			$N[$count][0]=~s/ //g;$N[$count][1]=~s/ //g;$N[$count][2]=~s/ //g;
		}
		elsif(substr($line,13,3)=~/CA / && substr($line,16,1)=~/[A| ]/)
		{
			$CA[$count][0]= substr($line,30,8);
			$CA[$count][1]= substr($line,38,8);
			$CA[$count][2]= substr($line,46,8);
			
			$CA[$count][0]=~s/ //g;$CA[$count][1]=~s/ //g;$CA[$count][2]=~s/ //g;
			$chain=substr($line,21,1);
		}
		elsif(substr($line,13,3)=~/C / && substr($line,16,1)=~/[A| ]/)
		{
			$C[$count][0]= substr($line,30,8);
			$C[$count][1]= substr($line,38,8);
			$C[$count][2]= substr($line,46,8);
			
			$C[$count][0]=~s/ //g;$C[$count][1]=~s/ //g;$C[$count][2]=~s/ //g;
		}
		
  }
  
  $prev_res_no=$res_no;
#  print $line;
}
shift @N; shift @CA; shift @C;
$n_len=@N; $ca_len=@CA; $c_len=@C; $resi_len=@resi;
print "Total resudues: $count\tRESIDUES: $resi_len\nN: $n_len\tCA: $ca_len\tC: $c_len\n";


if($n_len==0 && ca_len==0 && $c_len==0){print "No backbone atoms in the pdb file\n";exit;}
if($n_len==0 && ca_len>0 && $c_len==0){print "Only CA atoms in the pdb file\n";exit;}

open(PHI_PSI, ">$out_path/$name\.tor")||print("$path/$name\.tor can not be written: $!");
#open(PENTA, ">$out_path/$name\.penta")||print("$0\:$out_path/$name\.penta can not be written: $!");
open(AASEQ, ">$out_path/$name\.aaseq")||print("$0\:$out_path/$name\.aaseq can not be written: $!");
open(PBSEQ, ">$out_path/$name\.pbseq")||print("$0\:$out_path/$name\.pbseq can not be written: $!");

print PHI_PSI "Residue No, Chain, Residue, CA-CA Distance, Phi, Psi\n";
$torsions[0][0]="none";   # First phi angle = None
for (my $i=0;$i<($count-1);$i++)  # Calculate torsions
{
	$dist_ca=100;
	@tmp_arr=();
	@tmp_arr=($CA[$i][0],$CA[$i][1],$CA[$i][2],$CA[$i+1][0],$CA[$i+1][1],$CA[$i+1][2]);
	$dist_ca=&dist(\@tmp_arr);   # Calculate CA-CA distance
	$dist_ca=sprintf("%.4f",$dist_ca);
	if($dist_ca>4)# if CA-CA distance >4 then next
	{
		$torsions[$i][1]="none"; $torsions[$i+1][0]="none"; 
		print PHI_PSI "$res_num[$i]\,$chain\,$resi[$i]\,$dist_ca\,$torsions[$i][0]\,$torsions[$i][1]\n";
	}   
	else
	{
		$torsions[$i][1]=&calc_torsions($N[$i][0],$N[$i][1],$N[$i][2],$CA[$i][0],$CA[$i][1],$CA[$i][2],$C[$i][0],$C[$i][1],$C[$i][2],$N[$i+1][0],$N[$i+1][1],$N[$i+1][2]);
		$torsions[$i+1][0]=&calc_torsions($C[$i][0],$C[$i][1],$C[$i][2],$N[$i+1][0],$N[$i+1][1],$N[$i+1][2],$CA[$i+1][0],$CA[$i+1][1],$CA[$i+1][2],$C[$i+1][0],$C[$i+1][1],$C[$i+1][2]);
		$torsions[$i][1]=sprintf("%.4f",$torsions[$i][1]);
		$torsions[$i+1][0]=sprintf("%.4f",$torsions[$i+1][0]);
		print PHI_PSI "$res_num[$i]\,$chain\,$resi[$i]\,$dist_ca\,$torsions[$i][0]\,$torsions[$i][1]\n";
	}
}
$torsions[$count-1][1]="none";  # Last psi angle = None
print PHI_PSI "$res_num[$count-1]\,$chain\,$resi[$count-1]\,0.0000,$torsions[$count-1][0]\,$torsions[$count-1][1]\nTER\n";

for(my $i=0;$i<$count-4;$i++)
{
	my $aaseq="";my $pos="";my @tmp_dihed=();
	$aaseq="$resi[$i]$resi[$i+1]$resi[$i+2]$resi[$i+3]$resi[$i+4]"; # pentapeptide seq
	$pos=$res_num[$i];  # position
	@tmp_dihed=($torsions[$i][1],$torsions[$i+1][0],$torsions[$i+1][1],$torsions[$i+2][0],$torsions[$i+2][1],$torsions[$i+3][0],$torsions[$i+3][1],$torsions[$i+4][0]);  # torsion angles corresponding to pentapepptide
	$flag_none=0;
	foreach my $d(@tmp_dihed){if($d eq "none" || $d eq ""){$flag_none=1;last;}}  # Updated by Swapnil MAHAJAN on 10th Aug. 2011. Added or $d eq "" to avoid null torsions.
	if($flag_none==0){($rmsda[$i+2],$pb[$i+2])=&calc_rmsda(\@tmp_dihed);}
	else{$rmsda[$i+2]="NA";$pb[$i+2]="X";}
	
	#print PENTA "$name\,$pb[$i+2]\,$aaseq\,$pos\,$torsions[$i][1],$torsions[$i+1][0],$torsions[$i+1][1],$torsions[$i+2][0],$torsions[$i+2][1],$torsions[$i+3][0],$torsions[$i+3][1],$torsions[$i+4][0]\,$rmsda[$i+2]\n";
}

$pb[$count-2]="Z";$pb[$count-1]="Z";
$rmsda[$count-2]="NA";$rmsda[$count-1]="NA";

$tp_aa=join("",@resi);
if($tp_aa=~/[A-W|Y-Z]/ig)
{
print AASEQ "\>$name\#$chain\n";
for (my $i=0; $i<$count;$i++)
{
  print AASEQ "$resi[$i]";
}
print AASEQ "\n";

print PBSEQ "\>$name\#$chain\n";
for (my $i=0; $i<$count;$i++)
{
  print PBSEQ "$pb[$i]";
}
print PBSEQ "\n";
}
else
{close PHI_PSI; `rm -f $out_path/$name\.tor`;}
###############################################################################################
# To calculate distance between 2 atoms
sub dist
{
	$pi=3.1415926535;
  my $distance=1000;
	my @arr_dist=();
	@arr_dist=@{$_[0]};
	$distance=sqrt(($arr_dist[3]-$arr_dist[0])**2+($arr_dist[4]-$arr_dist[1])**2+($arr_dist[5]-$arr_dist[2])**2);
	return $distance;
}

# To calculate phi psi angles of 4 atoms 
sub calc_torsions
{
  my (@a1,@a2,@a3,@a4,@b1,@b2,@b3)=();
	my ($num,$deno,$rad,$deg,$len)=0;
	$len = @_;
	if($len==0 || $len<12){return "ERR";}
	@a1=($_[0], $_[1], $_[2]);
	@a2=($_[3], $_[4], $_[5]);
	@a3=($_[6], $_[7], $_[8]);
	@a4=($_[9], $_[10], $_[11]);
	
	@b1 = @{vector(\@a1,\@a2)};
	@b2 = @{vector(\@a2,\@a3)};
	@b3 = @{vector(\@a3,\@a4)};	
	
	$rad= atan2(dotProd(constProd(mod(\@b2),\@b1),crossProd(\@b2,\@b3)),dotProd(crossProd(\@b1,\@b2),crossProd(\@b2,\@b3)));
	$deg=$rad*(180/$pi);
	#print "$rad\n$deg\n";
	
	return $deg;	
}

sub vector
{
	my @p1 = @{$_[0]};
	my @p2 = @{$_[1]};
	
	return [$p2[0]-$p1[0] , $p2[1]-$p1[1] , $p2[2]-$p1[2]];
}

sub mod
{
	my @v = @{$_[0]};
	return sqrt($v[0]**2 + $v[1]**2 + $v[2]**2);
}

sub constProd
{
	my $const = $_[0];
	my @v = @{$_[1]};
	
	return [$const*$v[0] , $const*$v[1] , $const*$v[2]];
}

sub dotProd
{
	my @v1 = @{$_[0]};
	my @v2 = @{$_[1]};
	return $v1[0]*$v2[0] + $v1[1]*$v2[1] + $v1[2]*$v2[2];
}

sub crossProd
{
	my @v1 = @{$_[0]};
	my @v2 = @{$_[1]};
	
	return [$v1[1]*$v2[2]-$v2[1]*$v1[2] , $v2[0]*$v1[2]-$v1[0]*$v2[2] , $v1[0]*$v2[1]-$v2[0]*$v1[1]];
}

# To calculate RMSDa and assign PB using standard torsions for PBs. 
sub calc_rmsda
{
  my %std_dihed=(
	"A"=>[41.14,75.53,13.92,-99.80,131.88,-96.27,122.08,-99.68],
	"B"=>[108.24,-90.12,119.54,-92.21,-18.06,-128.93,147.04,-99.90],
	"C"=>[-11.61,-105.66,94.81,-106.09,133.56,-106.93,135.97,-100.63],
	"D"=>[141.98,-112.79,132.20,-114.79,140.11,-111.05,139.54,-103.16],
	"E"=>[133.25,-112.37,137.64,-108.13,133.00,-87.30,120.54,77.40],
	"F"=>[116.40,-105.53,129.32,-96.68,140.72,-74.19,-26.65,-94.51],
	"G"=>[0.40,-81.83,4.91,	-100.59,85.50,-71.65,130.78,84.98],
	"H"=>[119.14,-102.58,130.83,-67.91,121.55,76.25,-2.95,-90.88],
	"I"=>[130.68,-56.92,119.26,77.85,10.42,-99.43,141.40,-98.01],
	"J"=>[114.32,-121.47,118.14,82.88,-150.05,-83.81,23.35,-85.82],
	"K"=>[117.16,-95.41,140.40,-59.35,-29.23,-72.39,-25.08,-76.16],
	"L"=>[139.20,-55.96,-32.70,-68.51,-26.09,-74.44,-22.60,-71.74],
	"M"=>[-39.62,-64.73,-39.52,-65.54,-38.88,-66.89,-37.76,-70.19],
	"N"=>[-35.34,-65.03,-38.12,-66.34,-29.51,-89.10,-2.91,77.90],
	"O"=>[-45.29,-67.44,-27.72,-87.27,5.13,77.49,30.71,-93.23],
	"P"=>[-27.09,-86.14,0.30,59.85,21.51,-96.30,132.67,-92.91]
	);
	
	my @penta_arr=();
	@dihed=@{$_[0]};
	
	my $min_diff=999999; my $min_pb="";

	foreach my $k(sort(keys(%std_dihed)))
	{
		my $diff=0;
		for(my $n=0;$n<@dihed;$n++)
		{
			my $diff1=0;
			$diff1=&ang_diff($dihed[$n],${$std_dihed{$k}}[$n]);
			$diff+=($diff1**2);
		}
		$diff=sqrt($diff/8);
		if($diff<$min_diff){$min_diff=$diff;$min_pb=$k;}
	}
	
	return($min_diff,$min_pb);
}



# This subroutine calculates the real difference in two angles i.e. (-179)-(+179)=-358 but real difference is 2.
sub ang_diff
{
  $secondAngle=$_[1]; $firstAngle=$_[0];
  my $difference = $secondAngle - $firstAngle;
  while ($difference < -180) {$difference += 360;}
  while ($difference > 180) {$difference -= 360;}
  return $difference;
}
