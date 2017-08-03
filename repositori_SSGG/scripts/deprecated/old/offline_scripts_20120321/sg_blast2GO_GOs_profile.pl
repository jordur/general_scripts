#!/usr/bin/perl -w

if(!$ARGV[0])
{
	die("\n\n\t***\tI:Input File (GO annotation of blast2Go,first column is sequences name, then  is GO_BP,GO_CC nad GO_MF\t***\n\n");
}

$input =$ARGV[0];

@explode=();
@explodeBP=();
@explodeCC=();
@explodeMF=();

%BP=();
%BPcount=();
%MF=();
%MFcount=();
%CC=();
%CCcount=();

open($BP_OUT,'>',"GO_Biological_Process.txt");
open($MF_OUT,'>',"GO_Molecular_Function.txt");
open($CC_OUT,'>',"GO_Cellular_Component.txt");

open(FILE,'<',$input);
while($line=<FILE>)
{
	chomp($line);
	@explode=split("\t",$line);
	@explodeBP=split(";",$explode[1]);
	for($i=0;$i<scalar(@explodeBP);$i++)
	{
		if($explodeBP[$i]=~/\w/)
		{
			$explodeBP[$i]=~s/^\s//;
			$explodeBP[$i]=~s/\s$//g;
			if(!$BP{$explodeBP[$i]})
			{
				$BP{$explodeBP[$i]}=$explode[0];
				$BPcount{$explodeBP[$i]}=1;
			}
			else
			{
				$BP{$explodeBP[$i]}=$BP{$explodeBP[$i]}.",".$explode[0];
				$BPcount{$explodeBP[$i]}=$BPcount{$explodeBP[$i]}+1;
			}
		}
	}	
	@explodeCC=split(";",$explode[2]);
	for($i=0;$i<scalar(@explodeCC);$i++)
	{
		if($explodeCC[$i]=~/\w/)
		{
			$explodeCC[$i]=~s/^\s//;
			$explodeCC[$i]=~s/\s$//g;

			if(!$CC{$explodeCC[$i]})
			{
				$CC{$explodeCC[$i]}=$explode[0];
				$CCcount{$explodeCC[$i]}=1;
			}
			else
			{
				$CC{$explodeCC[$i]}= $CC{$explodeCC[$i]}.",".$explode[0];
				$CCcount{$explodeCC[$i]}=$CCcount{$explodeCC[$i]}+1;
			}
		}
	}
	@explodeMF=split(";",$explode[3]);
	for($i=0;$i<scalar(@explodeMF);$i++)
	{
		if($explodeMF[$i]=~/\w/)
		{
			$explodeMF[$i]=~s/^\s//;
			$explodeMF[$i]=~s/\s$//g;
		
			if(!$MF{$explodeMF[$i]})
			{
				$MF{$explodeMF[$i]}=$explode[0];
				$MFcount{$explodeMF[$i]}=1;
			}
			else
			{
				$MF{$explodeMF[$i]}=$MF{$explodeMF[$i]}.",".$explode[0];
				$MFcount{$explodeMF[$i]}=$MFcount{$explodeMF[$i]}+1;
			}
		}
	}
}

print $BP_OUT "GO Biological Process\tGO Count\tORFs\n";
print $CC_OUT "GO Cellular Component\tGO Count\tORFs\n";
print $MF_OUT "GO Molecular Function\tGO Count\tORFs\n";

foreach $keys (keys %BP)
{
	if($keys!~/^GO /)
	{
		print $BP_OUT "$keys\t$BPcount{$keys}\t$BP{$keys}\n";
	}
}
foreach $keys (keys %CC)
{
	if($keys!~/^GO /)
	{
		print $CC_OUT "$keys\t$CCcount{$keys}\t$CC{$keys}\n";
	}
}
foreach $keys (keys %MF)
{
	if($keys!~/^GO /)
	{
		print $MF_OUT "$keys\t$MFcount{$keys}\t$MF{$keys}\n";
	}
}



