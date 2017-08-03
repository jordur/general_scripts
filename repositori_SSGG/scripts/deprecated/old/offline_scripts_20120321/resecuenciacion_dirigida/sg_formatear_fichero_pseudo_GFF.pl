#!/usr/bin/perl -w

open (GFF,"<",$ARGV[0]);

while ($lines=<GFF>)
{
	chomp ($lines);
	@spltlinea=split ("\t",$lines);
	#print $spltlinea[2],"\n";
	$char = substr ($spltlinea[14],0,1);
	@splitAlleles=split ("/",$spltlinea[9]);
	
	if ($spltlinea[2] eq "deletion")
	{
		
	#	&hemihomo($spltlinea[13]),"\n";
		print "$spltlinea[0]\t",$spltlinea[3]-1,"\t.\t$spltlinea[14]\t$char\t255\t.\tAC=;AF=",&hemihomo($spltlinea[13]),";AN=;DP=",$spltlinea[12],";Dels=;HRun=;HaplotypeScore=;MQ=;MQ0=;SF=1\tGT:AD:DP:GQ:PL\t0:0:0:0:0\n";
	}
	elsif (($spltlinea[2] eq "insertion_site") && ($splitAlleles[0] eq ""))
	{
		print "$spltlinea[0]\t$spltlinea[3]\t.\t$char\t",$char.$splitAlleles[1],"\t255\t.\tAC=;AF=",&hemihomo($spltlinea[13]),";AN=;DP=",$spltlinea[12],";Dels=;HRun=;HaplotypeScore=;MQ=;MQ0=;SF=1\tGT:AD:DP:GQ:PL\t0:0:0:0:0\n";  
	}

}

sub hemihomo
{
	$valor = shift ;
	#print $valor,"\n";
	if ($valor eq "HOMOZYGOUS")
	{
		return "1.00";
	}
	else
	{
		return "0.50";
	}

}
