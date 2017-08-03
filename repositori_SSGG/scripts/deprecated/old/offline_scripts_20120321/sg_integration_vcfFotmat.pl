#!/usr/bin/perl

##SIRVE PARA UNIFICAR EL CONTENIDO DEL VCF PARA REALIZAR LA ANOTACION DE LAS VARIANTES, USADO PARA PANELES DE GENES (CARDIO Y CANCER POR EJEMPLO), 454 O EXOMA

if($ARGV[0] eq "")
{
	print "\n\n\t\t** I:VCF input **\n\n";
}

$input=$ARGV[0];

$AC="AC=";
$AF="AF=";
$AN="AN=";
$DP="DP=";
$Dels="Dels=";
$Hrun="HRun=";
$HaplotypeScore="HaplotypeScore=";
$MQ="MQ=";
$MQ0="MQ0=";
$QD="QD=";
$entryDP=0;
#$columns8="GT:AD:DP:GQ:PL";
#$columns8="GT:DP:GQ:PL";

open(FILE,'<',$input);
while($line=<FILE>)
{
	chomp($line);
	@explode=split("\t",$line);

	if(length($explode[3]) == length($explode[4]))
	{
		$columns8="GT:AD:DP:GQ:PL";
	}
	else
	{
		$columns8="GT:DP:GQ:PL";
	}

	# To avoid problems in the annotation, the "REF" and "ALT" fields will be filled up with "-" in case that they were empty:
	if($explode[3] eq ""){
		$explode[3]="-";
	}
	if($explode[4] eq ""){
		$explode[4]="-";
	}

	print "$explode[0]\t$explode[1]\t$explode[2]\t$explode[3]\t$explode[4]\t$explode[5]\t$explode[6]";

	$AC="AC=";
	$AF="AF=";
	$AN="AN=";
	$DP="DP=";
	$Dels="Dels=";
	$Hrun="HRun=";
	$HaplotypeScore="HaplotypeScore=";
	$MQ="MQ=";
	$MQ0="MQ0=";
	$QD="QD=";
	$entryDP=0;

	@explode2=split(";",$explode[7]);
	for($i=0;$i<scalar(@explode2);$i++)
	{
		if($explode2[$i]=~/^AC/)
		{
			$AC=$explode2[$i];
		}
		if($explode2[$i]=~/^AF/)
		{
			$AF=$explode2[$i];
		}
		if($explode2[$i]=~/^AN/)
		{
			$AN=$explode2[$i];
		}
		if($explode2[$i]=~/^DP=/)
		{
			$DP=$explode2[$i];
			if($DP ne "DP=")
			{
				$entryDP=1;
			}
		}
		if($explode2[$i]=~/^Dels/)
		{
			$Dels=$explode2[$i];
		}	
		if($explode2[$i]=~/^HRun/)
		{
			$Hrun=$explode2[$i];
		}
		if($explode2[$i]=~/^HaplotypeScore/)
		{
			$HaplotypeScore=$explode2[$i];
		}
		if($explode2[$i]=~/^MQ=/)
		{
			$MQ=$explode2[$i];
		}
		if($explode2[$i]=~/^MQ0=/)
		{
			$MQ0=$explode2[$i];
		}
		if($explode2[$i]=~/^QD/)
		{
			$QD=$explode2[$i];
		}
	}
	if($entryDP==0)
	{
		@explode3=split(":",$explode[8]);
		@explode4=split(":",$explode[9]);
		for($i=0;$i<scalar(@explode3);$i++)
		{
			if($explode3[$i] eq "DP")
			{
				$DP=$DP."".$explode4[$i];
			}
		}
	}
	@explode8=split(":",$explode[8]);
	@explode9=split(":",$explode[9]);
	
	$subGT="0";
	$subAD="1,1";
	if($DP ne "0")
	{
		$subDP=$DP;
		$subDP=~s/DP=//g;
	}
	else
	{
		$subDP="0";
	}
	$subGQ="0";
	$subPL="0";

	for($i=0;$i<scalar(@explode8);$i++)
	{
		if($explode8[$i] eq "GT")
		{
			$subGT=$explode9[$i];
		}
		if($explode8[$i] eq "AD")
		{
			$subAD=$explode9[$i];
		}
		if($explode8[$i] eq "DP")
		{
			$subDP=$explode9[$i];
		}
		if($explode8[$i] eq "GQ")
                {
			$subGQ=$explode9[$i];
		}
		if($explode8[$i] eq "PL")
		{
			$subPL=$explode9[$i];
		}
	}
	if(length($explode[3]) == length($explode[4]))
	{
		$string9=$subGT.":".$subAD.":".$subDP.":".$subGQ.":".$subPL;
	}
	else
	{
		$string9=$subGT.":".$subDP.":".$subGQ.":".$subPL;
	}

	print "\t$AC;$AF;$AN;$DP;$Dels;$Hrun;$HaplotypeScore;$MQ;$MQ0\t$columns8\t$string9\n";
}
