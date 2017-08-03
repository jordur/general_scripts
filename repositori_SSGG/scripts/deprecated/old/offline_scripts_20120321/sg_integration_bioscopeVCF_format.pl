#!/usr/bin/perl 

##SIRVE PARA TRNASFORMAR UN GFF QUE ES LA SALIDA DEL BIOSCOPE SMALLINDELS EN UN VCF EN FORMATO INTEGRADO (CON GATK Y SAMTOOLS) PARA SU POSTERIOR ANOTACION POR EL ENSEMBL

if($ARGV[0] eq "")
{
	die("\n\n\t\t** I: gff Input (bioscope smallIndels outout **\n\n");
}

$input=$ARGV[0];

$chr=0;
$start=0;
$end=0;
$id=0;
$ref=0;
$alt=0;
$qualData="NODATA";
$null=".";
$DP=0;
$AF=0;
$AC="AC=";
$AN="AN=";
$Dels="Dels=";
$HRun="HRun=";
$HaplotypeScore="HaplotypeScore=";
$MQ="MQ=";
$MQ0="MQ0=";


$columnX="GT:AD:DP:GQ:PL";
$columnZ="0:0:0:0:0";

open(FILE,'<',$input);
while($line=<FILE>)
{
	chomp($line);
	$primer_car=substr ($line,0,1);
   if ( $primer_car ne "#")
   {
	@explode=split("\t",$line);

	$chr=0;
	$start=0;
	$end=0;
	$id=0;
	$ref=0;
	$alt=0;
	$qualData="NODATA";
	$null=".";
	$DP=0;
	$AF=0;
	$AC="AC=";
	$AN="AN=";
	$Dels="Dels=";
	$HRun="HRun=";
	$HaplotypeScore="HaplotypeScore=";
	$MQ="MQ=";
	$MQ0="MQ0=";

	$explode[0]=~s/chrchr/chr/g;
	
	$chr=$explode[0];
	$start=$explode[3];
	$end=$explode[4];
	@explode2=split(";",$explode[8]);

	$explode2[4]=~s/allele-call=//g;
	@explode3=split("/",$explode2[4]);

	if($explode3[0]=~/\w/)
	{
		$ref=$explode3[0];
	}
	else
	{
		$ref="";
	}
	if($explode3[1]=~/\w/)
	{
		$alt=$explode3[1];
	}
	else
	{
		$alt="";
	}
	$explode2[10]=~s/no_nonred_reads/DP/g;
	@explode4=split("=",$explode2[12]);
	if($explode4[1] eq "HOMOZYGOUS")
	{
		$AF="AF=1.0";
	}
	else
	{
		$AF="AF=0.50";
	}
	$DP=$explode2[10];

	print "$chr\t$start\t.\t$ref\t$alt\t$qualData\t$null\t$AC;$AF;$AN;$DP;$Dels;$HRun;$HaplotypeScore;$MQ;$MQ0\t$columnX\t$columnZ\n";
   }
}
