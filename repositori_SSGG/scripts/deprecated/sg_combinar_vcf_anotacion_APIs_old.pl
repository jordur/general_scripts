#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE:  \n COMO SE USA:\n";
    print "Input: 1) VCF 2) fichero con la anotacion \n";
    print "\n";

    exit(1);
}

if(scalar(@ARGV) < 2)
{
    usage();
}

#X       119602863	C/A   119602863       ENSG00000005893 ENST00000434600 LAMP2   lysosomal-associated  membrane protein 2 [Source:HGNC Symbol;Acc:6501]       INTRONIC        -       -       -       -  2.9000
#X       119602863       .       C       A       75.50   PASS    AC=1;AF=0.50;AN=2;DP=15;Dels=0.00;HRun=2;HaplotypeScore=2.5157;MQ=33.38;MQ0=0;QD=5.03;SB=-0.01;sumGLbyD=7.03        GT:AD:DP:GQ:PL  0/1:10,4:10:99:105,0,189

open(VCF,"<",@ARGV[0]);
open(ANOT,"<",@ARGV[1]);

while($lineas_vcf=<VCF>)
{
	chomp ($lineas_vcf);
	push (@vcf,$lineas_vcf);

}
$j=0;
#17      41256413        AAAA/AAA        41256413        ENSG00000012048 ENST00000309486 BRCA1   breast cancer 1, early onset [Source:HGNC Symbol;Acc:1100]      INTRONIC        -       -    --       -0.3200 ENST00000309486.4:c.-587-138    -587-135delAAAAinsTTT
#13      32890573        G/A     32890573        ENSG00000139618 ENST00000380152 BRCA2   breast cancer 2, early onset [Source:HGNC Symbol;Acc:1101]      5PRIME UTR-       -       -       -       -1.2800 ENST00000380152.3:c.-25C>A
#CACNB2  10      18828663        G       T       SNV     0       heterozygous    rs4747352       ENSG00000165995 ENST00000282343 calcium channel, voltage-dependent, beta 2 subunit [Source:HGNC Symbol;Acc:1402]      3PRIME_UTR      -       -       -       -2.2800 128     99
print "Gene_name\tChr\tPosition\tReference_genotype\tSample_genotype\tVariation_type\tVariation_length\tZygosity\tVariant_ID\tEnsembl_Gene_ID\tEnsembl_Transcript_ID\tGene_description\t\Consequence_on_transcripts\tcDNA_position\taa_position\taa_change\tGERP_conservarion_score\tHGVS_ID\tCondel_Score\tCoverage\tGenotype_quality\tMapping_quality\tInterpro_ID\tInterpro_Description\n";
while ($lineas_anot=<ANOT>)
{
	chomp ($lineas_anot);
	@split_anot=split ("\t",$lineas_anot);
	$split_anot[0]=~s/chr//;
	#@split_linea=split ("_",$split_anot[1]);
	$es1indel="SNV";
	if (&es_indel($split_anot[2]) eq "indel")
	{
		#print "INDEL--------------\n";
		$split_anot[1]=$split_anot[1]-1;	
		$es1indel="indel";
		
	}
	
#	print $lineas_anot,"\n";
	for ($i=$j;$i<=$#vcf;)
	{
		@split_vcf=split ("\t",$vcf[$i]);
		$split_vcf[0]=~s/chr//;
		$split_vcf[7]=~s/;DS;/;/;
		@split_col8 = split (";",$split_vcf[7]);
		$split_col8[1]=~s/AF=//;
		$split_col8[3]=~s/DP=//;
		$split_col8[7]=~s/MQ=//;
		#$split_col8[7]=~s/MQ=//;		
		@split_col10 = split (":",$split_vcf[9]);
#	print "-------VALORES-----------COORDENADA ANNOT $split_anot[1] --- COORDENADA VCF $split_vcf[1]\n";	
		if (($split_vcf[0] eq $split_anot[0]) && ($split_anot[1] == $split_vcf[1]))
		{
#			print "ENTRA\n";
			$split_anot[13]=~s/\./,/;
			$split_col10[3]=~s/\./,/;
			$split_col8[7]=~s/\./,/;	
			if ($es1indel eq "SNV")
			{
				print $split_anot[6],"\t",$split_anot[0],"\t",$split_anot[1],"\t",$split_vcf[3],"\t", $split_vcf[4],"\t",&longitud_variacion ($split_anot[2]),"\t",&homoheterozigoto($split_col8[1]),"\t",$split_anot[12],"\t",$split_anot[4],"\t",$split_anot[5],"\t",$split_anot[7],"\t",$split_anot[8],"\t",$split_anot[9],"\t",$split_anot[10],"\t",$split_anot[11],"\t",$split_anot[13],"\t",$split_anot[14],"\t",$split_anot[15],"\t",$split_col8[3],"\t", $split_col10[3],"\t",$split_col8[7],"\t",$split_anot[16],"\t",$split_anot[17],"\n"; 
			}
			else
			{
				print $split_anot[6],"\t",$split_anot[0],"\t",$split_anot[1],"\t",$split_vcf[3],"\t", $split_vcf[4],"\t",&longitud_variacion ($split_anot[2]),"\t",&homoheterozigoto($split_col8[1]),"\t",$split_anot[12],"\t",$split_anot[4],"\t",$split_anot[5],"\t",$split_anot[7],"\t",$split_anot[8],"\t",$split_anot[9],"\t",$split_anot[10],"\t",$split_anot[11],"\t",$split_anot[13],"\t",$split_anot[14],"\t",$split_anot[15],"\t",$split_col8[3],"\t", $split_col10[2],"\t",$split_col8[7],"\t",$split_anot[16],"\t",$split_anot[17],"\n";
			}
			$j=$i;
			last;
		}
		elsif ($split_anot[1] > $split_vcf[1])
		{	
			$i++;
			$j=$i;
		}
		else
		{
			last;
		}
	}		


}

sub es_indel
{
	$_=shift;
        @split_variation = split ("/",$_);
        chomp($split_variation[0]);
        chomp ($split_variation[1]);
        if (($split_variation[0] eq "-") || ($split_variation[1] eq "-") || (length($split_variation[0]) != length ($split_variation[1])))
	#if ((length($split_variation[0]) != length ($split_variation[1])))
	{		
	#	print "es indel\n";
		return "indel";
	}
	else
	{
		return "SNV";
	}
	
}

sub longitud_variacion 
{
	$_=shift;
	@split_variation = split ("/",$_);
	chomp($split_variation[0]);
	chomp ($split_variation[1]);
	$lalong = abs(length($split_variation[0])-length($split_variation[1]));
	#print "-----LONGITUD-------$lalong\n";
	if ($split_variation[0] eq "-")
	#if ($lalong > 0 )
	{
		return "indel\t",(length ($split_variation[1]));
	}	
	elsif ($split_variation[1] eq "-")
	{
		return "indel\t",(length ($split_variation[0]));
	}
	elsif (length ($split_variation[1]) eq length ($split_variation[0]))
	{
		return "SNV\t-";
	}
	elsif (($split_variation[0] ne "-") && ($split_variation[1] ne "-") && ($lalong > 0) )
        {
                return "indel\t",$lalong;
        }
}

sub homoheterozigoto 
{
	$_=shift;
	if ($_ eq "0.50")
	{
		return "heterozygous";
	}
	else
	{
		return "homozygous";
	}


}
