#!/usr/bin/perl

if($ARGV[0] eq "" || $ARGV[1] eq "" || $ARGV[2] eq "" || $ARGV[3] eq "")
{
	die("\n\n\t\tI:ORFs.fasta (En la cabecera se debe indicar el contigs al que pertenece la secuencia. Es necesario que las secuencias del mismo contig esten consecutivas EJM:\n\t\">contig00095\"\n\t\tII:Directory output\t\t\n\t\tIII: Resultado del blastx\n\t\tIV: ORFs_predict\n\n");
}

$input=$ARGV[0];

%entry=();
$string="";
%contig=();
%blastx=();

$system=`mkdir $ARGV[1]`;

open(FILE,'<',$input);
while($line=<FILE>)
{
        chomp($line);

	if($line=~/^>/)
	{
		$line2=$line;
		$line2=~s/>//g;
		if(!$entry{$line2})
		{
			$contig{$line2}=1;
			$string=$ARGV[1]."/Artemis_".$line2.".fasta";
			$entry{$line2}=1;

			open(FILEOUT,'>',$string);

			print FILEOUT "$line\n";
		}
		else
		{
			$string=$ARGV[1]."/Artemis_".$line2.".fasta";
			print FILEOUT "$line\n";
		}
	}
	else
	{
		print FILEOUT "$line\n";
	}
} 
open(FILE2,'<',$ARGV[2]);
while($line=<FILE2>)
{
        chomp($line);
	@explode=split("\t",$line);
	if(!$blastx{$explode[0]})
	{
		$blastx{$explode[0]}=$explode[2];
	}
}

%hash_f=();

open(FILE3,'<',$ARGV[3]);
while($line=<FILE3>)
{
	chmod ($line);

	$start=0;
	$end=0;

        @explode=split("\t",$line);
        @explode2=split("_",$explode[0]);

	if($explode[1]>$explode[2])
	{
		$start=$explode[2];
		$end=$explode[1];
	}
	else
	{
		$start=$explode[1];
		$end=$explode[2];
	}

	if($explode[3]<0)
	{
		$strand="-";
	}
	else
	{
		$strand="+";
	}
	$frame=abs($explode[3]);

	if(!$hash_f{$explode2[0]})
	{
		$hash_f{$explode2[0]}="$explode[0]\tartemis\tCDS\t$start\t$end\t$frame\t$strand\t.\tID=$blastx{$explode[0]}\n";	
	}
	else
	{
		$hash_f{$explode2[0]}=$hash_f{$explode2[0]}."$explode[0]\tartemis\tCDS\t$start\t$end\t$frame\t$strand\t.\tID=$blastx{$explode[0]}\n";

	}
}
foreach $keys (keys %hash_f)
{
	$string=$ARGV[1]."/Artemis_".$keys."_annotation.gff";
	open(OUT,'>',$string);
	print OUT "$hash_f{$keys}";
}


