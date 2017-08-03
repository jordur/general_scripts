#!/usr/bin/perl w-

##INTEGRACION EN UNA UNICA TABLA DE TODOS LOS VALORES DE CONTAJE DE LECTURAS DE DIFERENTEES MUESTRAS. INPUT UN FICHERO CON LA RUTA DE LOS DIFERENTES FICHEROS DE LAS MUESTRAS.EL PROGRAMA ESTA DISEÑADO PARA COJER EL OUTPUT DE MIRDEEP2 (MIRAR BF411 PARA VER UN EJEMPLO

if(scalar(@ARGV) <1)
{
	die;
}

$input=$ARGV[0];

%programName=();
%score=();
%maturemiRNAcount=();
%novel=();
%coordinate=();
%identificate=();
%chr=();
%coodinate=();
%sequence=();


@files=();

$count=0;
$entry=0;

$headStringScore="";
$headStringCount="";
$lineScore="";
$lineCount="";

open(FILE,'<',$input);
while($line=<FILE>)
{
        chomp($line);
	$files[$count]=$line;
	$count++;
}

for($i=0;$i<scalar(@files);$i++)
{
	open(FILE,'<',$files[$i]);
	while($line=<FILE>)
	{
		chomp($line);

		@explode=split("\t",$line);
		if($line=~/^chr/)
		{
			$entry=0;

			if($explode[9] eq "-")
			{
				$explode[16]=~s/\.\./_/g;


				@explode2=split(":",$explode[16]);
				@explode3=split("_",$explode2[1]);

				foreach $keys (keys %coordinate)
				{
					if($chr{$keys} eq $explode2[0])
					{
						if($explode3[0] >= $coordinate{$keys}[0] && $explode3[1] <= $coordinate{$keys}[1])
						{
							$entry=1;
							
							$score{$keys}{$files[$i]}=$explode[1];
							$maturemiRNAcount{$keys}{$files[$i]}=$explode[5];
						}
						elsif($explode3[0] <= $coordinate{$keys}[0] && $explode3[0] <= $coordinate{$keys}[1] && $explode3[1] > $coordinate{$keys}[0] && $explode3[1] <= $coordinate{$keys}[1])
						{
							$entry=1;

							$score{$keys}{$files[$i]}=$explode[1];
							$maturemiRNAcount{$keys}{$files[$i]}=$explode[5];
						}
						elsif($explode3[0] >= $coordinate{$keys}[0] && $explode3[0] <= $coordinate{$keys}[1] && $explode3[1] > $coordinate{$keys}[0] && $explode3[1] >= $coordinate{$keys}[1])
						{
							$entry=1;

							$score{$keys}{$files[$i]}=$explode[1];
							$maturemiRNAcount{$keys}{$files[$i]}=$explode[5];
						}
						elsif($explode3[0] >= $coordinate{$keys}[0] && $explode3[0] <= $coordinate{$keys}[1] && $explode3[1] >= $coordinate{$keys}[0] && $explode3[1] <= $coordinate{$keys}[1])
						{
							$entry=1;

							$score{$keys}{$files[$i]}=$explode[1];
							$maturemiRNAcount{$keys}{$files[$i]}=$explode[5];
						}
					}
				}
				if($entry==0)
				{
					$programName{$explode[16]}{$files[$i]}=$explode[16];
					$identificate{$explode[16]}=$explode[16];

					$score{$explode[16]}{$files[$i]}=$explode[1];
					$maturemiRNAcount{$explode[16]}{$files[$i]}=$explode[5];

					$chr{$explode[16]}=$explode2[0];
					$coordinate{$explode[16]}[0]=$explode3[0];
					$coordinate{$explode[16]}[1]=$explode3[1];

					$sequenceString="";

					@explodeS=split("",$explode[13]);
					for($j=0;$j<8;$j++)
					{
						$sequenceString=$sequenceString."".$explodeS[$j];
					}
					$sequence{$explode[16]}=$sequenceString;
				}
			}
			else
			{
				$programName{$explode[9]}{$files[$i]}=$explode[9];
				$identificate{$explode[9]}=$explode[9];
				$score{$explode[9]}{$files[$i]}=$explode[1];
				$maturemiRNAcount{$explode[9]}{$files[$i]}=$explode[5];

				$sequenceString="";
                                @explodeS=split("",$explode[13]);

                                for($j=0;$j<8;$j++)
                                {
                                	$sequenceString=$sequenceString."".$explodeS[$j];
                                }
                                $sequence{$explode[9]}=$sequenceString;
			}
		}
	}
}
for($i=0;$i<scalar(@files);$i++)
{
	@explode=split("/",$files[$i]);
	$scalar=scalar(@explode)-2;
	if($i==0)
	{
		$headStringScore="Score $explode[$scalar]";
		$headStringCount="Count $explode[$scalar]";
	}
	else
	{
		$headStringScore=$headStringScore."\tScore ".$explode[$scalar];
		$headStringCount=$headStringCount."\tCount ".$explode[$scalar];
	}
}

print "Name miRNA\t$headStringScore\t$headStringCount\n";

foreach $keys (keys %identificate)
{
	$lineScore="";
	$lineCount="";

	for($i=0;$i<scalar(@files);$i++)
	{
		if($score{$keys}{$files[$i]})
		{
			if($i==0)
			{
				$lineScore=$score{$keys}{$files[$i]};
			}
			else
			{
				$lineScore=$lineScore."\t".$score{$keys}{$files[$i]};
			}
		}
		else
		{
			if($i==0)
			{
				$lineScore="0";
			}
			else
			{
				$lineScore=$lineScore."\t0";
			}
		}
		if($maturemiRNAcount{$keys}{$files[$i]})
		{
			if($i==0)
			{
				$lineCount=$maturemiRNAcount{$keys}{$files[$i]};
			}
			else
			{
				$lineCount=$lineCount."\t".$maturemiRNAcount{$keys}{$files[$i]};
			}
		}
		else
		{
			if($i==0)
			{
				$lineCount="0";
			}
			else
			{
				$lineCount=$lineCount."\t0"
			}
		}
	}
	print "$identificate{$keys}\t$sequence{$keys}\t$lineScore\t$lineCount\n";
}
