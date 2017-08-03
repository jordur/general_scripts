#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Extrae informacion de la Biobase en un fichero de SNV e Indels (Separado por cromosoma y ordenado por su posicion \n COMO SE USA:\n";
    print "Input: 1) fichero de SNVs e Indels annotados (por cromosoma y ordenado)\n";
    print "       2) Biobase por Cromosoma\n";
    
    print " Output: Fichero de SNVs e Indels annotados (por cromosoma y ordenado) con informacion de la Biobase \n";
    print "\n";
    

    exit(1);
}

if(scalar(@ARGV) < 2)
{
    usage();
}

my %hash=();

my @explode1=();
my @explode2=();
my @annotationExplode=();
my @finalColumn=();
my @f=();
my @SNVspecific=();
my @SNVspecific2=();
my @changeIndels=();
my @changeIndelsTo=();

my $line="";
my $line2="";
my $chr="";

my $position=0;
my $positionA=0;
my $positionB=0;
my $start=0;
my $end=0;
my $scalar=0;
my $count=0;
my $SNV="NO";
my $entry=0;
my $once=0;
my $fistLine=0;

my $identification="";
my $change="";
my $changeTo="";
my $change2="";
my $changeTo2="";


open(FILE1,"<","$ARGV[1]");
@f=<FILE1>;

open(FILE,"<",$ARGV[0]);
while($line=<FILE>)
{
	chomp($line);

	$change="";
	$change2="";
	$changeTo="";
	$changeTo2="";
	
	$once=0;

	if($fistLine==0)
	{
		#print "$line\tDisease\tHGVS\tPMID\n";
		$fistLine=1;
		
	}
	else
	{
		$SNVs="NO";
		$entry=0;

		@annotationExplode=split("\t",$line);

	
		if($annotationExplode[5] eq "SNV") #IS SNVs
        	{
      	 	        $position=$annotationExplode[2];
       		        $SNVs="YES";
        	}
		else
		{
#			$positionA=$annotationExplode[2]+1; # VCF TO BIOBASE
			$positionA=$annotationExplode[2];


			#$positionB=($annotationExplode[2]+$annotationExplode[6])+1;
			$positionB=($annotationExplode[2]+$annotationExplode[6]);
			$SNVs="NO";

		}
		$change=$annotationExplode[3];
		$changeTo=$annotationExplode[4];

		$change2=&complementario($change);
		$changeTo2=&complementario($changeTo);

		$line=~s/SNV\t-/SNV\t0/;

		for($k=$count;$k<scalar(@f);$k++)
		{
			chomp($f[$k]);
	
			$biostart=0;
			$bioend=0;

			@explode1=split("\t",$f[$k]);
			@explode2=split(":",$explode1[4]);

			if($explode2[1] !~/-/)
			{
				$biostart=$explode2[1];
				$bioend=$explode2[1];
			}
			else
			{
				@explode3=split("-",$explode2[1]);
				$biostart=$explode3[0];
				$bioend=$explode3[1];
			}
			$strand=$explode2[2];

			if($SNVs eq "YES")
			{
				if($position >= $biostart && $position <= $bioend)	
				{
					if($explode1[8] eq "-" || $explode1[8] eq "null") ##ESA COLUMNA ESTA VACIA
					{
						$entry=1;
						$count=$k;
						
						if($once==0)
						{
							$line=~s/\t$annotationExplode[8]\t/\t$annotationExplode[8]\t$explode1[1]\t/;
							$once=1;
						}
						else
						{
							 $line=~s/$annotationExplode[9]/$annotationExplode[9],$explode1[1]/;
						}

						$line=~s/SNV\t0/SNV\t-/;

						if($explode1[11]!~/[A-Z]/)##PUBMED IDENTIFICATION
						{
							print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";
						}
						else
						{
							print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
						}
					}
					else
					{
						$explode1[8]=~s/\[/\t/;
						$explode1[8]=~s/\]/\t/;
						$explode1[8]=~s/\//-/;
						@SNVspecific=split("\t",$explode1[8]);

						if($SNVspecific[1]=~/-/)
						{
							@SNVspecific2=split("-",$SNVspecific[1]);

							if($strand eq "+")
							{
								if($change  eq $SNVspecific2[0] && $changeTo eq $SNVspecific2[1])
								{
									$entry=1;
									$count=$k;
	
									if($once==0)
									{
										$line=~s/\t$annotationExplode[8]\t/\t$annotationExplode[8]\t$explode1[1]\t/;
										$once=1;
									}
									else
									{
										$line=~s/$annotationExplode[9]/$annotationExplode[9],$explode1[1]/;
									}

									$line=~s/SNV\t0/SNV\t-/;

									if($explode1[11]!~/A-Z]/) ## PUBMED IDENTIFICATION
									{
										print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";
									}
									else
									{
										print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
									}
								}
								else ###BECAUSE BIOBASE DATABASE AREN'T STANDARD HGVS, IF SOME BASE IS EQUAL TO REFERENCE IS OK
								{
									if($change  eq $SNVspecific2[1] && $changeTo eq $SNVspecific2[0])
									{
										$entry=1;
										$count=$k;

										if($once==0)
										{
											$line=~s/\t$annotationExplode[8]\t/\t$annotationExplode[8]\t$explode1[1]\t/;
											$once=1;
										}
										else
										{
											$line=~s/$annotationExplode[9]/$annotationExplode[9],$explode1[1]/;
										}

										$line=~s/SNV\t0/SNV\t-/;

										if($explode1[11]!~/A-Z]/) ## PUBMED IDENTIFICATION
										{
											print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";
										}
										else
										{
											print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
										}
									}
								}
							}
							else
							{
								if($change2 eq $SNVspecific2[0] && $changeTo2 eq $SNVspecific2[1])#COMPLEMENTARY STRAND
								{
									$entry=1;
									$count=$k;

									if($once==0)
									{
										$line=~s/\t$annotationExplode[8]\t/\t$annotationExplode[8]\t$explode1[1]\t/;
										$once=1;
									}
									else
									{
										$line=~s/$annotationExplode[9]/$annotationExplode[9],$explode1[1]/;
									}

									$line=~s/SNV\t0/SNV\t-/;

									if($explode[11]!~/A-Z]/)
									{
										print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";
									}
									else
									{
										print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
									}
								}
								else ###BECAUSE BIOBASE DATABASE AREN'T STANDARD HGVS, IF SOME BASE IS EQUAL TO REFERENCE IS OK
								{
									if($change2 eq $SNVspecific2[1] && $changeTo2 eq $SNVspecific2[0])
									{
										$entry=1;
										$count=$k;
										
										if($once==0)
										{
											$line=~s/\t$annotationExplode[8]\t/\t$annotationExplode[8]\t$explode1[1]\t/;
											$once=1;
										}
										else
										{
											$line=~s/$annotationExplode[9]/$annotationExplode[9],$explode1[1]/;
										}
										
										$line=~s/SNV\t0/SNV\t-/;

										if($explode[11]!~/A-Z]/)
										{
	                                                                               print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";
										}
										else
										{
											 print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
										}
									}
								}
							}
						}
					}
				}
			}
			else
			{
				if($explode1[scalar(@explode1)-1]=~/Indels/)
				{
					$start_indel=0;
					$end_indel=0;
					$entryIndels=0;
					if($explode2[1]=~/-/)
					{
						@explode3=split("-",$explode2[1]);
						$start_indel=$explode3[0];
						$end_indel=$explode3[1];
					}
					else
					{
						$start_indel=$explode2[1]; ## OFTEN DELECTIONS
					}
					
					if($end_indel!=0)
					{
						if($positionA ==  $start_indel || $positionA == $start_indel -1 )
						{
							$entryIndels=1;
						}
					}
					else
					{
						if($start_indel == $positionA || $start_indel -1 == $positionA)  ## SOMETIMES THE BIOBASE DELECTION WITH ONLY START COORDINATE HAVE REFERENCE TO NEXT NUCLEOTIDE
						{
							$entryIndels=1;
						}
					}
				}			

				if($explode1[18] eq "Indels" && $entryIndels==1)
				{
					if($explode1[8] eq "-" || $explode1[8] eq "null") ##COLUMNA VACIA
					{
						$entry=1;
						$count=$k;

						if($once==0)
                                                {
	                                                $line=~s/\t$annotationExplode[8]\t/\t$annotationExplode[8]\t$explode1[1]\t/;
                                                        $once=1;
                                                }
                                                else
                                                {
        	                                        $line=~s/$annotationExplode[9]/$annotationExplode[9],$explode1[1]/;
                                                }

						if($explode[11]!~/A-Z]/)
						{
							print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";
						}
						else
						{
							 print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
						}
					}
					else
					{
						$explode1[8]=~s/\[/\t/;
						$explode1[8]=~s/\]/\t/;

						@indelSpecific_tmp=split("\t",$explode1[8]);
						@indelSpecific=split("/",$indelSpecific_tmp[1]);
						if($indelSpecific[0] eq "-")
						{
							$len_0=0;
						}
						else
						{
							$len_0=length($indelSpecific[0]);
						}
	
						if($indelSpecific[1] eq "-")
						{
							$len_1=0;
						}
						else
						{
							$len_1=length($indelSpecific[1]);
						}
						$indelSpecific[0]=~s/-//;
						$indelSpecific[1]=~s/-//;
						if(($len_0-$len_1)!=0)
						{
							$indel_type2=$indelSpecific[0];
							$indel_type2=~s/$indelSpecific[1]//;

							if(($len_0-$len_1)>0)
							{

								$variation=$change;
								$variation=~s/$changeTo//;


								if($strand eq "-")
								{
									$variation=&complementario($variation);
								}

								if($indel_type2 eq $variation)
								{
									$entry=1;
									$count=$k;
									
									if($once==0)
									{
										$line=~s/\t$annotationExplode[8]\t/\t$annotationExplode[8]\t$explode1[1]\t/;
										$once=1;
									}
									else
									{
										$line=~s/$annotationExplode[9]/$annotationExplode[9],$explode1[1]/;
									}
									if($explode[11]!~/A-Z]/)
									{
										print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";
									}
									else
									{
										print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
									}
								}
							}
							if($len_0-$len_1 < 0)
							{
								$indel_type2==$indelSpecific[1];
								$indel_type2=~s/$indelSpecific[0]//;

								$variation=$changeTo;
								$variation=~s/$change//;
								if($strand eq "-")
								{
									$variation=&complementario($variation);
								}
								if($indel_type2 eq $variation)
								{
									$entry=1;
									$count=$k;

									if($once==0)
									{
										$line=~s/\t$annotationExplode[8]\t/\t$annotationExplode[8]\t$explode1[1]\t/;
										$once=1;
									}
									else
									{
										$line=~s/$annotationExplode[9]/$annotationExplode[9],$explode1[1]/;
									}
									if($explode[11]!~/A-Z]/)
									{
										print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";
									}
									else
									{
										print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
									}
								}
							}
						}
						else
						{
							$entry=1;
							$count=$k;
							if($once==0)
							{
								$line=~s/\t$annotationExplode[8]\t/\t$annotationExplode[8]\t$explode1[1]\t/;
								$once=1;
							}
							else
							{
								$line=~s/$annotationExplode[9]/$annotationExplode[9],$explode1[1]/;
							}

							if($explode[11]!~/A-Z]/)
							{
								print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";
							}
							else
							{
								print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
							}
						}
					}
				}	
			}
		}
		if($entry==0)
		{
			$line=~s/\t$annotationExplode[8]\t/\t$annotationExplode[8]\t-\t/;
			print "$line\t-\t-\t-\n";
		}
	}
}
sub complementario
{
	$y="";
	$x=shift;
	@split=split("",$x);

	for(my $ii=0;$ii<scalar(@split);$ii++)
	{
		if($split[$ii] eq "A")
		{
			$y="T".$y;
		}	
		if($split[$ii] eq "T")
		{
			$y="A".$y;
		}
		if($split[$ii] eq "C")
		{
			$y="G".$y;
		}
		if($split[$ii] eq "G")
		{
			$y="C".$y;
		}
	}
	return($y);
}
