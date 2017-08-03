#!/usr/bin/perl -w
use POSIX;


if($ARGV[0] eq "" || $ARGV[1] eq "")
{
        die("\n\n\t\tI:Input data,type F3 qual(one o more)\n\n");
}
else
{
	if($ARGV[1] ==1)
	{
		preprocess_file($ARGV[0]);
	}
	if($ARGV[1] == 2)
	{
		random_file($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);
	}
	if($ARGV[1] == 3)
	{
		merge_file();
	}
}
sub merge_file
{
	open(OUT_IND_POS,'>',"indt_per_positions_final.stats");
	open(OUT_IND_MEAN,'>'."indt_number_final.stats");
	open(OUT_QUAL_MEAN,'>',"qual_mean_final.stats");
	open(OUT_QUAL_POS,'>',"qual_per_positions_final.stats");
	open(OUT_SEQ_STATS,'>',"sequence_positions_final.stats");

	my $reads_total=0;
	my $entry=0;
	my $mean_ind_pos=0;

	my @indetermination_per_position=();
	my @indetermination_mean=();
	my @quality_mean=();
	my @quality_positions=();
	my $MAX=();
	my $MIN=();
	my $MEAN=();
	my $scalar_quality_mean=0;
	my $scalar_indetermination_mean=0;
	my $scalar_quality_position=0;
	my $scalar_indetermination_per_positions=0;

	my $systemD=`ls -latr *_sequence_positions.stats`;
	@lines=split("\n",$systemD);
	for($i=0;$i<scalar(@lines);$i++)
	{
		@explode=split(" ",$lines[$i]);
		print "\t$explode[8] *** process\t";

		open(FILE,'<',$explode[8]);
		while($line=<FILE>)
		{
			chomp($line);
			@explode_line=split("\t",$line);
			if(!$MAX[$explode_line[0]])
			{
				$MAX[$explode_line[0]]=$explode_line[2];
			}
			else
			{
				if($explode_line[2] >= $MAX[$explode_line[0]])
				{
					$MAX[$explode_line[0]]=$explode_line[2];
				}
			}

			if(!$MIN[$explode_line[0]])
			{
				$MIN[$explode_line[0]]=$explode_line[3];
			}
			else
			{
				if($explode_line[3] < $MIN[$explode_line[0]])
				{
					$MIN[$explode_line[0]]=$explode_line[3];
				}
			}

			if(!$MEAN[$explode_line[0]])
			{
				$MEAN[$explode_line[0]]=$explode_line[1];
			}
			else
			{
				$MEAN[$explode_line[0]]=$MEAN[$explode_line[0]]+$explode_line[1];
			}
		}
		print "OK!\n";
	}	
	for($t=1;$t<scalar(@MEAN);$t++)
	{
		$MEAN[$t]=$MEAN[$t]/(scalar(@lines));
		print OUT_SEQ_STATS  "$t\t$MEAN[$t]\t$MAX[$t]\t$MIN[$t]\n";
	}


	my $system=`ls -latr *indt_per_positions.stats`;
	@lines=split("\n",$system);
	for($i=0;$i<scalar(@lines);$i++)
	{
		$entry=0;

		@explode=split(" ",$lines[$i]);

		print "\t$explode[8] *** process\t";

		open(FILE,'<',$explode[8]);
		while($line=<FILE>)
		{
			chomp($line);
			@explode_line=split("\t",$line);
			if($entry==0)
			{
				$reads_total=$reads_total+$explode_line[2];
				$entry=1;
			}

			if(!$indetermination_per_position[$explode_line[0]])
			{
				$indetermination_per_position[$explode_line[0]]=$explode_line[1];
			}
			else
			{
				$indetermination_per_position[$explode_line[0]]=$indetermination_per_position[$explode_line[0]]+$explode_line[1];
			}
			$scalar_indetermination_per_positions=$scalar_indetermination_per_positions+$explode_line[1];
		}
		print "OK!\n";
	}
	print "\n";

	my $systemA=`ls -latr *_indt_number.stats`;	
	@lines=split("\n",$systemA);
	for($i=0;$i<scalar(@lines);$i++)
	{
		@explode=split(" ",$lines[$i]);

		print "\t$explode[8] *** process\t";

		open(FILE,'<',$explode[8]);
		while($line=<FILE>)
		{
			chomp($line);
			@explode_line=split("\t",$line);
			if(!$indetermination_mean[$explode_line[0]])
			{
				$indetermination_mean[$explode_line[0]]=$explode_line[1];
			}
			else
			{
				$indetermination_mean[$explode_line[0]]=$indetermination_mean[$explode_line[0]]+$explode_line[1];
			}
			$scalar_indetermination_mean=$scalar_indetermination_mean+$explode_line[1];
		}
		print "OK!\n";
	}
	print "\n";

	my $systemB=`ls -latr *_quality_mean.stats`;
	@lines=split("\n",$systemB);
	for($i=0;$i<scalar(@lines);$i++)
	{
		@explode=split(" ",$lines[$i]);

		print "\t$explode[8] *** process\t";
		
		open(FILE,'<',$explode[8]);
		while($line=<FILE>)
		{
			chomp($line);
			@explode_line=split("\t",$line);
			if(!$quality_mean[$explode_line[0]])
			{
				$quality_mean[$explode_line[0]]=$explode_line[1];
			}
			else
			{
				$quality_mean[$explode_line[0]]=$quality_mean[$explode_line[0]]+$explode_line[1];
			}
			$scalar_quality_mean=$scalar_quality_mean+$explode_line[1];
		}
		print "OK!\n";
	}
	print "\n";

	my $systemC=`ls -latr *_quality_positions.stats`;
	@lines=split("\n",$systemC);
	for($i=0;$i<scalar(@lines);$i++)
	{
		@explode=split(" ",$lines[$i]);
		
		print "\t$explode[8] *** process\t";
		open(FILE,'<',$explode[8]);
		while($line=<FILE>)
		{
			chomp($line);
			@explode_line=split("\t",$line);
			if(!$quality_positions[$explode_line[0]])
			{
				$quality_positions[$explode_line[0]]=$explode_line[1];
			}
			else
			{
				$quality_positions[$explode_line[0]]=$quality_positions[$explode_line[0]]+$explode_line[1];
			}
			$scalar_quality_position=$scalar_quality_position+$explode_line[1];
		}
		print "OK!\n";
	}

	
	for($j=1;$j<scalar(@indetermination_per_position);$j++)
	{
		if($scalar_indetermination_per_positions!=0)
		{
			$indetermination_per_position[$j]=($indetermination_per_position[$j]/$scalar_indetermination_per_positions)*100;
		}
		if(!$indetermination_per_position[$j])
		{
			print OUT_IND_POS  "$j\t0\t$scalar_indetermination_per_positions\n";
		}
		else
		{
			print OUT_IND_POS "$j\t$indetermination_per_position[$j]\t$scalar_indetermination_per_positions\n";
		}
	}

	for($k=0;$k<(scalar(@indetermination_mean)+1);$k++)
	{
		if($indetermination_mean[$k])
		{
			$indetermination_mean[$k]=($indetermination_mean[$k]/$scalar_indetermination_mean)*100;
			print OUT_IND_MEAN "$k\t$indetermination_mean[$k]\n";
		}
		else
		{
			print OUT_IND_MEAN "$k\t0\n";
		}
	}
	for($l=0;$l<scalar(@quality_mean);$l++)
	{
		$quality_mean[$l]=($quality_mean[$l]/$scalar_quality_mean)*100;
		if(!$quality_mean[$l])
		{
			print OUT_QUAL_MEAN "$l\t0\t$scalar_quality_mean\n";
		}
		else
		{
			print OUT_QUAL_MEAN  "$l\t$quality_mean[$l]\t$scalar_quality_mean\n";
		}
	}
	for($m=0;$m<scalar(@quality_positions);$m++)
	{
		if($quality_positions[$m])
		{
			if($scalar_quality_position > 0)
			{
				$quality_positions[$m]=($quality_positions[$m]/$scalar_quality_position)*100;
			}

			if(!$quality_positions[$m])
			{
				print OUT_QUAL_POS "$m\t0\t$scalar_quality_position\n";
			}
			else
			{
				print OUT_QUAL_POS  "$m\t$quality_positions[$m]\t$scalar_quality_position\n";
			}
		}
	}
}


sub random_file
{
	my $input=$_[0];
	my $random=$_[2];
	my $scalar_file=$_[3];
	my $identification_indetermination=$_[4];
	my $indetermination_threshold=0;
#
	if($identification_indetermination!~/y/)
	{
		$indetermination_threshold=-1;
	}
	my $output_indetermination_number=$_[0]."_indt_number.stats";
	my $output_indetermination_positions=$_[0]."_indt_per_positions.stats";
	my $output_quality_mean=$_[0]."_quality_mean.stats";
	my $output_quality_positions=$_[0]."_quality_positions.stats";
	my $output_random=$_[0]."_random_lines.ghost";
	my $output_sequence=$_[0]."_sequence_positions.stats";
	open(OUT_IND_NUM,'>',$output_indetermination_number);
	open(OUT_IND_POS,'>',$output_indetermination_positions);
	open(OUT_QUAL_MEAN,'>',$output_quality_mean);
	open(OUT_QUAL_POS,'>',$output_quality_positions);
	open(OUT_SEQ_STATS,'>',$output_sequence);

	my @mean_indetermination=();
	my @mean_quality=();
	my @indetermination_per_positions=();
	my @quality_per_position=();
	my @size_quality_per_positions=();
	my @quality_per_positions_bar=();
	my @MAX=();
	my @MIN=();
	my @rand=();

	my $rand_entry=0;
	my $rand_count=0;
	my $sub_entry=0;
	my $count_indetermination=0;
	my $max_indetermination=0;
	my $count=0;
	my $length_fasta=0;
	my $min_length_fasta=10000;
	my $mean_quality_variable=0;
	my $max_quality=0;

	open(OUT_RANDOM,'>',$output_random);

        while($rand_entry==0)
        {
                $rand_number=int(rand($scalar_file));
                if(!$rand[$rand_number])
                {
                        $rand[$rand_number]=1;
                        $rand_count++;
                }
                if($rand_count==$random)
                {
                        $rand_entry=1;
                }
        }
 	open(FILE,'<',$input);
	while($line=<FILE>)
	{
		chomp($line);
		if($line!~/^>/ && $line!~/^#/)
		{
			if($rand[$count])
			{
				$sub_entry=0;
				$count_indetermination=0;
				$mean_quality_variable=0;

				$line=~s/  / /g;

				print OUT_RANDOM "$line\n";

				@explode=split(" ",$line);

				$length_fasta=scalar(@explode)+1;
				if($length_fasta < $min_length_fasta)
				{
					$min_length_fasta = $length_fasta;
				}

				for($i=0;$i<scalar(@explode);$i++)
				{
					if($explode[$i] == $indetermination_threshold) ## ESTOY VALE SOLO PARA RESECUENCIACION CON SOLID, PARA JUNIOR ADAPTARLO DESPUES CON VARIABLES
					{
						$count_indetermination=$count_indetermination+1;
						if(!$indetermination_per_positions[$i])
						{
							$indetermination_per_positions[$i]=1;
						}
						else
						{
							$indetermination_per_positions[$i]=$indetermination_per_positions[$i]+1;
						}
					}
					else
					{
						$sub_entry=1;
						$mean_quality_variable=$mean_quality_variable+$explode[$i];
					}
					if(!$quality_per_position[$i])
					{
						$quality_per_position[$i]=$explode[$i];
						$size_quality_per_positions[$i]=1;
					}
					else
					{
						$quality_per_position[$i]=$quality_per_position[$i]+$explode[$i];
						$size_quality_per_positions[$i]++;
					}
					if($explode[$i]!="-1")
					{	
						if(!$quality_per_positions_bar[$explode[$i]])
						{
							$quality_per_positions_bar[$explode[$i]]=1;
						}
						else
						{
							$quality_per_positions_bar[$explode[$i]]=$quality_per_positions_bar[$explode[$i]]+1;
						}
					}
					
					if(!$MAX[$i])
					{
						$MAX[$i] = $explode[$i];
					}
					else
					{
						if($explode[$i] > $MAX[$i])
						{
							$MAX[$i] = $explode[$i];
						}
					}
					if(!$MIN[$i])
					{
						$MIN[$i] = $explode[$i];
					}
					else
					{
						if($explode[$i] < $MIN[$i])
						{
							$MIN[$i] = $explode[$i];
						}
					}
				}
				if($count_indetermination > $max_indetermination)
				{
					$max_indetermination=$count_indetermination;
				}
				if(!$mean_indetermination[$count_indetermination])
				{
					$mean_indetermination[$count_indetermination]=1;
				}
				else
				{
					$mean_indetermination[$count_indetermination]=$mean_indetermination[$count_indetermination]+1;
				}
				#$mean_indetermination[$count_indetermination]=$mean_indetermination[$count_indetermination]/(scalar(@explode)-1);
				$mean_quality_variable=$mean_quality_variable/(scalar(@explode)-1);

				if($mean_quality_variable > $max_quality)
				{
					$max_quality=$mean_quality_variable;
				}


				if(!$mean_quality[$mean_quality_variable])
				{
					$mean_quality[$mean_quality_variable]=1;
				}
				else
				{
					$mean_quality[$mean_quality_variable]=$mean_quality[$mean_quality_variable]+1;
				}
			}
			$count++;
		}
	}
	


	for($j=0;$j<($max_indetermination+1);$j++)
	{
		if(!$mean_indetermination[$j])
		{
			print OUT_IND_NUM "$j\t0\t$random\n";
		}
		else
		{
			print OUT_IND_NUM "$j\t$mean_indetermination[$j]\t$random\n";
		}
	}
	for($k=0;$k<($length_fasta)-1;$k++)
	{
		$kk=$k+1;
		if(!$indetermination_per_positions[$k])
		{
			print OUT_IND_POS "$kk\t0\t$random\n";
		}
		else
		{
			print OUT_IND_POS "$kk\t$indetermination_per_positions[$k]\t$random\n";
		}
	}
	for($l=0;$l<($max_quality +1);$l++)
	{
		if(!$mean_quality[$l])
		{
			print OUT_QUAL_MEAN "$l\t0\n";
		}
		else
		{
			print OUT_QUAL_MEAN "$l\t$mean_quality[$l]\n";
		}
	}
	for($s=0;$s<(@quality_per_positions_bar);$s++)
	{
		$ss=$s+1;
		if(!$quality_per_positions_bar[$s])
		{
			print OUT_QUAL_POS "$ss\t0\n";
		}
		else
		{
#			$quality_per_position[$s]=$quality_per_position[$s]/$size_quality_per_positions[$s];
			print OUT_QUAL_POS "$ss\t$quality_per_positions_bar[$s]\n";
		}
	}
	for($r=0;$r<scalar($length_fasta)-1;$r++)
	{
		$rr=$r+1;
		if(!$quality_per_position[$r])
		{
			print OUT_SEQ_STATS "$rr\t0\t$MAX[$r]\t$MIN[$r]\n";
		}
		else
		{
			$quality_per_position[$r] = $quality_per_position[$r]/$size_quality_per_positions[$r];
			print OUT_SEQ_STATS "$rr\t$quality_per_position[$r]\t$MAX[$r]\t$MIN[$r]\n";
		}
	}
}

sub preprocess_file
{
	my $input =$_[0];
	$string ="";

	open(FILE,'<',$input);
        while($line=<FILE>)
        {
		chomp($line);

		if($line!~/^#/)
		{
			if($line=~/^>/)
			{
				if($string ne "")
                        	{
					$string=~s/  / /g;
					print  "$string\n";
					$string="";
				}
				print "$line\n";
			}
			else
			{
				if($string eq "")
				{
					$string=$line;
				}
				else
				{
					$string=$string." ".$line;
                	         }
			}
		}
	}
	if($string ne "")
	{
		$string=~s/ / /g;
		print  "$string\n";
	}
}
