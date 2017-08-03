#!/usr/bin/perl -w

if($ARGV[0] eq "" || $ARGV[1] eq "")
{
	die("\n\n\t\tI:Csfasta  input\n\t\tII: Random sampling (ejm 10000000) \n\n");
}
else
{

	$input=$ARGV[0];
	$random=$ARGV[1];

	@sequencesObtainFromRandom=();
	@percentage=();
	@indeter_positions=();
	
	$entry=0;
	$count=0;
	$count_indet=0;
	$percent=0;
	$size=0;
	$number_indetermination=0;

	$system=`grep -c ">" $input`;
	@explode=split(" ",$system);

	while($entry==0)
	{
		$rand=int(rand($explode[0]));
		if(!$sequencesObtainFromRandom[$rand])
		{
			$sequencesObtainFromRandom[$rand]=1;
			$count++;
		}
		if($count==$random)
		{
			$entry=1;
		}
	}

	$count=0;
	$entry=0;

	open(FILE,'<',$input);
	while($line=<FILE>)
	{
	        chomp($line);
		if($line=~/^>/)
		{
			$count++;
			if($sequencesObtainFromRandom[$count])
			{
				if($sequencesObtainFromRandom[$count]==1)
				{
					$entry=1;
				}
			}
			else
			{
				$entry=0;
			}
		}
		else
		{
			if($entry==1)
			{
#print "$line\n";
				$count_indet=0;
				@explode=split("",$line);

				$size=scalar(@explode);

				for($j=0;$j<scalar(@explode);$j++)
				{
					if($explode[$j] eq ".")
					{
						if(!$indeter_positions[$j])
						{
							$indeter_positions[$j]=1;
						}
						else
						{
							$indeter_positions[$j]=$indeter_positions[$j]+1;
						}
						$count_indet++;
						$number_indetermination++;
					}
				}
				if(!$percentage[$count_indet])
				{
					$percentage[$count_indet]=1;
				}
				else
				{
					$percentage[$count_indet]=$percentage[$count_indet]+1;
				}
			}
		}
 	}

	print "** Indetermination by Positions\n";

	for($s=1;$s<$size;$s++)
	{
		if($s==1)
		{
			for($ss=1;$ss<$size;$ss++)
			{
				if($ss==1)
				{
					print "$ss";
				}
				else
				{
					print "\t$ss";
				}
			}
			print "\n";
			if(!$indeter_positions[$s])
			{
				print "0";
			}
			else
			{
				print "$indeter_positions[$s]";
			}
		}
		else
		{
			if(!$indeter_positions[$s])
			{
				print "\t0";
			}
			else
			{
				print "\t$indeter_positions[$s]";
			}
		}
	}

        for($s=1;$s<$size;$s++)
        {
                if($s==1)
                {
                        print "\n";
                        if(!$indeter_positions[$s])
                        {
                                print "0";
                        }
                        else
                        {
				$percent=$indeter_positions[$s]/$number_indetermination;
                                print "$percent";
                        }
                }
                else
                {
                        if(!$indeter_positions[$s])
                        {
                                print "\t0";
                        }
                        else
                        {
				$percent=$indeter_positions[$s]/$number_indetermination;
                                print "\t$percent";
                        }
                }
        }



	print "\n\n**  Number Reads by Number of Indetermination\n";
	print "Indetermination Number\t#Reads\t\%Reads\n";
	for($z=0;$z<scalar(@percentage);$z++)
	{
		if($percentage[$z])
		{
			$percent=$percentage[$z]/$random;
			print "$z\t$percentage[$z]\t$percent\n";
		}
		else
		{
			print "$z\t0\t0\n";
		}
	}

}
