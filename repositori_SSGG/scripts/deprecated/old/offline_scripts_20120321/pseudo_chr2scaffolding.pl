#!/usr/bin/perl -w

if($ARGV[0] eq "" || $ARGV[1] eq "" )
{
	die("\n\n\t\tI: Input File,coordinate of scaffolding in pseudo_chr in Order:Scafoldding_Name\tlength_fasta_Scaffolding\tPseudo_Chr_number\n\t\tII:Sam File, Coordinate Sort(Ejm Picard Tools)\n\n");
}

$input=$ARGV[0];
$sam_input=$ARGV[1];

%hash_length=();
%hash_chr=();

open(FILE,'<',$input)|| die("\t\t\nI can't open \"$input\" file\n\n");
while($line=<FILE>)
{
        chomp($line);
	@explode=split("\t",$line); 

	if(!$hash_chr{$explode[2]})
	{
		$hash_chr{$explode[2]}=$explode[0];
	}
	else
	{
		$hash_chr{$explode[2]}=$hash_chr{$explode[2]}.",".$explode[0];
		
	}
	$hash_length{$explode[0]}=$explode[1];
}

open(SAM,'<',$sam_input)|| die("\t\t\nI can't open \"$sam_input\" file\n\n");
while($line=<SAM>)
{
        chomp($line);
	$total=0;
	$coorA="*";
	$chrA=0;
	$coorB=0;
	$chrB="*";
	$entryA=0;
	$entryB=0;
	$entry_equal=0;


        @explode=split("\t",$line);

	if($explode[2] ne "*" && $explode[3]!=0)
	{

		@explode2=split(",",$hash_chr{$explode[2]});


		for($i=0;$i<scalar(@explode2);$i++)
		{
			$j=$i-1;
			if($entryA==0)
			{
				if($j!=-1)
				{
					if($explode[3] >= $total && $explode[3] <= ($hash_length{$explode2[$i]}+$total))
					{
						$coorA=$explode[3]-$total;
						$chrA=$explode2[$i];
						$entryA=1;
					}
				}
				else
				{
					if($explode[3] <= $hash_length{$explode2[$i]})
					{
						$coorA=$explode[3];
						$chrA=$explode2[$i];
						$entryA=1;
					}
				}
			}
			$total=$total+$hash_length{$explode2[$i]};
		}
	}
	if($explode[6] ne "*" && $explode[7] !=0)
	{
		if($explode[6] eq "=")
		{
			$explode[6]=$explode[2];
			$entry_equal=1;
		}
	
	
		if($explode[6] eq "*")
		{
			$line2=$line;
			$coorB=$explode[7];
			$chrB=$explode[6];
		}
		else
		{
			$total=0;
	
			@explode2=split(",",$hash_chr{$explode[6]});
			for($i=0;$i<scalar(@explode2);$i++)
			{
				$j=$i-1;
				if($entryB==0)
				{
					if($j!=-1)
					{
						if($explode[7] >= $total && $explode[7] <= ($hash_length{$explode2[$i]}+$total))
						{
							$coorB=$explode[7]-$total;
							$chrB=$explode2[$i];
							$entryB=1;
						}
					}
					else
					{
						if($explode[7] <= $hash_length{$explode2[$i]})
						{
							$coorB=$explode[7];
							$chrB=$explode2[$i];
							$entryB=1;
						}
					}
				}
				$total=$total+$hash_length{$explode2[$i]};
			}		
		}
	}
	$line2=$line;

	if($explode[2] ne "*" && $explode[3]!=0)
	{
		if($explode[2] ne "*")
		{
			$line2=~s/$explode[2]/$chrA/;
		}
		$line2=~s/$explode[3]/$coorA/;
	}

	if($chrA eq $chrB)
	{
		$line2=~s/$explode[6]/"="/;
		$line2=~s/$explode[7]/$coorB/;
	}
	else
	{
		if($explode[6] ne "*" && $explode[7] !=0)
		{
			$line2=~s/$explode[7]/$coorB/;
			if($explode[6] ne "*")
			{
				if($entry_equal==1)
				{
					$line2=~s/=/$chrB/;
				}
				else
				{
					$line2=~s/$explode[6]/$chrB/;
				}
			}
		}
	}

	@explode2=split("\t",$line2);

	if($explode2[2] ne "*" && $explode2[2] ne "=")
	{
		if($explode2[6] ne "*" && $explode2[6] ne "=")
		{
			if($explode2[2] ne $explode2[6])
			{
				$line2=~s/\t$explode[8]\t/\t0\t/;
			}
		}

	}

	print "$line2\n";

}
