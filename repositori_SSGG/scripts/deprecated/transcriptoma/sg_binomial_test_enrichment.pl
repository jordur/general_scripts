#!/usr/bin/perl

if(scalar(@ARGV)<3)
{
	die("\n\n\t\t*** I:Input: List Annotate Gene:Column 1 GeneID, Column 2 Descriptions,Column 3 Category,Column 4 Functional ID,Column 5, Other Descriptions.... The essential Columns are 1 and 4. Ejem Output of Blast2GO ***\n\n\t\t*** II:Reference File, Same format of I ***\n\n\t\t *** III:Threshold ***\n\n");
}
else
{
	%refID=();
	%description=();
	%numGeneRef=();

	%listID=();
	%Category=();
	%numGeneList=();
	%orfList=();
	%orfRef=();

	%Identify=();

	@explode=();

	$refCount=0;
	$listCount=0;
	$probability=0;

	$threshold=$ARGV[2];

	$refFile=$ARGV[1];

	open($BP_OUT,'>',"BP_Enrichment_Report.txt");
	open($MF_OUT,'>',"MF_Enrichment_Report.txt");
	open($CC_OUT,'>',"CC_Enrichment_Report.txt");
	open($REPORT_OUT,'>',"warning.err");

	open(FILE,"<",$refFile);
	while ($line = <FILE>)
	{
		chomp ($line);

		@explode=split("\t",$line);

		chomp($explode[0]);
		$exploe[0]=~s/\n//g;

		$numGeneRef{$explode[0]}=1;

		if(!$refID{$explode[3]})
		{
			$refID{$explode[3]}=1;
			$category{$explode[3]}=$explode[2];
			$description{$explode[3]}=$explode[4];

			$orfRef{$explode[3]}=$explode[0];
		}
		else
		{
			$refID{$explode[3]}=$refID{$explode[3]}+1;
			if($explode[0]=~/\w/)
			{
				$orfRef{$explode[3]}=$orfRef{$explode[3]}.",".$explode[0];
			}
		}

		$Identify{$explode[3]}=1;
	}

	foreach $keys (keys %numGeneRef)
	{
		$refCount=$refCount+1;
	}

	$listFile=$ARGV[0];
	open(FILE2,"<",$listFile);
	while ($line = <FILE2>)
	{
		chomp($line);
		
		@explode=split("\t",$line);

		if(!$listID{$explode[3]})
		{
			$listID{$explode[3]}=1;
			$category{$explode[3]}=$explode[2];
			$description{$explode[3]}=$explode[4];

			$orfList{$explode[3]}=$explode[0];

		}
		else
		{
			$listID{$explode[3]}=$listID{$explode[3]}+1;
			if($explode[0]=~/\w/)
			{
				$orfList{$explode[3]}=$orfList{$explode[3]}.",".$explode[0];
			}
		}
		$Identify{$explode[3]}=1;
		$numGeneList{$explode[0]}=1;
	}

	foreach $keys (keys %numGeneList)
	{
		 $listCount=$listCount+1;
	}

	print $BP_OUT "GO_IDar\tCategory\tDescription\tpValue\tTerm_count_in_reference\tTerm_count_in_sample\tORF_Reference\tORF_List\n";
	print $CC_OUT "GO_IDar\tCategory\tDescription\tpValue\tTerm_count_in_reference\tTerm_count_in_sample\tORF_Reference\tORF_List\n";
	print $MF_OUT "GO_IDar\tCategory\tDescription\tpValue\tTerm_count_in_reference\tTerm_count_in_sample\tORF_Reference\tORF_List\n";

	foreach $keys (keys %Identify)
	{
		if(!$refID{$keys})
		{
			$refID{$keys}=0;
		}
		if(!$listID{$keys})
		{
			$listID{$keys}=0;
		}
		$probability=$refID{$keys}/$refCount;
#		print "$keys\t ref $refID{$keys}\t list $listID{$keys}\t geneList $listCount\t background prob $probability\n";

		
		$system=`Rscript /share/apps/scripts/transcriptoma/binomial.R $listID{$keys} $listCount $probability`;

		@stringSystem=split(",",$system);
		for(my $kk=0;$kk<scalar(@stringSystem);$kk++)
		{
       			 if($stringSystem[$kk]=~/ p-value/)
       			 {
		                @subStringSystem=split("\n",$stringSystem[$kk]);
                		$subStringSystem[0]=~s/ //gi;
                		$subStringSystem[0]=~s/\n//gi;
                		$subStringSystem[0]=~s/p-value=//gi;
				$subStringSystem[0]=~s/p-value<//gi;

				chomp($subStringSystem[0]);
				chomp($orfRef{$keys});
				chomp($orfList{$keys});

				if($subStringSystem[0]<=$threshold && $refID{$keys}!=$listID{$keys})
				{

#print "$keys\t$category{$keys}\t$description{$keys}\t$subStringSystem[0]\t(ref $refID{$keys} list $listID{$keys} $probability)\n";
					if($category{$keys} eq "P")
					{
						print $BP_OUT "$keys\t$category{$keys}\t$description{$keys}\t$subStringSystem[0]\t$refID{$keys}\t$listID{$keys}\t$orfRef{$keys}\t$orfList{$keys}\n";
					}
					elsif($category{$keys} eq "F")
					{
                                                print $MF_OUT "$keys\t$category{$keys}\t$description{$keys}\t$subStringSystem[0]\t$refID{$keys}\t$listID{$keys}\t$orfRef{$keys}\t$orfList{$keys}\n";
					}
					elsif($category{$keys} eq "C")
					{
						print $CC_OUT "$keys\t$category{$keys}\t$description{$keys}\t$subStringSystem[0]\t$refID{$keys}\t$listID{$keys}\t$orfRef{$keys}\t$orfList{$keys}\n";
					}
					else
					{
						print $REPORT_OUT "$keys\t$category{$keys}\t$description{$keys}\t$subStringSystem[0]\t$orfRef{$keys}\t$orfList\t$refID{$keys}\t$$listID{$keys}\n";
					}
				}
			}
		}
	}


}
