#!/usr/bin/perl

sub usage 
{
    print "\nPARA QUE SIRVE: \n COMO SE USA:\n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(FASTA,"<",@ARGV[0]);
while  (my $line=<FASTA>)
{
        chomp($line);
        push (@fastafile,$line);
}


$counter=0;
for (my $i=$#fastafile; $i>=0; $i--)
{
	$first_character = substr ($fastafile[$i],0,1);
	if ($first_character eq ">")
	{
	
		if ($fastafile[$i] eq $fastafile[$i-2])
		{
			if ($counter==0)
			{
				print $fastafile[$i],"\n";
				print $fastafile[$i+1]; 
				$counter=$counter+1;
			}
			else
			{
				print $fastafile[$i+1];
			}
	}
		elsif (($fastafile[$i] ne $fastafile[$i-2]) && ($counter>0))
		{
			print $fastafile[$i+1],"\n";
			$counter=0;
			
		} 
	}


}
