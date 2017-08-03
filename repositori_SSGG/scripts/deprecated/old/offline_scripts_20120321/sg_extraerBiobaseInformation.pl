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
my $fistLine=0;

my $identification="";




open(FILE1,"<","$ARGV[1]");
@f=<FILE1>;

open(FILE,"<",$ARGV[0]);
while($line=<FILE>)
{
	chomp($line);
	
	if($fistLine==0)
	{
		print "$line\tDisease\tHGVS\tPMID\n";
		$fistLine=1;
		
	}
	else
	{
		$SNVs="NO";
		$entry=0;

		@annotationExplode=split("\t",$line);
	
		if($annotationExplode[6] == 0) #IS SNVs
        	{
      	 	        $position=$annotationExplode[2];
       		        $SNVs="YES";
        	}
		else
		{
			$positionA=$annotationExplode[2];
			$positionB=$annotationExplode[2]+$annotationExplode[6];
			$SNVs="NO";
		}


		for($k=$count;$k<scalar(@f);$k++)
		{
			chomp($f[$k]);
	
			@explode1=split("\t",$f[$k]);
			$explode1[4]=~s/chr[0-9]{1,2}://g;

#print "$explode1[4]\n";

			$explode1[4]=~s/chr[A-Z]://g;
	
			@explode2=split("-",$explode1[4]);
		
			if($SNVs eq "YES")
			{
				if($position>=$explode2[0] && $position <=$explode2[1])	
				{

#print "$position\t$explode2[0]\t$explode2[1]\n";

					$entry=1;
					$count=$k;	
	
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
			else
			{
				if($explode2[0] <=$positionA && $explode2[1] >=$positionB)
				{
       		                           $entry=1;
       	       	         	           $count=$k;
                     		           if($explode[11]!~/A-Z]/)
                              		   {
                                     		    print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";

                              		   }
                               	 	   else
				       	{
 	                              		 print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
					}
				}
				elsif($explode2[0]>=$positionA && $explode2[0] <=$positionB && $explode2[1]>=$positionA && $explode2[1] >=$positionB)
				{
					$entry=1;
					$count=$k;
               		                if($explode[11]!~/A-Z]/)
               	        	        {
                	               	          print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";
	
        	                        }
                        	        else
					{
						print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
					}
				}
				elsif($explode2[0]<=$positionA && $explode2[0]<=$positionB && $explode2[1] >=$positionA && $explode2[1]<=$positionB)
				{
                        	         $entry=1;
                               		 $count=$k;
                               		 if($explode[11]!~/A-Z]/)
                                	 {
                                         	 print "$line\t$explode1[7]\t$explode1[9]\thttp://www.ncbi.nlm.nih.gov/pubmed/$explode1[11]\n";

                                	 }
                              	         else
 				       	{

                               			 print "$line\t$explode1[7]\t$explode1[9]\t$explode1[11]\n";
			 		}
				}			
				elsif($explode2[0]>=$positionA && $explode2[1] <=$positionB)
				{
        	                        $entry=1;
                	                $count=$k;
	
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
		if($entry==0)
		{
			print "$line\t-\t-\t-\n";
		}
	}
}

