#!/usr/bin/perl -w

use strict;

if($ARGV[0] eq "" || $ARGV[1] eq "")
{
	print "\n\nINPUTS:\n\t\tI:Input File (output blastx.parseado ejm /data/results/Solid0065/BF5_Trocar/GATC_454/mira/ensamblajes)\n\t\tII:% of fragment (ejm 50, if the  alignmented segment is >= of % isn't annotated\n"; 
}
else
{
	my $input=$ARGV[0];
	my $percent=$ARGV[1];
	
	my @explode=();
	my @groupName=();
	my @groupStart=();
	my @groupEnd=();
	my @groupLine=();
	my @prohibitStart=();
	my @prohibitEnd=();


	my $count=0;
	my $fraction_conserved=-1;
	my $startQuery=0;
	my $endQuery=0;
	my $nameID=0;
	my $subcount=0;
	my $ini=0;
	my $end=0;
	my $size=0;
	my $size2=0;
	my $sizeAllow=0;
	my $mainEntry=0;
	my $entry=1;
	my $countprohibit=0;
	my $diferenceLeft=0;
	my $sum=0;
	
	my $line="";
	my $group="";


	open(FILE,'<',$input);
	while($line=<FILE>)
	{
		chomp($line);

		@explode=split("\t",$line);

		if($count==0)
		{
			for(my $i=0;$i<scalar(@explode);$i++)
			{
				if($explode[$i] eq "%_conserved")
				{
#					print "$i\t$explode[$i]\n";	#TESTED OK....
					$fraction_conserved=$i;
				}
				if($explode[$i] eq "Query_start")
				{
#					print "$i\t$explode[$i]\n";	#TESTED OK....
					$startQuery=$i;
				}
				if($explode[$i] eq "Query_end")
				{
#					print "$i\t$explode[$i]\n";	#TESTED OK....
					$endQuery=$i;
				}
				if($explode[$i] eq "Query_name")	#TESTED OK....
				{
#					print "$i\t$explode[$i]\n";
					$nameID=$i;
				}
			}
		
#print "$fraction_conserved\t$startQuery\t$endQuery\t$nameID\n";


			if($fraction_conserved==-1)
			{
				if($ARGV[3] eq "")
				{
					print "\n\n\t\t** The column fraction _conserved is missing!!, please write the correct header columns **\n\n\n";
					die;
				}
				else
				{
	                	       for(my $i=0;$i<scalar(@explode);$i++)
                        		{
					 	if($explode[$i] eq $ARGV[3])
                               		 	{
							$fraction_conserved=$i;
						}
					}
				}
			}
		}
		else
		{
			if($group ne $explode[$nameID])
			{
				if(@groupName)
				{
					@prohibitStart=();
					@prohibitEnd=();

					$countprohibit=0;

					#print "$groupName[0]\t$groupStart[0]\t$groupEnd[0]\n";
					print "$groupLine[0]\n";	
					for(my $j=1;$j<scalar(@groupName);$j++) #BECAUSE THE FIRST ALIGNAMENT IS TARGET
					{			
						$ini=$groupStart[$j];
       		      				$end=$groupEnd[$j];
       					        $size=$end-$ini;
             					$sizeAllow=($size*$percent)/100;
   
						$mainEntry=0;
 
             					if($ini>=$groupStart[0])
             					{
                     					if($ini<$groupEnd[0])
                     					{
                            					if($end<=$groupEnd[0])
                             					{
                                     					#NOTHING
                             					}
                             					else
                             					{
                                     					$size2=$size-($end-$groupEnd[0]);
                                     					if($size2<$sizeAllow)
                                     					{
										$mainEntry=1;
                                             					#print "\t$groupName[$j]\t$groupStart[$j]\t$groupEnd[$j]\n";
										#print "$groupLine[$j]\n";
                                     					}
                             					}
                     					}
                     					else
                    	 				{
								$mainEntry=1;
                             					#print "\t$groupName[$j]\t$groupStart[$j]\t$groupEnd[$j]\n";
								#print "$groupLine[$j]\n";
                     					}
             					}
             					else
             					{
               		      				if($end<$groupStart[0])
               		      				{
								$mainEntry=1;
               	        	      				 #print "\t$groupName[$j]\t$groupStart[$j]\t$groupEnd[$j]\n";
								 #print "$groupLine[$j]\n";
                     					}
                     					if($end >=$groupStart[0] && $end<=$groupEnd[0])
                     					{
                             					$size2=$size-($groupStart[0]-$ini);
                             					if($size2<$sizeAllow)
                             					{
									$mainEntry=1;
                                      					#print "\t$groupName[$j]\t$groupStart[$j]\t$groupEnd[$j]\n";
									#print "$groupLine[$j]\n";
                             					}
                     					}
                     					if($end>$groupEnd[0])
                     					{
                             					$size2=$size;
                             					if($size2<$sizeAllow)
                             					{
									$mainEntry=1;
                                     					#print "\t$groupName[$j]\t$groupStart[$j]\t$groupEnd[$j]\n";
									#print "$groupLine[$j]\n";
								}
                             				}
                     				}
						if($mainEntry==1)
						{
							$entry=1;
							$diferenceLeft=0;
							$sum=0;			
				
							if(@prohibitStart)
							{
								for(my $k=0;$k<scalar(@prohibitStart);$k++)
								{

#print "$ini==$prohibitStart[$k]\t$end==$prohibitStart[$k]\n";


									if($ini==$prohibitStart[$k])
									{			
										if($end==$prohibitStart[$k])
										{

#print "$ini==$prohibitStart[$k]\t$end==$prohibitStart[$k]\n";
											$entry=0;
										}
									}

									if($ini>=$prohibitStart[$k] && $end<=$prohibitStart[$k])
									{
										$entry=0;
									}
									else
									{
										if($ini<$prohibitStart[$k])
										{
											$diferenceLeft=$prohibitStart[$k]-$ini;
											if($end>$prohibitStart[$k])
											{
												if($end<=$prohibitEnd[$k])
												{
													if($diferenceLeft <= $sizeAllow)
													{
														$entry=0;
													}
												}
												else
												{
													$sum=($prohibitStart[$k]-$ini)+($end-$prohibitEnd[$k]);
													if($sum<=$sizeAllow)
													{
														$entry=0;	
													}
												}
											}
										}
										if($ini>$prohibitStart[$k])
										{
											$sum=$prohibitEnd[$k]-$ini;
											if($ini<$prohibitEnd[$k])
											{
												if($sum>=$sizeAllow)
												{
													$entry=0;
												}
											}
										}
									}
								}
								if($entry==1)
								{
									$prohibitStart[$countprohibit]=$ini;
									$prohibitEnd[$countprohibit]=$end;
									$countprohibit=$countprohibit+1;
								#	print "$groupLine[$j]\n";
								}
								else
								{
									 $prohibitStart[$countprohibit]=$ini;
                                                                         $prohibitEnd[$countprohibit]=$end;
                                                                         $countprohibit=$countprohibit+1;

								}
							}
							else
							{
								$prohibitStart[$countprohibit]=$ini;
								$prohibitEnd[$countprohibit]=$end;
								$countprohibit=$countprohibit+1;
								print "$groupLine[$j]\n";
							}
						}	
	           			}

					@groupName=();
					@groupStart=();
					@groupEnd=();
					@groupLine=();
					$subcount=0;
					
					if($explode[$fraction_conserved]>=65)
					{
						$groupName[$subcount]="$explode[$nameID]";
						$groupStart[$subcount]=$explode[$startQuery];
						$groupEnd[$subcount]=$explode[$endQuery];
						$groupLine[$subcount]=$line;
						$subcount++;
					}
				}
				else
				{
					if($explode[$fraction_conserved]>=65)
					{
						$groupName[$subcount]="$explode[$nameID]";
						$groupStart[$subcount]=$explode[$startQuery];
						$groupEnd[$subcount]=$explode[$endQuery];
						$groupLine[$subcount]=$line;
						$subcount++;
					}
				}
				$group=$explode[$nameID];
			}
			else
			{
				$groupName[$subcount]="$explode[$nameID]";
				$groupStart[$subcount]=$explode[$startQuery];
				$groupEnd[$subcount]=$explode[$endQuery];
				$groupLine[$subcount]=$line;
				$subcount++;
			}
		}
		$count=$count+1;
	}
if(@groupName)
{
		@prohibitStart=();
		@prohibitEnd=();


		$mainEntry=0;

		#print "FINALL\t$groupName[0]\t$groupStart[0]\t$groupEnd[0]\n";
		print "$groupLine[0]\n";
		for(my $j=1;$j<scalar(@groupName);$j++) #BECAUSE THE FIRST ALIGNAMENT IS TARGET
       	 	{
			$mainEntry=0;

     		   $ini=$groupStart[$j];
                   $end=$groupEnd[$j];
                   $size=$end-$ini;
                   $sizeAllow=($size* $percent)/100;
                   if($ini>=$groupStart[0])
                   {
      		             if($ini<$groupEnd[0])
                             {
                	             if($end<=$groupEnd[0])
                                     {
                        	             #NOTHING
                                     }
                                     else
                                     {
                                             $size2=$size-($end-$groupEnd[0]);
                                             if($size2<$sizeAllow)
                                             {
                                 	             #print "\t$groupName[$j]\t$groupStart[$j]\t$groupEnd[$j]\n";
						     #print "$groupLine[$j]\n";
						     $mainEntry=1;
                                             }
                                     }
                              }
                              else
                              {
                      	              #print "\t$groupName[$j]\t$groupStart[$j]\t$groupEnd[$j]\n";
				      #print "$groupLine[$j]\n";
				      $mainEntry=1;
                              }
                   }
                   else
                   {
                              if($end<$groupStart[0])
                              {
                                      #print "\t$groupName[$j]\t$groupStart[$j]\t$groupEnd[$j]\n";
				      #print "$groupLine[$j]\n";
				      $mainEntry=1;
                              }
                              if($end >=$groupStart[0] && $end<=$groupEnd[0])
                              {
                                      $size2=$size-($groupStart[0]-$ini);
                                      if($size2<$sizeAllow)
                                      {
              	                                     #print "\t$groupName[$j]\t$groupStart[$j]\t$groupEnd[$j]\n";
						     #print "$groupLine[$j]\n";
						     $mainEntry=1;
                                      }
                              }
                              if($end>$groupEnd[0])
                              {
                                      $size2=$size;
                                      if($size2<$sizeAllow)
                                      {
                                                    # print "\t$groupName[$j]\t$groupStart[$j]\t$groupEnd[$j]\n";
						    # print "$groupLine[$j]\n";
						    $mainEntry=1;
                                      }
                              }
		     }
                                                 if($mainEntry==1)
                                                 {
                                                         $entry=1;
                                                         $diferenceLeft=0;
                                                         $sum=0;

                                                         if(@prohibitStart)
                                                         {
                                                                 for(my $k=0;$k<scalar(@prohibitStart);$k++)
                                                                 {
                                                                         if($ini>=$prohibitStart[$k] && $end<=$prohibitStart[$k])
                                                                         {
                                                                                 $entry=0;
                                                                         }
                                                                         else
                                                                         {
                                                                                 if($ini<$prohibitStart[$k])
                                                                                 {
                                                                                         $diferenceLeft=$prohibitStart[$k]-$ini;
                                                                                         if($end>$prohibitStart[$k])
                                                                                         {
                                                                                                 if($end<=$prohibitEnd[$k])
                                                                                                 {
                                                                                                         if($diferenceLeft <= $sizeAllow)
                                                                                                         {
                                                                                                                 $entry=0;
                                                                                                         }
                                                                                                 }
                                                                                                 else
                                                                                                 {
                                                                                                         $sum=($prohibitStart[$k]-$ini)+($end-$prohibitEnd    [$k]);
                                                                                                         if($sum<=$sizeAllow)
                                                                                                         {
                                                                                                                 $entry=0;
                                                                                                         }
                                                                                                 }
                                                                                         }
                                                                                 }
                                                                                 if($ini>$prohibitStart[$k])
                                                                                 {
                                                                                         $sum=$prohibitEnd[$k]-$ini;
                                                                                         if($ini<$prohibitEnd[$k])
                                                                                         {
                                                                                                 if($sum>=$sizeAllow)
                                                                                                 {
                                                                                                         $entry=0;
                                                                                                 }
                                                                                         }
                                                                                 }
                                                                         }
                                                                 }
                                                                 if($entry==1)
                                                                 {
                                                                         $prohibitStart[$countprohibit]=$ini;
                                                                         $prohibitEnd[$countprohibit]=$end;
                                                                         $countprohibit=$countprohibit+1;
                                                                         print "$groupLine[$j]\n";
                                                                 }
        						}
							else
                                                        {
                                                                 $prohibitStart[$countprohibit]=$ini;
                                                                 $prohibitEnd[$countprohibit]=$end;
                                                                 $countprohibit=$countprohibit+1;
                                                                 print "$groupLine[$j]\n";
                                                         }

						}
                    }
       }           
}

